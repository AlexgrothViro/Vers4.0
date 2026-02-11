#!/usr/bin/env python3
import argparse
import json
import os
import re
import shutil
import threading
import time
import uuid
import zipfile
from datetime import datetime
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from subprocess import PIPE, STDOUT, Popen, run
from urllib.parse import parse_qs, urlparse

REPO_ROOT = Path(__file__).resolve().parents[1]
DASHBOARD_DIR = REPO_ROOT / "dashboard"
LOG_DIR = REPO_ROOT / "logs"
RUNS_DIR = REPO_ROOT / "results" / "runs"
TARGETS_FILE = REPO_ROOT / "config" / "targets.json"

MAX_OUTPUT_LINES = 4000

jobs = {}
jobs_lock = threading.Lock()


def sanitize_token(text):
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", (text or "unknown").strip())
    return value[:80] or "unknown"


def iso_now(epoch=None):
    dt = datetime.fromtimestamp(epoch or time.time())
    return dt.isoformat(timespec="seconds")


def read_tail(path, lines=30):
    try:
        content = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except FileNotFoundError:
        return "(log não encontrado)"
    return "\n".join(content[-lines:])


def list_targets():
    if not TARGETS_FILE.exists():
        return []
    try:
        data = json.loads(TARGETS_FILE.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return []
    return data if isinstance(data, list) else []

def list_db_profiles():
    try:
        completed = run(
            ["bash", "-lc", "scripts/13_db_manager.sh list --json"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        if completed.returncode != 0:
            return []
        data = json.loads(completed.stdout or "[]")
        return data if isinstance(data, list) else []
    except Exception:
        return []


def validate_fastq_name(name):
    return name.endswith(".fastq") or name.endswith(".fastq.gz")


def normalize_sample_from_filename(filename):
    name = filename
    for suffix in ["_R1.fastq.gz", "_R2.fastq.gz", "_R1.fastq", "_R2.fastq"]:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return ""


def list_samples():
    raw_dir = REPO_ROOT / "data" / "raw"
    if not raw_dir.exists():
        return []
    sample_map = {}
    for item in raw_dir.iterdir():
        if not item.is_file() or not validate_fastq_name(item.name):
            continue
        sample = normalize_sample_from_filename(item.name)
        if not sample:
            continue
        side = "R1" if "_R1." in item.name else "R2"
        sample_map.setdefault(sample, set()).add(side)
    return sorted(sample for sample, sides in sample_map.items() if {"R1", "R2"}.issubset(sides))


def tool_versions():
    tools = {
        "python": ["python3", "--version"],
        "blastn": ["blastn", "-version"],
        "bowtie2": ["bowtie2", "--version"],
        "velveth": ["velveth", "--help"],
        "spades": ["spades.py", "--version"],
        "make": ["make", "--version"],
    }
    versions = {}
    for name, cmd in tools.items():
        try:
            completed = run(cmd, cwd=REPO_ROOT, capture_output=True, text=True, check=False)
            first = (completed.stdout or completed.stderr or "").splitlines()
            versions[name] = first[0].strip() if first else "available"
        except FileNotFoundError:
            versions[name] = "not found"
    return versions


def clamp_output(lines):
    if len(lines) <= MAX_OUTPUT_LINES:
        return lines
    return lines[-MAX_OUTPUT_LINES:]


def find_blast_path(sample):
    blast_dir = REPO_ROOT / "results" / "blast"
    if not sample or not blast_dir.exists():
        return None
    preferred = sorted(blast_dir.glob(f"{sample}*_vs_db.tsv"), key=lambda p: p.stat().st_mtime, reverse=True)
    return preferred[0] if preferred else None


def find_report_path(sample):
    report = REPO_ROOT / "results" / "reports" / f"{sample}_summary.md"
    return report if sample and report.exists() else None


def snapshot_run_artifacts(metadata):
    if metadata.get("action") != "pipeline":
        return metadata
    sample = metadata.get("sample")
    if not sample:
        return metadata

    ts = datetime.fromtimestamp(metadata.get("end_epoch", time.time())).strftime("%Y%m%d_%H%M%S")
    run_dir = RUNS_DIR / f"{ts}_{sanitize_token(sample)}"
    run_dir.mkdir(parents=True, exist_ok=True)

    paths = metadata.setdefault("paths", {})
    log_path = Path(paths.get("log", "")) if paths.get("log") else None
    if log_path and log_path.exists():
        target = run_dir / log_path.name
        shutil.copy2(log_path, target)
        paths["run_log"] = str(target.relative_to(REPO_ROOT))

    blast_path = find_blast_path(sample)
    if blast_path and blast_path.exists():
        target = run_dir / blast_path.name
        shutil.copy2(blast_path, target)
        paths["run_blast"] = str(target.relative_to(REPO_ROOT))

    report_path = find_report_path(sample)
    if report_path and report_path.exists():
        target = run_dir / report_path.name
        shutil.copy2(report_path, target)
        paths["run_report"] = str(target.relative_to(REPO_ROOT))

    metadata["run_dir"] = run_dir.name
    run_json = run_dir / "run.json"
    run_json.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    paths["run_json"] = str(run_json.relative_to(REPO_ROOT))
    return metadata


def parse_pipeline_details(params):
    assembler = (params.get("assembler") or os.environ.get("ASSEMBLER") or "velvet").lower()
    kmer = str(params.get("kmer") or "31")
    sample = params.get("sample")
    threads = str(os.environ.get("THREADS", "4"))
    return sample, assembler, kmer, threads


def build_command(action, params):
    if action == "check_env":
        return ["make", "test-env"], {}
    if action == "demo":
        return ["make", "demo"], {}
    if action == "import_sample":
        sample = params.get("sample")
        r1 = params.get("r1")
        r2 = params.get("r2")
        copy = params.get("copy") is True
        if not sample or not r1 or not r2:
            raise ValueError("Campos obrigatórios ausentes: sample, r1, r2")
        cmd = ["bash", "scripts/00_import_sample.sh", "--sample", sample, "--r1", r1, "--r2", r2]
        if copy:
            cmd.append("--copy")
        return cmd, {}
    if action == "build_db":
        target = params.get("target")
        query = (params.get("query") or "").strip()
        taxid = (params.get("taxid") or "").strip()
        if not target:
            raise ValueError("Campo obrigatório ausente: target")
        cmd = ["bash", "scripts/10_build_viral_db.sh", "--target", target]
        if query:
            cmd += ["--query", query]
        elif taxid:
            cmd += ["--taxid", taxid]
        return cmd, {}
    if action == "pipeline":
        sample = params.get("sample")
    if not sample:
        raise ValueError("Campo obrigatório ausente: sample")

    assembler = (params.get("assembler") or "velvet").lower()
    kmer = params.get("kmer") or "31"

    # >>> NOVO: DB selecionável via UX
    db = (params.get("db") or "").strip()

    cmd = ["bash", "scripts/20_run_pipeline.sh", "--sample", sample, "--kmer", str(kmer)]
    env = {}

    if assembler in {"velvet", "spades"}:
        env["ASSEMBLER"] = assembler

    if db:
        env["DB"] = db   # usado pelo 13_db_manager.sh e/ou pipeline
        # opcional: também passar por argumento se teu 20_run_pipeline.sh suportar
        # cmd += ["--db", db]

    return cmd, env

    raise ValueError("Ação inválida")


def run_job(job_id, action, params):
    env = os.environ.copy()
    try:
        cmd, extra_env = build_command(action, params)
        env.update(extra_env)
    except ValueError as exc:
        with jobs_lock:
            jobs[job_id].update({"status": "error", "output": [f"[ERRO] {exc}\n"], "returncode": 1})
        return

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    now = datetime.now().strftime("%Y%m%d_%H%M%S")
    sample = params.get("sample") or params.get("target") or action
    log_name = f"ux_dashboard_{now}_{sanitize_token(action)}_{sanitize_token(sample)}.log"
    log_path = LOG_DIR / log_name

    start_epoch = time.time()
    pipeline_sample, assembler, kmer, threads = parse_pipeline_details(params)

    with jobs_lock:
        jobs[job_id].update({"status": "running", "command": " ".join(cmd), "log_path": str(log_path.relative_to(REPO_ROOT))})

    process = Popen(cmd, cwd=REPO_ROOT, env=env, stdout=PIPE, stderr=STDOUT, text=True, bufsize=1)

    output_lines = []
    with log_path.open("w", encoding="utf-8") as log_file:
        # prepara DB automaticamente quando existir DB no env
        if env.get("DB"):
            log_file.write(f"[INFO] Preparando DB via 13_db_manager.sh (DB={env['DB']})...\n")
            log_file.flush()

            completed = run(
                ["bash", "scripts/13_db_manager.sh", "setup"],
                cwd=REPO_ROOT,
                env=env,
                stdout=PIPE,
                stderr=STDOUT,
                text=True,
                check=False,
            )
            if completed.stdout:
                log_file.write(completed.stdout)
                log_file.flush()
                output_lines.extend(completed.stdout.splitlines(True))
                output_lines = clamp_output(output_lines)
                with jobs_lock:
                    jobs[job_id]["output"] = output_lines

        if process.stdout:
            for line in process.stdout:
                log_file.write(line)
                log_file.flush()
                output_lines.append(line)
                output_lines = clamp_output(output_lines)
                with jobs_lock:
                    jobs[job_id]["output"] = output_lines
    returncode = process.wait()
    end_epoch = time.time()

    metadata = {
        "id": job_id,
        "action": action,
        "sample": pipeline_sample or params.get("sample"),
        "assembler": assembler if action == "pipeline" else None,
        "kmer": kmer if action == "pipeline" else None,
        "threads": threads if action == "pipeline" else None,
        "command": " ".join(cmd),
        "start": iso_now(start_epoch),
        "end": iso_now(end_epoch),
        "start_epoch": start_epoch,
        "end_epoch": end_epoch,
        "exit_code": returncode,
        "params": params,
        "versions": tool_versions() if action == "pipeline" else {},
        "paths": {
            "log": str(log_path.relative_to(REPO_ROOT)),
            "report": str(find_report_path(pipeline_sample).relative_to(REPO_ROOT)) if find_report_path(pipeline_sample) else None,
            "blast": str(find_blast_path(pipeline_sample).relative_to(REPO_ROOT)) if find_blast_path(pipeline_sample) else None,
        },
    }
    metadata = snapshot_run_artifacts(metadata)

    with jobs_lock:
        jobs[job_id].update(
            {
                "returncode": returncode,
                "status": "done" if returncode == 0 else "error",
                "output": output_lines,
                "finished_at": end_epoch,
                "run": metadata,
                "tail": read_tail(log_path, lines=30) if returncode != 0 else "",
            }
        )


def list_run_history():
    runs = []
    for run_json in RUNS_DIR.glob("*/run.json"):
        try:
            data = json.loads(run_json.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            continue
        data["run_dir"] = run_json.parent.name
        runs.append(data)
    runs.sort(key=lambda item: item.get("end_epoch", 0), reverse=True)
    return runs


def resolve_history_file(run_dir, file_type):
    run_path = RUNS_DIR / run_dir
    run_json = run_path / "run.json"
    if not run_json.exists():
        return None
    data = json.loads(run_json.read_text(encoding="utf-8"))
    key_map = {"report": "run_report", "log": "run_log", "blast": "run_blast"}
    relpath = data.get("paths", {}).get(key_map.get(file_type, ""))
    if not relpath:
        return None
    target = REPO_ROOT / relpath
    if not target.exists() or run_path not in target.parents:
        return None
    return target


def json_response(handler, status, payload):
    body = json.dumps(payload).encode("utf-8")
    handler.send_response(status)
    handler.send_header("Content-Type", "application/json")
    handler.send_header("Content-Length", str(len(body)))
    handler.end_headers()
    handler.wfile.write(body)


def text_response(handler, status, text):
    body = text.encode("utf-8")
    handler.send_response(status)
    handler.send_header("Content-Type", "text/plain; charset=utf-8")
    handler.send_header("Content-Length", str(len(body)))
    handler.end_headers()
    handler.wfile.write(body)


def serve_file(handler, filepath, content_type):
    if not filepath.exists():
        handler.send_error(HTTPStatus.NOT_FOUND, "Arquivo não encontrado")
        return
    body = filepath.read_bytes()
    handler.send_response(HTTPStatus.OK)
    handler.send_header("Content-Type", content_type)
    handler.send_header("Content-Length", str(len(body)))
    handler.end_headers()
    handler.wfile.write(body)


def import_uploaded_files(sample, r1_name, r1_data, r2_name, r2_data):
    if not validate_fastq_name(r1_name) or not validate_fastq_name(r2_name):
        raise ValueError("Extensão inválida (use .fastq ou .fastq.gz)")
    if not r1_data:
        raise ValueError("arquivo vazio: R1")
    if not r2_data:
        raise ValueError("arquivo vazio: R2")

    raw_dir = REPO_ROOT / "data" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    out_r1 = raw_dir / f"{sample}_R1.fastq.gz"
    out_r2 = raw_dir / f"{sample}_R2.fastq.gz"

    import gzip

    if r1_name.endswith(".gz"):
        out_r1.write_bytes(r1_data)
    else:
        with gzip.open(out_r1, "wb") as fh:
            fh.write(r1_data)
    if r2_name.endswith(".gz"):
        out_r2.write_bytes(r2_data)
    else:
        with gzip.open(out_r2, "wb") as fh:
            fh.write(r2_data)


def import_zip_sample(sample, zip_bytes):
    if not zip_bytes:
        raise ValueError("arquivo vazio")

    with zipfile.ZipFile(Path(REPO_ROOT / "tmp_upload.zip"), "w") as _:
        pass
    temp_zip = REPO_ROOT / f"tmp_upload_{uuid.uuid4().hex}.zip"
    temp_zip.write_bytes(zip_bytes)

    try:
        with zipfile.ZipFile(temp_zip) as zf:
            names = [n for n in zf.namelist() if not n.endswith("/")]
            fastqs = [n for n in names if validate_fastq_name(n.lower())]
            r1_candidates = [n for n in fastqs if "r1" in n.lower()]
            r2_candidates = [n for n in fastqs if "r2" in n.lower()]
            if not r1_candidates:
                raise ValueError("R1 não encontrado no .zip")
            if not r2_candidates:
                raise ValueError("R2 não encontrado no .zip")

            r1_name = r1_candidates[0]
            r2_name = r2_candidates[0]
            import_uploaded_files(sample, Path(r1_name).name, zf.read(r1_name), Path(r2_name).name, zf.read(r2_name))
    finally:
        temp_zip.unlink(missing_ok=True)


class DashboardHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        parsed = urlparse(self.path)
        if parsed.path == "/":
            return serve_file(self, DASHBOARD_DIR / "index.html", "text/html; charset=utf-8")
        if parsed.path == "/api/dbs":
            return json_response(self, HTTPStatus.OK, {"dbs": list_db_profiles()})
        if parsed.path == "/styles.css":
            return serve_file(self, DASHBOARD_DIR / "styles.css", "text/css; charset=utf-8")
        if parsed.path == "/app.js":
            return serve_file(self, DASHBOARD_DIR / "app.js", "application/javascript; charset=utf-8")
        if parsed.path == "/api/samples":
            return json_response(self, HTTPStatus.OK, {"samples": list_samples()})
        if parsed.path == "/api/targets":
            return json_response(self, HTTPStatus.OK, {"targets": list_targets()})
        if parsed.path == "/api/history":
            return json_response(self, HTTPStatus.OK, {"runs": list_run_history()})
        ...

        if parsed.path == "/api/targets":
            return json_response(self, HTTPStatus.OK, {"targets": list_targets()})
        if parsed.path == "/api/history":
            return json_response(self, HTTPStatus.OK, {"runs": list_run_history()})
        if parsed.path == "/api/history/file":
            query = parse_qs(parsed.query)
            target = resolve_history_file((query.get("run") or [""])[0], (query.get("type") or [""])[0])
            if not target:
                return self.send_error(HTTPStatus.NOT_FOUND, "Arquivo do histórico não encontrado")
            ctype = "text/markdown; charset=utf-8" if target.suffix == ".md" else "text/plain; charset=utf-8"
            return serve_file(self, target, ctype)
        if parsed.path.startswith("/api/job/"):
            job_id = parsed.path.rsplit("/", 1)[-1]
            with jobs_lock:
                job = jobs.get(job_id)
            if not job:
                return json_response(self, HTTPStatus.NOT_FOUND, {"error": "Job não encontrado"})
            return json_response(self, HTTPStatus.OK, {"id": job_id, "status": job.get("status"), "output": "".join(job.get("output", [])), "returncode": job.get("returncode"), "command": job.get("command"), "log_path": job.get("log_path"), "run": job.get("run"), "tail": job.get("tail", "")})
        return self.send_error(HTTPStatus.NOT_FOUND, "Rota inválida")

    def do_POST(self):
        parsed = urlparse(self.path)
        if parsed.path == "/api/import-upload":
            return self.handle_upload_import()

        content_length = int(self.headers.get("Content-Length", 0))
        if content_length <= 0:
            return text_response(self, HTTPStatus.BAD_REQUEST, "Body vazio")
        raw_body = self.rfile.read(content_length)
        try:
            payload = json.loads(raw_body)
        except json.JSONDecodeError:
            return text_response(self, HTTPStatus.BAD_REQUEST, "JSON inválido")

        if parsed.path == "/api/history/rerun":
            run_json = RUNS_DIR / str(payload.get("run_dir")) / "run.json"
            if not run_json.exists():
                return text_response(self, HTTPStatus.NOT_FOUND, "run.json não encontrado")
            data = json.loads(run_json.read_text(encoding="utf-8"))
            return self.start_job(data.get("action", "pipeline"), data.get("params") or {})

        if parsed.path != "/api/run":
            return self.send_error(HTTPStatus.NOT_FOUND, "Rota inválida")
        return self.start_job(payload.get("action"), payload.get("params") or {})

    def handle_upload_import(self):
        import cgi

        form = cgi.FieldStorage(fp=self.rfile, headers=self.headers, environ={"REQUEST_METHOD": "POST", "CONTENT_TYPE": self.headers.get("Content-Type", "")})
        sample = (form.getvalue("sample") or "").strip()
        if not sample:
            return text_response(self, HTTPStatus.BAD_REQUEST, "Campo obrigatório: sample")
        try:
            zip_file = form["zipfile"] if "zipfile" in form else None
            r1_file = form["r1file"] if "r1file" in form else None
            r2_file = form["r2file"] if "r2file" in form else None

            if zip_file is not None and getattr(zip_file, "filename", ""):
                import_zip_sample(sample, zip_file.file.read())
            elif r1_file is not None and r2_file is not None and getattr(r1_file, "filename", "") and getattr(r2_file, "filename", ""):
                import_uploaded_files(sample, r1_file.filename, r1_file.file.read(), r2_file.filename, r2_file.file.read())
            else:
                return text_response(self, HTTPStatus.BAD_REQUEST, "Envie R1/R2 ou arquivo .zip")
        except (ValueError, zipfile.BadZipFile) as exc:
            return text_response(self, HTTPStatus.BAD_REQUEST, str(exc))
        return json_response(self, HTTPStatus.OK, {"message": f"Amostra importada: {sample}"})

    def start_job(self, action, params):
        if action not in {"check_env", "demo", "import_sample", "build_db", "pipeline"}:
            return text_response(self, HTTPStatus.BAD_REQUEST, "Ação inválida")
        job_id = uuid.uuid4().hex
        with jobs_lock:
            jobs[job_id] = {"status": "queued", "created_at": time.time(), "action": action, "output": [], "returncode": None}
        threading.Thread(target=run_job, args=(job_id, action, params), daemon=True).start()
        return json_response(self, HTTPStatus.OK, {"job_id": job_id})

    def log_message(self, format, *args):
        return


def main():
    parser = argparse.ArgumentParser(description="Painel de uso para o pipeline Picornavirus-quali2026")
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument("--port", default=8000, type=int)
    args = parser.parse_args()

    server = ThreadingHTTPServer((args.host, args.port), DashboardHandler)
    print("Painel ativo em:")
    print(f"  http://localhost:{args.port}")
    print(f"  http://{args.host}:{args.port}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nEncerrando painel...")


if __name__ == "__main__":
    main()
