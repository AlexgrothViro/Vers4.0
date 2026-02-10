#!/usr/bin/env python3
import argparse
import json
import os
import re
import shutil
import threading
import time
import uuid
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


def list_samples():
    raw_dir = REPO_ROOT / "data" / "raw"
    if not raw_dir.exists():
        return []
    samples = set()
    for item in raw_dir.glob("*_R1.fastq.gz"):
        samples.add(item.name.replace("_R1.fastq.gz", ""))
    return sorted(samples)


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
    action = metadata.get("action")
    if action != "pipeline":
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
        cmd = [
            "bash",
            "scripts/00_import_sample.sh",
            "--sample",
            sample,
            "--r1",
            r1,
            "--r2",
            r2,
        ]
        if copy:
            cmd.append("--copy")
        return cmd, {}
    if action == "pipeline":
        sample = params.get("sample")
        if not sample:
            raise ValueError("Campo obrigatório ausente: sample")
        assembler = (params.get("assembler") or "velvet").lower()
        kmer = params.get("kmer") or "31"
        cmd = [
            "bash",
            "scripts/20_run_pipeline.sh",
            "--sample",
            sample,
            "--kmer",
            str(kmer),
        ]
        env = {}
        if assembler in {"velvet", "spades"}:
            env["ASSEMBLER"] = assembler
        return cmd, env
    raise ValueError("Ação inválida")


def run_job(job_id, action, params):
    env = os.environ.copy()
    try:
        cmd, extra_env = build_command(action, params)
        env.update(extra_env)
    except ValueError as exc:
        with jobs_lock:
            job = jobs[job_id]
            job["status"] = "error"
            job["output"] = [f"[ERRO] {exc}\n"]
            job["returncode"] = 1
        return

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    now = datetime.now().strftime("%Y%m%d_%H%M%S")
    sample = params.get("sample") or action
    log_name = f"ux_dashboard_{now}_{sanitize_token(action)}_{sanitize_token(sample)}.log"
    log_path = LOG_DIR / log_name

    start_epoch = time.time()
    pipeline_sample, assembler, kmer, threads = parse_pipeline_details(params)

    with jobs_lock:
        jobs[job_id]["status"] = "running"
        jobs[job_id]["command"] = " ".join(cmd)
        jobs[job_id]["log_path"] = str(log_path.relative_to(REPO_ROOT))

    process = Popen(
        cmd,
        cwd=REPO_ROOT,
        env=env,
        stdout=PIPE,
        stderr=STDOUT,
        text=True,
        bufsize=1,
    )

    output_lines = []
    with log_path.open("w", encoding="utf-8") as log_file:
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
        jobs[job_id]["returncode"] = returncode
        jobs[job_id]["status"] = "done" if returncode == 0 else "error"
        jobs[job_id]["output"] = output_lines
        jobs[job_id]["finished_at"] = end_epoch
        jobs[job_id]["run"] = metadata
        jobs[job_id]["tail"] = read_tail(log_path, lines=30) if returncode != 0 else ""


def list_run_history():
    runs = []
    for run_json in RUNS_DIR.glob("*/run.json"):
        try:
            data = json.loads(run_json.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            continue
        run_dir = run_json.parent.name
        data["run_dir"] = run_dir
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


class DashboardHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        parsed = urlparse(self.path)
        if parsed.path == "/":
            return serve_file(self, DASHBOARD_DIR / "index.html", "text/html; charset=utf-8")
        if parsed.path == "/styles.css":
            return serve_file(self, DASHBOARD_DIR / "styles.css", "text/css; charset=utf-8")
        if parsed.path == "/app.js":
            return serve_file(self, DASHBOARD_DIR / "app.js", "application/javascript; charset=utf-8")
        if parsed.path == "/api/samples":
            return json_response(self, HTTPStatus.OK, {"samples": list_samples()})
        if parsed.path == "/api/history":
            return json_response(self, HTTPStatus.OK, {"runs": list_run_history()})
        if parsed.path == "/api/history/file":
            query = parse_qs(parsed.query)
            run_dir = (query.get("run") or [""])[0]
            file_type = (query.get("type") or [""])[0]
            target = resolve_history_file(run_dir, file_type)
            if not target:
                return self.send_error(HTTPStatus.NOT_FOUND, "Arquivo do histórico não encontrado")
            ctype = "text/plain; charset=utf-8"
            if target.suffix == ".md":
                ctype = "text/markdown; charset=utf-8"
            return serve_file(self, target, ctype)
        if parsed.path.startswith("/api/job/"):
            job_id = parsed.path.rsplit("/", 1)[-1]
            with jobs_lock:
                job = jobs.get(job_id)
            if not job:
                return json_response(self, HTTPStatus.NOT_FOUND, {"error": "Job não encontrado"})
            return json_response(
                self,
                HTTPStatus.OK,
                {
                    "id": job_id,
                    "status": job.get("status"),
                    "output": "".join(job.get("output", [])),
                    "returncode": job.get("returncode"),
                    "command": job.get("command"),
                    "log_path": job.get("log_path"),
                    "run": job.get("run"),
                    "tail": job.get("tail", ""),
                },
            )
        return self.send_error(HTTPStatus.NOT_FOUND, "Rota inválida")

    def do_POST(self):
        parsed = urlparse(self.path)
        if parsed.path == "/api/history/rerun":
            content_length = int(self.headers.get("Content-Length", 0))
            if content_length <= 0:
                return text_response(self, HTTPStatus.BAD_REQUEST, "Body vazio")
            raw_body = self.rfile.read(content_length)
            try:
                payload = json.loads(raw_body)
            except json.JSONDecodeError:
                return text_response(self, HTTPStatus.BAD_REQUEST, "JSON inválido")
            run_dir = payload.get("run_dir")
            run_json = RUNS_DIR / str(run_dir) / "run.json"
            if not run_json.exists():
                return text_response(self, HTTPStatus.NOT_FOUND, "run.json não encontrado")
            data = json.loads(run_json.read_text(encoding="utf-8"))
            action = data.get("action", "pipeline")
            params = data.get("params") or {}
            return self.start_job(action, params)

        if parsed.path != "/api/run":
            return self.send_error(HTTPStatus.NOT_FOUND, "Rota inválida")

        content_length = int(self.headers.get("Content-Length", 0))
        if content_length <= 0:
            return text_response(self, HTTPStatus.BAD_REQUEST, "Body vazio")
        raw_body = self.rfile.read(content_length)

        try:
            payload = json.loads(raw_body)
        except json.JSONDecodeError:
            return text_response(self, HTTPStatus.BAD_REQUEST, "JSON inválido")

        action = payload.get("action")
        params = payload.get("params") or {}

        return self.start_job(action, params)

    def start_job(self, action, params):
        if action not in {"check_env", "demo", "import_sample", "pipeline"}:
            return text_response(self, HTTPStatus.BAD_REQUEST, "Ação inválida")

        job_id = uuid.uuid4().hex
        with jobs_lock:
            jobs[job_id] = {
                "status": "queued",
                "created_at": time.time(),
                "action": action,
                "output": [],
                "returncode": None,
            }

        thread = threading.Thread(target=run_job, args=(job_id, action, params), daemon=True)
        thread.start()

        return json_response(self, HTTPStatus.OK, {"job_id": job_id})

    def log_message(self, format, *args):
        return


def main():
    parser = argparse.ArgumentParser(description="Painel de uso para o pipeline Picornavirus-quali2026")
    parser.add_argument("--host", default="0.0.0.0", help="Host para bind (padrão: 0.0.0.0)")
    parser.add_argument("--port", default=8000, type=int, help="Porta do painel (padrão: 8000)")
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
