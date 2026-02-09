#!/usr/bin/env python3
import argparse
import json
import os
import threading
import time
import uuid
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from subprocess import Popen, PIPE, STDOUT
from urllib.parse import urlparse

REPO_ROOT = Path(__file__).resolve().parents[1]
DASHBOARD_DIR = REPO_ROOT / "dashboard"
LOG_DIR = REPO_ROOT / "logs"

MAX_OUTPUT_LINES = 4000

jobs = {}
jobs_lock = threading.Lock()


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

    with jobs_lock:
        jobs[job_id]["status"] = "running"
        jobs[job_id]["command"] = " ".join(cmd)

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
    if process.stdout:
        for line in process.stdout:
            output_lines.append(line)
            output_lines = clamp_output(output_lines)
            with jobs_lock:
                jobs[job_id]["output"] = output_lines

    returncode = process.wait()

    with jobs_lock:
        jobs[job_id]["returncode"] = returncode
        jobs[job_id]["status"] = "done" if returncode == 0 else "error"
        jobs[job_id]["output"] = output_lines
        jobs[job_id]["finished_at"] = time.time()


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
                },
            )
        return self.send_error(HTTPStatus.NOT_FOUND, "Rota inválida")

    def do_POST(self):
        parsed = urlparse(self.path)
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
