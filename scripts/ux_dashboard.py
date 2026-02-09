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
from subprocess import Popen, STDOUT
from urllib.parse import parse_qs, urlparse

ROOT = Path(__file__).resolve().parents[1]
LOG_DIR = ROOT / "logs"
LOG_DIR.mkdir(exist_ok=True)

JOBS = {}
JOBS_LOCK = threading.Lock()

HTML_PAGE = """<!DOCTYPE html>
<html lang=\"pt-BR\">
<head>
  <meta charset=\"UTF-8\" />
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />
  <title>Painel Picornavirus-quali2026</title>
  <style>
    :root {
      color-scheme: light;
      --bg: #f7f7fb;
      --card: #ffffff;
      --text: #1d1d21;
      --muted: #5c5f6b;
      --accent: #2f6fed;
      --accent-dark: #1b4bbf;
      --border: #e0e3eb;
    }
    body {
      font-family: "Segoe UI", system-ui, -apple-system, sans-serif;
      background: var(--bg);
      color: var(--text);
      margin: 0;
      padding: 24px;
    }
    header {
      display: flex;
      flex-direction: column;
      gap: 8px;
      margin-bottom: 20px;
    }
    h1 {
      font-size: 24px;
      margin: 0;
    }
    .subtitle {
      color: var(--muted);
      font-size: 14px;
    }
    .grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
      gap: 16px;
    }
    .card {
      background: var(--card);
      border: 1px solid var(--border);
      border-radius: 12px;
      padding: 16px;
      box-shadow: 0 10px 20px rgba(0, 0, 0, 0.04);
      display: flex;
      flex-direction: column;
      gap: 12px;
    }
    .card h2 {
      font-size: 16px;
      margin: 0;
    }
    label {
      font-size: 13px;
      color: var(--muted);
      display: flex;
      flex-direction: column;
      gap: 6px;
    }
    input, select, textarea {
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 8px 10px;
      font-size: 14px;
      font-family: inherit;
    }
    textarea {
      min-height: 160px;
      resize: vertical;
      background: #111827;
      color: #f9fafb;
      font-family: "Fira Mono", ui-monospace, SFMono-Regular, Menlo, monospace;
    }
    button {
      background: var(--accent);
      color: white;
      border: none;
      border-radius: 8px;
      padding: 10px 12px;
      font-size: 14px;
      cursor: pointer;
      transition: background 0.2s ease;
    }
    button:hover {
      background: var(--accent-dark);
    }
    .secondary {
      background: #eef2ff;
      color: #1b3f9b;
    }
    .pill {
      display: inline-flex;
      align-items: center;
      gap: 6px;
      font-size: 12px;
      background: #f0f4ff;
      color: #1b3f9b;
      padding: 4px 8px;
      border-radius: 999px;
      width: fit-content;
    }
    .row {
      display: flex;
      gap: 12px;
      flex-wrap: wrap;
    }
    .status {
      font-size: 13px;
      color: var(--muted);
    }
    .log-actions {
      display: flex;
      justify-content: space-between;
      align-items: center;
    }
    .muted {
      color: var(--muted);
      font-size: 12px;
    }
  </style>
</head>
<body>
  <header>
    <h1>Painel Picornavirus-quali2026</h1>
    <div class=\"subtitle\">Execute o pipeline localmente no WSL/Ubuntu e acompanhe logs em tempo real.</div>
    <span class=\"pill\">WSL + navegador Windows (localhost)</span>
  </header>

  <div class=\"grid\">
    <section class=\"card\">
      <h2>1) Verificar ambiente</h2>
      <p class=\"muted\">Valida dependências essenciais (Velvet, BLAST, Bowtie2, Python, etc.).</p>
      <button onclick=\"runCheckEnv()\">Checar dependências</button>
    </section>

    <section class=\"card\">
      <h2>2) Gerar DEMO reprodutível</h2>
      <p class=\"muted\">Gera FASTQ demo em <code>data/raw</code> (DEMO_R1/R2).</p>
      <button onclick=\"runDemo()\">Gerar DEMO</button>
    </section>

    <section class=\"card\">
      <h2>3) Importar amostra real</h2>
      <label>
        Nome da amostra
        <input id=\"sample-name\" placeholder=\"ex: 81554_S150\" />
      </label>
      <label>
        Caminho do R1 (.fastq.gz)
        <input id=\"sample-r1\" placeholder=\"/mnt/c/Users/.../R1.fastq.gz\" />
      </label>
      <label>
        Caminho do R2 (.fastq.gz)
        <input id=\"sample-r2\" placeholder=\"/mnt/c/Users/.../R2.fastq.gz\" />
      </label>
      <label>
        <span class=\"row\">
          <input id=\"sample-copy\" type=\"checkbox\" />
          Copiar arquivos para <code>data/raw</code> (recomendado)
        </span>
      </label>
      <button onclick=\"runImport()\">Importar amostra</button>
    </section>

    <section class=\"card\">
      <h2>4) Rodar pipeline</h2>
      <label>
        Amostra
        <input id=\"pipeline-sample\" placeholder=\"ex: 81554_S150\" />
      </label>
      <label>
        K-mer (Velvet)
        <input id=\"pipeline-kmer\" type=\"number\" min=\"15\" max=\"127\" value=\"31\" />
      </label>
      <label>
        Assembler
        <select id=\"pipeline-assembler\">
          <option value=\"velvet\">Velvet (padrão)</option>
          <option value=\"spades\">SPAdes (requer config.env)</option>
        </select>
      </label>
      <label>
        <span class=\"row\">
          <input id=\"pipeline-skip-host\" type=\"checkbox\" />
          Pular filtro de hospedeiro
        </span>
      </label>
      <button onclick=\"runPipeline()\">Executar pipeline</button>
    </section>
  </div>

  <section class=\"card\" style=\"margin-top: 20px;\">
    <div class=\"log-actions\">
      <h2>Logs da última tarefa</h2>
      <button class=\"secondary\" onclick=\"refreshLog()\">Atualizar</button>
    </div>
    <div class=\"status\" id=\"job-status\">Nenhuma tarefa em execução.</div>
    <textarea id=\"job-log\" readonly></textarea>
  </section>

  <script>
    let currentJobId = null;
    let pollTimer = null;

    function setStatus(text) {
      document.getElementById('job-status').textContent = text;
    }

    async function startJob(endpoint, payload) {
      const response = await fetch(endpoint, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload || {})
      });
      const data = await response.json();
      if (!response.ok) {
        setStatus(data.error || 'Erro ao iniciar tarefa.');
        return;
      }
      currentJobId = data.job_id;
      setStatus(`Tarefa iniciada: ${data.name} (ID ${currentJobId})`);
      await fetchLog();
      if (pollTimer) {
        clearInterval(pollTimer);
      }
      pollTimer = setInterval(fetchLog, 2500);
    }

    async function fetchLog() {
      if (!currentJobId) {
        return;
      }
      const response = await fetch(`/api/job/${currentJobId}`);
      const data = await response.json();
      if (!response.ok) {
        setStatus(data.error || 'Não foi possível ler o log.');
        return;
      }
      const status = data.status;
      setStatus(`Status: ${status} | Comando: ${data.command}`);
      document.getElementById('job-log').value = data.tail || '';
      if (status === 'completed' || status === 'failed') {
        clearInterval(pollTimer);
      }
    }

    function refreshLog() {
      fetchLog();
    }

    function runCheckEnv() {
      startJob('/api/check');
    }

    function runDemo() {
      startJob('/api/demo');
    }

    function runImport() {
      const sample = document.getElementById('sample-name').value.trim();
      const r1 = document.getElementById('sample-r1').value.trim();
      const r2 = document.getElementById('sample-r2').value.trim();
      const copy = document.getElementById('sample-copy').checked;
      startJob('/api/import', { sample, r1, r2, copy });
    }

    function runPipeline() {
      const sample = document.getElementById('pipeline-sample').value.trim();
      const kmer = document.getElementById('pipeline-kmer').value.trim();
      const assembler = document.getElementById('pipeline-assembler').value;
      const skipHost = document.getElementById('pipeline-skip-host').checked;
      startJob('/api/pipeline', { sample, kmer, assembler, skip_host: skipHost });
    }
  </script>
</body>
</html>
"""


def _tail_file(path, max_lines=200):
    try:
        content = path.read_text(encoding="utf-8", errors="replace")
    except FileNotFoundError:
        return ""
    lines = content.splitlines()
    return "\n".join(lines[-max_lines:])


def _run_job(job_id, command, env):
    job = JOBS[job_id]
    log_path = job["log_path"]
    start_time = time.time()
    with log_path.open("w", encoding="utf-8") as handle:
        handle.write(f"[INFO] Iniciando tarefa {job_id}\n")
        handle.write(f"[INFO] Comando: {command}\n\n")
        handle.flush()
        process = Popen(
            ["bash", "-lc", command],
            cwd=ROOT,
            stdout=handle,
            stderr=STDOUT,
            env=env,
        )
        exit_code = process.wait()
        handle.write(f"\n[INFO] Finalizado com código {exit_code}\n")
    duration = time.time() - start_time
    with JOBS_LOCK:
        job["status"] = "completed" if exit_code == 0 else "failed"
        job["exit_code"] = exit_code
        job["duration"] = round(duration, 1)


def _start_job(name, command, extra_env=None):
    job_id = uuid.uuid4().hex[:10]
    log_path = LOG_DIR / f"ux_dashboard_{job_id}.log"
    base_env = os.environ.copy()
    if extra_env:
        base_env.update(extra_env)
    with JOBS_LOCK:
        JOBS[job_id] = {
            "id": job_id,
            "name": name,
            "command": command,
            "log_path": log_path,
            "status": "running",
            "exit_code": None,
            "duration": None,
        }
    thread = threading.Thread(target=_run_job, args=(job_id, command, base_env), daemon=True)
    thread.start()
    return JOBS[job_id]


class DashboardHandler(BaseHTTPRequestHandler):
    def _send_json(self, payload, status=HTTPStatus.OK):
        data = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _send_html(self, html):
        data = html.encode("utf-8")
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _read_json(self):
        length = int(self.headers.get("Content-Length", "0"))
        if length == 0:
            return {}
        raw = self.rfile.read(length)
        return json.loads(raw.decode("utf-8"))

    def do_GET(self):
        parsed = urlparse(self.path)
        if parsed.path == "/":
            self._send_html(HTML_PAGE)
            return
        if parsed.path.startswith("/api/job/"):
            job_id = parsed.path.split("/")[-1]
            with JOBS_LOCK:
                job = JOBS.get(job_id)
            if not job:
                self._send_json({"error": "Job não encontrado."}, status=HTTPStatus.NOT_FOUND)
                return
            tail = _tail_file(job["log_path"])
            payload = {
                "id": job["id"],
                "name": job["name"],
                "status": job["status"],
                "command": job["command"],
                "exit_code": job["exit_code"],
                "duration": job["duration"],
                "tail": tail,
            }
            self._send_json(payload)
            return
        if parsed.path == "/api/jobs":
            with JOBS_LOCK:
                jobs = list(JOBS.values())[-20:]
            payload = [
                {
                    "id": job["id"],
                    "name": job["name"],
                    "status": job["status"],
                    "command": job["command"],
                }
                for job in jobs
            ]
            self._send_json(payload)
            return
        self.send_error(HTTPStatus.NOT_FOUND, "Not Found")

    def do_POST(self):
        parsed = urlparse(self.path)
        try:
            data = self._read_json()
        except json.JSONDecodeError:
            self._send_json({"error": "JSON inválido."}, status=HTTPStatus.BAD_REQUEST)
            return

        if parsed.path == "/api/check":
            job = _start_job("check-env", "make test-env")
            self._send_json({"job_id": job["id"], "name": job["name"]})
            return
        if parsed.path == "/api/demo":
            job = _start_job("demo", "make demo")
            self._send_json({"job_id": job["id"], "name": job["name"]})
            return
        if parsed.path == "/api/import":
            sample = (data.get("sample") or "").strip()
            r1 = (data.get("r1") or "").strip()
            r2 = (data.get("r2") or "").strip()
            copy_flag = data.get("copy") is True
            if not sample or not r1 or not r2:
                self._send_json({"error": "Sample, R1 e R2 são obrigatórios."}, status=HTTPStatus.BAD_REQUEST)
                return
            copy_arg = " --copy" if copy_flag else ""
            command = (
                "bash scripts/00_import_sample.sh"
                f" --sample {sample} --r1 '{r1}' --r2 '{r2}'{copy_arg}"
            )
            job = _start_job("import-sample", command)
            self._send_json({"job_id": job["id"], "name": job["name"]})
            return
        if parsed.path == "/api/pipeline":
            sample = (data.get("sample") or "").strip()
            kmer = (data.get("kmer") or "").strip()
            assembler = (data.get("assembler") or "velvet").strip()
            skip_host = data.get("skip_host") is True
            if not sample:
                self._send_json({"error": "Sample é obrigatório."}, status=HTTPStatus.BAD_REQUEST)
                return
            args = ["bash scripts/20_run_pipeline.sh", f"--sample {sample}"]
            if kmer:
                args.append(f"--kmer {kmer}")
            if skip_host:
                args.append("--skip-host-filter")
            command = " ".join(args)
            job = _start_job(
                "pipeline",
                command,
                extra_env={"ASSEMBLER": assembler},
            )
            self._send_json({"job_id": job["id"], "name": job["name"]})
            return

        self.send_error(HTTPStatus.NOT_FOUND, "Not Found")


def main():
    parser = argparse.ArgumentParser(description="Painel UX do Picornavirus-quali2026")
    parser.add_argument("--host", default="0.0.0.0", help="Host para bind (default: 0.0.0.0)")
    parser.add_argument("--port", type=int, default=8787, help="Porta do servidor")
    args = parser.parse_args()

    server = ThreadingHTTPServer((args.host, args.port), DashboardHandler)
    print(f"[INFO] Painel rodando em http://{args.host}:{args.port}")
    print("[INFO] Use Ctrl+C para encerrar")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n[INFO] Encerrando servidor...")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
