const statusEl = document.getElementById("job-status");
const actionEl = document.getElementById("job-action");
const outputEl = document.getElementById("job-output");
const samplesDatalist = document.getElementById("samples");
const finalStatusEl = document.getElementById("final-status");
const historyListEl = document.getElementById("history-list");
const pipelineProgressEl = document.getElementById("pipeline-progress");
const reportViewerEl = document.getElementById("report-viewer");
const reportContentEl = document.getElementById("report-content");

const stageConfig = {
  host: { marker: "[3/6]", label: "Remoção de hospedeiro" },
  assembly: { marker: "[4/6]", label: "Montagem" },
  blast: { marker: "[5/6]", label: "BLAST" },
};

const stageIcon = {
  pending: "⏳",
  running: "🔄",
  done: "✅",
  error: "❌",
  skipped: "⏭️",
};

const setStatus = (status, action) => {
  statusEl.textContent = status;
  if (action) actionEl.textContent = action;
};

const setOutput = (text) => {
  outputEl.textContent = text || "";
  outputEl.scrollTop = outputEl.scrollHeight;
};

const appendOutput = (text) => {
  outputEl.textContent += text;
  outputEl.scrollTop = outputEl.scrollHeight;
};

const setFinalStatus = (html, tone = "") => {
  finalStatusEl.className = `final-status ${tone}`.trim();
  finalStatusEl.innerHTML = html || "";
};

const setReportContent = (content) => {
  reportContentEl.textContent = content || "";
  reportViewerEl.hidden = !content;
};

const initPipelineProgress = () => {
  pipelineProgressEl.querySelectorAll("li[data-stage]").forEach((item) => {
    item.dataset.state = "pending";
    item.querySelector(".stage-icon").textContent = stageIcon.pending;
  });
};

const showPipelineProgress = (show) => {
  pipelineProgressEl.hidden = !show;
  if (show) initPipelineProgress();
};

const applyStageState = (stage, state) => {
  const item = pipelineProgressEl.querySelector(`li[data-stage="${stage}"]`);
  if (!item) return;
  item.dataset.state = state;
  item.querySelector(".stage-icon").textContent = stageIcon[state] || stageIcon.pending;
};

const updatePipelineProgressFromOutput = (output, jobStatus) => {
  const hasHost = output.includes(stageConfig.host.marker);
  const hasAssembly = output.includes(stageConfig.assembly.marker);
  const hasBlast = output.includes(stageConfig.blast.marker);
  const skippedHost = output.includes("[AVISO] Índice do hospedeiro não encontrado");

  applyStageState("host", skippedHost ? "skipped" : (hasHost ? "done" : "pending"));
  applyStageState("assembly", hasAssembly ? "done" : (hasHost || skippedHost ? "running" : "pending"));
  applyStageState("blast", hasBlast ? "running" : "pending");

  if (jobStatus === "done") {
    applyStageState("blast", "done");
  }
  if (jobStatus === "error") {
    if (hasBlast) {
      applyStageState("blast", "error");
    } else if (hasAssembly) {
      applyStageState("assembly", "error");
    } else if (hasHost || skippedHost) {
      applyStageState("host", "error");
    }
  }
};

const fetchSamples = async () => {
  try {
    const response = await fetch("/api/samples");
    if (!response.ok) return;
    const data = await response.json();
    samplesDatalist.innerHTML = "";
    data.samples.forEach((sample) => {
      const option = document.createElement("option");
      option.value = sample;
      samplesDatalist.appendChild(option);
    });
  } catch (err) {
    console.error("Falha ao listar amostras", err);
  }
};

const openRunArtifact = (runDir, fileType) => {
  window.open(`/api/history/file?run=${encodeURIComponent(runDir)}&type=${encodeURIComponent(fileType)}`, "_blank");
};

const loadReportInline = async (runDir) => {
  if (!runDir) {
    setReportContent("");
    return;
  }
  try {
    const response = await fetch(`/api/history/file?run=${encodeURIComponent(runDir)}&type=report`);
    if (!response.ok) {
      setReportContent("");
      return;
    }
    setReportContent(await response.text());
  } catch (err) {
    console.error("Falha ao carregar resumo", err);
    setReportContent("");
  }
};

const rerunHistory = async (runDir) => {
  const response = await fetch("/api/history/rerun", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ run_dir: runDir }),
  });
  if (!response.ok) {
    const text = await response.text();
    alert(`Falha ao reexecutar: ${text}`);
    return;
  }
  const { job_id: jobId } = await response.json();
  setStatus("Executando...", `rerun:${runDir}`);
  setOutput("");
  setFinalStatus("");
  setReportContent("");
  showPipelineProgress(true);
  pollJob(jobId, `rerun:${runDir}`);
};

const formatDate = (value) => (value ? value.replace("T", " ") : "-");

const renderHistory = (runs) => {
  if (!runs.length) {
    historyListEl.innerHTML = "<p>Nenhuma execução registrada ainda.</p>";
    return;
  }

  historyListEl.innerHTML = "";
  runs.forEach((run) => {
    const card = document.createElement("article");
    card.className = "history-item";
    const sample = run.sample || "(sem sample)";
    const status = run.exit_code === 0 ? "SUCESSO" : "FALHA";
    const statusClass = run.exit_code === 0 ? "ok" : "error";

    card.innerHTML = `
      <header>
        <h3>${sample}</h3>
        <span class="badge ${statusClass}">${status}</span>
      </header>
      <p><strong>Início:</strong> ${formatDate(run.start)} • <strong>Fim:</strong> ${formatDate(run.end)}</p>
      <p><strong>Assembler:</strong> ${run.assembler || "-"} • <strong>k-mer:</strong> ${run.kmer || "-"} • <strong>Threads:</strong> ${run.threads || "-"}</p>
      <div class="actions">
        <button data-open="report">Abrir report</button>
        <button data-open="log">Abrir log</button>
        <button data-open="blast">Abrir blast</button>
        <button data-rerun="1">Reexecutar</button>
      </div>
    `;

    card.querySelectorAll("button[data-open]").forEach((button) => {
      button.addEventListener("click", () => openRunArtifact(run.run_dir, button.dataset.open));
    });
    card.querySelector("button[data-rerun]").addEventListener("click", () => rerunHistory(run.run_dir));

    historyListEl.appendChild(card);
  });
};

const fetchHistory = async () => {
  try {
    const response = await fetch("/api/history");
    if (!response.ok) return;
    const data = await response.json();
    renderHistory(data.runs || []);
  } catch (err) {
    console.error("Falha ao carregar histórico", err);
  }
};

const runAction = async (action, params = {}) => {
  setStatus("Executando...", action);
  setOutput("");
  setFinalStatus("");
  setReportContent("");
  showPipelineProgress(action === "pipeline");

  const response = await fetch("/api/run", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ action, params }),
  });

  if (!response.ok) {
    const errorText = await response.text();
    setStatus("Erro ao iniciar", action);
    appendOutput(errorText || "Falha ao iniciar a execução.\n");
    return;
  }

  const { job_id: jobId } = await response.json();
  pollJob(jobId, action);
};

const pollJob = (jobId, action) => {
  const interval = setInterval(async () => {
    const response = await fetch(`/api/job/${jobId}`);
    if (!response.ok) {
      setStatus("Erro ao consultar", action);
      clearInterval(interval);
      return;
    }
    const data = await response.json();
    const output = data.output || "";
    setOutput(output);

    if ((action === "pipeline" || action.startsWith("rerun:")) && data.status !== "queued") {
      showPipelineProgress(true);
      updatePipelineProgressFromOutput(output, data.status);
    } else {
      showPipelineProgress(false);
    }

    if (data.status === "running" || data.status === "queued") {
      setStatus("Executando...", action);
      return;
    }

    if (data.status === "done") {
      setStatus("Concluído", action);
      const runDir = data.run?.run_dir || "";
      const reportLink = runDir
        ? `<a href="/api/history/file?run=${encodeURIComponent(runDir)}&type=report" target="_blank">Abrir report</a>`
        : "Report indisponível";
      setFinalStatus(`<strong>SUCESSO</strong> ✅ ${reportLink}`, "ok");
      await loadReportInline(runDir);
    } else {
      setStatus("Falhou", action);
      const tail = data.tail || "Sem log disponível.";
      setFinalStatus(`<strong>FALHA</strong> ❌<pre>${tail}</pre>`, "error");
    }

    clearInterval(interval);
    fetchSamples();
    fetchHistory();
  }, 1200);
};

const bindButtons = () => {
  document.querySelectorAll("button[data-action]").forEach((button) => {
    button.addEventListener("click", () => runAction(button.dataset.action));
  });

  const importForm = document.getElementById("import-form");
  importForm.addEventListener("submit", (event) => {
    event.preventDefault();
    const formData = new FormData(importForm);
    runAction("import_sample", {
      sample: formData.get("sample"),
      r1: formData.get("r1"),
      r2: formData.get("r2"),
      copy: formData.get("copy") === "on",
    });
  });

  const pipelineForm = document.getElementById("pipeline-form");
  pipelineForm.addEventListener("submit", (event) => {
    event.preventDefault();
    const formData = new FormData(pipelineForm);
    runAction("pipeline", {
      sample: formData.get("sample"),
      assembler: formData.get("assembler"),
      kmer: formData.get("kmer"),
    });
  });

  document.querySelectorAll(".tab").forEach((tab) => {
    tab.addEventListener("click", () => {
      document.querySelectorAll(".tab").forEach((item) => item.classList.remove("active"));
      document.querySelectorAll(".tab-panel").forEach((panel) => panel.classList.remove("active"));
      tab.classList.add("active");
      document.getElementById(tab.dataset.tab).classList.add("active");
      if (tab.dataset.tab === "tab-historico") fetchHistory();
    });
  });
};

window.addEventListener("load", () => {
  bindButtons();
  fetchSamples();
  fetchHistory();
  showPipelineProgress(false);
  setReportContent("");
});
