const statusEl = document.getElementById("job-status");
const actionEl = document.getElementById("job-action");
const outputEl = document.getElementById("job-output");
const samplesDatalist = document.getElementById("samples");
const finalStatusEl = document.getElementById("final-status");
const historyListEl = document.getElementById("history-list");
const pipelineProgressEl = document.getElementById("pipeline-progress");
const reportPreviewEl = document.getElementById("report-preview");
const reportContentEl = document.getElementById("report-content");

const PIPELINE_STEPS = ["host_filter", "assembly", "blast"];

const escapeHtml = (text = "") => text
  .replaceAll("&", "&amp;")
  .replaceAll("<", "&lt;")
  .replaceAll(">", "&gt;");

const setStatus = (status, action) => {
  statusEl.textContent = status;
  if (action) {
    actionEl.textContent = action;
  }
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

const setReportPreview = (text = "") => {
  if (!text) {
    reportPreviewEl.classList.add("hidden");
    reportContentEl.textContent = "";
    return;
  }
  reportContentEl.textContent = text;
  reportPreviewEl.classList.remove("hidden");
};

const resetPipelineProgress = (action) => {
  if (action !== "pipeline") {
    pipelineProgressEl.classList.add("hidden");
    return;
  }
  pipelineProgressEl.classList.remove("hidden");
  PIPELINE_STEPS.forEach((step) => {
    const item = pipelineProgressEl.querySelector(`[data-step="${step}"]`);
    if (!item) return;
    item.classList.remove("running", "done", "error");
    item.querySelector(".step-icon").textContent = "⚪";
  });
};

const markStep = (step, state) => {
  const item = pipelineProgressEl.querySelector(`[data-step="${step}"]`);
  if (!item) return;
  item.classList.remove("running", "done", "error");
  item.classList.add(state);

  const icon = item.querySelector(".step-icon");
  if (state === "running") {
    icon.textContent = "🟡";
  } else if (state === "done") {
    icon.textContent = "✅";
  } else if (state === "error") {
    icon.textContent = "❌";
  } else {
    icon.textContent = "⚪";
  }
};

const updatePipelineProgressFromOutput = (output, finished = false, failed = false) => {
  if (pipelineProgressEl.classList.contains("hidden")) return;

  const hasHost = /\[3\/6\].*Filtrando hospedeiro|\[AVISO\].*Índice do hospedeiro/i.test(output);
  const hasAssembly = /\[4\/6\].*Montagem|run_assembly_router|01_run_velvet/i.test(output);
  const hasBlast = /\[5\/6\].*BLAST|Resultado salvo em:/i.test(output);

  if (hasHost) {
    markStep("host_filter", hasAssembly || finished ? "done" : "running");
  }
  if (hasAssembly) {
    markStep("assembly", hasBlast || finished ? "done" : "running");
  }
  if (hasBlast) {
    markStep("blast", finished ? "done" : "running");
  }

  if (failed) {
    if (hasBlast) {
      markStep("blast", "error");
    } else if (hasAssembly) {
      markStep("assembly", "error");
    } else {
      markStep("host_filter", "error");
    }
  }
};

const fetchSamples = async () => {
  try {
    const response = await fetch("/api/samples");
    if (!response.ok) {
      return;
    }
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
  setReportPreview("");
  resetPipelineProgress("pipeline");
  pollJob(jobId, "pipeline");
};

const formatDate = (value) => {
  if (!value) return "-";
  return value.replace("T", " ");
};

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
    if (!response.ok) {
      return;
    }
    const data = await response.json();
    renderHistory(data.runs || []);
  } catch (err) {
    console.error("Falha ao carregar histórico", err);
  }
};

const fetchReportPreview = async (runDir) => {
  if (!runDir) return;
  try {
    const response = await fetch(`/api/history/file?run=${encodeURIComponent(runDir)}&type=report`);
    if (!response.ok) {
      setReportPreview("Report indisponível para esta execução.");
      return;
    }
    const text = await response.text();
    setReportPreview(text);
  } catch (err) {
    console.error("Falha ao carregar relatório", err);
    setReportPreview("Falha ao carregar report.");
  }
};

const runAction = async (action, params = {}) => {
  setStatus("Executando...", action);
  setOutput("");
  setFinalStatus("");
  setReportPreview("");
  resetPipelineProgress(action);

  const response = await fetch("/api/run", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
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
    updatePipelineProgressFromOutput(output);

    if (data.status === "running") {
      setStatus("Executando...", action);
      return;
    }

    if (data.status === "done") {
      setStatus("Concluído", action);
      updatePipelineProgressFromOutput(output, true, false);
      const report = data.run?.paths?.run_report;
      const runDir = data.run?.run_dir || "";
      const reportLink = report && runDir
        ? `<a href="/api/history/file?run=${encodeURIComponent(runDir)}&type=report" target="_blank">Abrir report</a>`
        : "Report indisponível";
      setFinalStatus(`<strong>SUCESSO</strong> ✅ ${reportLink}`, "ok");
      if (action === "pipeline" && runDir) {
        fetchReportPreview(runDir);
      }
    } else {
      setStatus("Falhou", action);
      updatePipelineProgressFromOutput(output, false, true);
      const tail = data.tail || "Sem log disponível.";
      setFinalStatus(`<strong>FALHA</strong> ❌<pre>${escapeHtml(tail)}</pre>`, "error");
    }

    clearInterval(interval);
    fetchSamples();
    fetchHistory();
  }, 1200);
};

const bindButtons = () => {
  document.querySelectorAll("button[data-action]").forEach((button) => {
    button.addEventListener("click", () => {
      const action = button.dataset.action;
      runAction(action);
    });
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
      if (tab.dataset.tab === "tab-historico") {
        fetchHistory();
      }
    });
  });
};

window.addEventListener("load", () => {
  bindButtons();
  fetchSamples();
  fetchHistory();
  resetPipelineProgress("");
  setReportPreview("");
});
