// Defensive element guards - get elements safely
const getEl = (id) => {
  const el = document.getElementById(id);
  if (!el) console.warn(`Missing element: ${id}`);
  return el;
};

const statusEl = getEl("job-status");
const actionEl = getEl("job-action");
const outputEl = getEl("job-output");
const samplesDatalist = getEl("samples");
const sampleSelectEl = getEl("sample-select");
const finalStatusEl = getEl("final-status");
const historyListEl = getEl("history-list");
const pipelineProgressEl = getEl("pipeline-progress");
const reportViewerEl = getEl("report-viewer");
const reportContentEl = getEl("report-content");
const dbTargetEl = getEl("db-target");
const dbStatusEl = getEl("db-status");
const configFormEl = getEl("config-form");
const configStatusEl = getEl("config-status");
const rebuildEnvBtnEl = getEl("rebuild-env-btn");
const rebuildStatusEl = getEl("rebuild-status");
const envFileStatusEl = getEl("env-file-status");
const envMtimeEl = getEl("env-mtime");
const envPathEl = getEl("env-path");

// Constants
const JOB_POLL_INTERVAL_MS = 1200;

// State to remember last selected DB
let currentDB = {
  target: null,
  query: null,
  taxid: null,
  ncbi_db: null
};

// Early abort if essential elements are missing
if (!statusEl || !outputEl || !finalStatusEl) {
  console.error("Critical UI elements missing. Dashboard cannot initialize.");
}

const stageConfig = {
  host: { marker: "[3/6]", label: "Remoção de hospedeiro" },
  assembly: { marker: "[4/6]", label: "Montagem" },
  blast: { marker: "[5/6]", label: "BLAST" },
};

const stageIcon = { pending: "⏳", running: "🔄", done: "✅", error: "❌", skipped: "⏭️" };

const setStatus = (status, action) => {
  if (statusEl) statusEl.textContent = status;
  if (action && actionEl) actionEl.textContent = action;
};

const updateDBStatus = () => {
  if (!dbStatusEl) return;
  if (currentDB.target) {
    const displayText = currentDB.target === "ptv" ? "Teschovirus A (PTV)" :
                       currentDB.target === "psv" ? "Sapelovirus A" :
                       currentDB.target === "evg" ? "Enterovirus G" :
                       currentDB.target === "teschovirus_a" ? "Teschovirus A (PTV)" :
                       currentDB.target === "sapelovirus_a" ? "Sapelovirus A" :
                       currentDB.target === "enterovirus_g" ? "Enterovirus G" :
                       currentDB.target;
    dbStatusEl.textContent = `DB ativo: ${displayText}`;
    dbStatusEl.className = "db-status active";
  } else {
    dbStatusEl.textContent = "DB: nenhum selecionado";
    dbStatusEl.className = "db-status";
  }
};

const setOutput = (text) => {
  if (!outputEl) return;
  outputEl.textContent = text || "";
  outputEl.scrollTop = outputEl.scrollHeight;
};

const setFinalStatus = (html, tone = "") => {
  if (!finalStatusEl) return;
  finalStatusEl.className = `final-status ${tone}`.trim();
  finalStatusEl.innerHTML = html || "";
};

const setReportContent = (content) => {
  if (!reportContentEl || !reportViewerEl) return;
  reportContentEl.textContent = content || "";
  reportViewerEl.hidden = !content;
};

const initPipelineProgress = () => {
  if (!pipelineProgressEl) return;
  pipelineProgressEl.querySelectorAll("li[data-stage]").forEach((item) => {
    item.dataset.state = "pending";
    item.querySelector(".stage-icon").textContent = stageIcon.pending;
  });
};

const showPipelineProgress = (show) => {
  if (!pipelineProgressEl) return;
  pipelineProgressEl.hidden = !show;
  if (show) initPipelineProgress();
};

const applyStageState = (stage, state) => {
  if (!pipelineProgressEl) return;
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

  if (jobStatus === "done") applyStageState("blast", "done");
  if (jobStatus === "error") {
    if (hasBlast) applyStageState("blast", "error");
    else if (hasAssembly) applyStageState("assembly", "error");
    else if (hasHost || skippedHost) applyStageState("host", "error");
  }
};

const fetchTargets = async () => {
  if (!dbTargetEl) return;
  try {
    const response = await fetch("/api/targets");
    if (!response.ok) return;
    const data = await response.json();
    dbTargetEl.innerHTML = "";
    (data.targets || []).forEach((target) => {
      const option = document.createElement("option");
      option.value = target.key;
      option.textContent = `${target.display_name} (${target.key})`;
      dbTargetEl.appendChild(option);
    });
  } catch (err) {
    console.error("Falha ao listar targets", err);
  }
};

const fetchSamples = async () => {
  try {
    const response = await fetch("/api/samples");
    if (!response.ok) return;
    const data = await response.json();
    if (samplesDatalist) samplesDatalist.innerHTML = "";
    if (sampleSelectEl) sampleSelectEl.innerHTML = '<option value="">Selecione uma amostra</option>';
    (data.samples || []).forEach((sample) => {
      if (samplesDatalist) {
        const option = document.createElement("option");
        option.value = sample;
        samplesDatalist.appendChild(option);
      }

      if (sampleSelectEl) {
        const selectOption = document.createElement("option");
        selectOption.value = sample;
        selectOption.textContent = sample;
        sampleSelectEl.appendChild(selectOption);
      }
    });
  } catch (err) {
    console.error("Falha ao listar amostras", err);
  }
};

const openRunArtifact = (runDir, fileType) => {
  window.open(`/api/history/file?run=${encodeURIComponent(runDir)}&type=${encodeURIComponent(fileType)}`, "_blank");
};

const loadReportInline = async (runDir) => {
  if (!runDir) return setReportContent("");
  try {
    const response = await fetch(`/api/history/file?run=${encodeURIComponent(runDir)}&type=report`);
    if (!response.ok) return setReportContent("");
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
  if (!response.ok) return alert(`Falha ao reexecutar: ${await response.text()}`);
  const { job_id: jobId } = await response.json();
  setStatus("Executando...", `rerun:${runDir}`);
  setOutput("");
  setFinalStatus("");
  showPipelineProgress(true);
  pollJob(jobId, `rerun:${runDir}`);
};

const renderHistory = (runs) => {
  if (!historyListEl) return;
  if (!runs.length) {
    historyListEl.innerHTML = "<p>Nenhuma execução registrada ainda.</p>";
    return;
  }
  historyListEl.innerHTML = "";
  runs.forEach((run) => {
    const card = document.createElement("article");
    card.className = "history-item";
    const statusClass = run.exit_code === 0 ? "ok" : "error";
    card.innerHTML = `
      <header><h3>${run.sample || "(sem sample)"}</h3><span class="badge ${statusClass}">${run.exit_code === 0 ? "SUCESSO" : "FALHA"}</span></header>
      <p><strong>Início:</strong> ${(run.start || "-").replace("T", " ")} • <strong>Fim:</strong> ${(run.end || "-").replace("T", " ")}</p>
      <div class="actions">
        <button data-open="report">Abrir report</button>
        <button data-open="log">Abrir log</button>
        <button data-open="blast">Abrir blast</button>
        <button data-rerun="1">Reexecutar</button>
      </div>`;

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
    renderHistory((await response.json()).runs || []);
  } catch (err) {
    console.error("Falha ao carregar histórico", err);
  }
};

const loadEnvironmentStatus = async () => {
  try {
    const response = await fetch("/api/config/environment");
    if (!response.ok) {
      if (envFileStatusEl) envFileStatusEl.textContent = "Erro ao carregar";
      return;
    }
    const data = await response.json();
    
    if (envFileStatusEl) {
      envFileStatusEl.textContent = data.has_environment_yml ? "Encontrado ✅" : "Não encontrado ❌";
    }
    if (envMtimeEl) {
      envMtimeEl.textContent = data.environment_yml_mtime ? data.environment_yml_mtime.replace("T", " ") : "N/A";
    }
    if (envPathEl) {
      envPathEl.textContent = data.environment_yml_path || "N/A";
    }
  } catch (err) {
    console.error("Falha ao carregar status do ambiente", err);
    if (envFileStatusEl) envFileStatusEl.textContent = "Erro ao carregar";
  }
};

const loadConfigEnv = async () => {
  if (!configFormEl) return;
  
  try {
    const response = await fetch("/api/config/env");
    if (!response.ok) {
      console.error("Falha ao carregar configuração");
      return;
    }
    const data = await response.json();
    const config = data.config || {};
    
    // Populate form fields
    Object.keys(config).forEach((key) => {
      const input = configFormEl.querySelector(`[name="${key}"]`);
      if (input) {
        input.value = config[key] || "";
      }
    });
  } catch (err) {
    console.error("Falha ao carregar configuração", err);
  }
};

const saveConfigEnv = async (event) => {
  event.preventDefault();
  if (!configFormEl || !configStatusEl) return;
  
  const formData = new FormData(configFormEl);
  const config = {};
  
  // Collect non-empty values
  for (const [key, value] of formData.entries()) {
    if (value.trim()) {
      config[key] = value.trim();
    }
  }
  
  try {
    const response = await fetch("/api/config/env", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ config }),
    });
    
    const result = await response.json();
    
    if (response.ok && result.success) {
      configStatusEl.className = "final-status ok";
      configStatusEl.innerHTML = `<strong>SUCESSO</strong> ✅ ${result.message}`;
    } else {
      configStatusEl.className = "final-status error";
      configStatusEl.innerHTML = `<strong>ERRO</strong> ❌ ${result.error || "Falha ao salvar configuração"}`;
    }
  } catch (err) {
    configStatusEl.className = "final-status error";
    configStatusEl.innerHTML = `<strong>ERRO</strong> ❌ ${err.message}`;
  }
};

const rebuildEnvironment = async () => {
  if (!rebuildEnvBtnEl || !rebuildStatusEl) return;
  
  rebuildStatusEl.className = "final-status";
  rebuildStatusEl.innerHTML = "";
  
  try {
    const response = await fetch("/api/config/environment/rebuild", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({}),
    });
    
    if (!response.ok) {
      const result = await response.json();
      rebuildStatusEl.className = "final-status error";
      rebuildStatusEl.innerHTML = `<strong>ERRO</strong> ❌ ${result.error || "Falha ao iniciar recriação"}`;
      return;
    }
    
    const { job_id: jobId } = await response.json();
    rebuildStatusEl.className = "final-status ok";
    rebuildStatusEl.innerHTML = "<strong>Recriando ambiente...</strong> ⏳ Aguarde a conclusão.";
    
    // Poll the job
    const interval = setInterval(async () => {
      const jobResponse = await fetch(`/api/job/${jobId}`);
      if (!jobResponse.ok) {
        clearInterval(interval);
        rebuildStatusEl.className = "final-status error";
        rebuildStatusEl.innerHTML = "<strong>ERRO</strong> ❌ Falha ao consultar status do job";
        return;
      }
      
      const jobData = await jobResponse.json();
      
      if (jobData.status === "running" || jobData.status === "queued") {
        return;
      }
      
      clearInterval(interval);
      
      if (jobData.status === "done") {
        rebuildStatusEl.className = "final-status ok";
        rebuildStatusEl.innerHTML = "<strong>SUCESSO</strong> ✅ Ambiente recriado com sucesso";
        loadEnvironmentStatus();
      } else {
        rebuildStatusEl.className = "final-status error";
        rebuildStatusEl.innerHTML = `<strong>ERRO</strong> ❌<pre>${jobData.tail || "Falha ao recriar ambiente"}</pre>`;
      }
    }, JOB_POLL_INTERVAL_MS);
  } catch (err) {
    rebuildStatusEl.className = "final-status error";
    rebuildStatusEl.innerHTML = `<strong>ERRO</strong> ❌ ${err.message}`;
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
    setStatus("Erro ao iniciar", action);
    setOutput((await response.text()) || "Falha ao iniciar");
    return;
  }
  pollJob((await response.json()).job_id, action);
};

const pollJob = (jobId, action) => {
  const interval = setInterval(async () => {
    const response = await fetch(`/api/job/${jobId}`);
    if (!response.ok) {
      setStatus("Erro ao consultar", action);
      return clearInterval(interval);
    }
    const data = await response.json();
    setOutput(data.output || "");

    if ((action === "pipeline" || action.startsWith("rerun:")) && data.status !== "queued") {
      showPipelineProgress(true);
      updatePipelineProgressFromOutput(data.output || "", data.status);
    } else showPipelineProgress(false);

    if (data.status === "running" || data.status === "queued") return setStatus("Executando...", action);

    if (data.status === "done") {
      setStatus("Concluído", action);
      const runDir = data.run?.run_dir || "";
      const reportLink = runDir ? `<a href="/api/history/file?run=${encodeURIComponent(runDir)}&type=report" target="_blank">Abrir report</a>` : "";
      setFinalStatus(`<strong>SUCESSO</strong> ✅ ${reportLink}`, "ok");
      await loadReportInline(runDir);
    } else {
      setStatus("Falhou", action);
      setFinalStatus(`<strong>FALHA</strong> ❌<pre>${data.tail || "Sem log disponível."}</pre>`, "error");
    }

    clearInterval(interval);
    fetchSamples();
    fetchHistory();
  }, JOB_POLL_INTERVAL_MS);
};

const uploadImport = async (formData) => {
  setStatus("Executando...", "upload_import");
  setOutput("");
  setFinalStatus("");

  const response = await fetch("/api/import-upload", { method: "POST", body: formData });
  if (!response.ok) {
    const msg = await response.text();
    setFinalStatus(`<strong>FALHA</strong> ❌<pre>${msg}</pre>`, "error");
    setStatus("Falhou", "upload_import");
    return;
  }
  const data = await response.json();
  setStatus("Concluído", "upload_import");
  setFinalStatus(`<strong>SUCESSO</strong> ✅ ${data.message}`, "ok");
  fetchSamples();
};

const bindButtons = () => {
  document.querySelectorAll("button[data-action]").forEach((button) => {
    button.addEventListener("click", () => runAction(button.dataset.action));
  });

  document.getElementById("db-form").addEventListener("submit", (event) => {
    event.preventDefault();
    const formData = new FormData(event.target);
    const target = formData.get("target");
    const query = formData.get("query");
    const taxid = formData.get("taxid");
    
    // Remember the selected DB configuration
    currentDB = {
      target: target,
      query: query,
      taxid: taxid,
      ncbi_db: null
    };
    
    runAction("build_db", {
      target: target,
      query: query,
      taxid: taxid,
    });
    
    updateDBStatus();
  });

  document.getElementById("import-form").addEventListener("submit", (event) => {
    event.preventDefault();
    const formData = new FormData(event.target);
    runAction("import_sample", {
      sample: formData.get("sample"),
      r1: formData.get("r1"),
      r2: formData.get("r2"),
      copy: formData.get("copy") === "on",
    });
  });

  document.getElementById("upload-form").addEventListener("submit", (event) => {
    event.preventDefault();
    uploadImport(new FormData(event.target));
  });

  document.getElementById("pipeline-form").addEventListener("submit", (event) => {
    event.preventDefault();
    const formData = new FormData(event.target);
    const params = {
      sample: formData.get("sample_select") || formData.get("sample"),
      assembler: formData.get("assembler"),
      kmer: formData.get("kmer"),
    };
    
    // Include DB configuration if one was selected
    if (currentDB.target) {
      params.db = currentDB.target;
      if (currentDB.query) {
        params.db_query = currentDB.query;
      }
      if (currentDB.ncbi_db) {
        params.ncbi_db = currentDB.ncbi_db;
      }
    }
    
    runAction("pipeline", params);
  });

  document.querySelectorAll(".tab").forEach((tab) => {
    tab.addEventListener("click", () => {
      document.querySelectorAll(".tab").forEach((item) => item.classList.remove("active"));
      document.querySelectorAll(".tab-panel").forEach((panel) => panel.classList.remove("active"));
      tab.classList.add("active");
      document.getElementById(tab.dataset.tab).classList.add("active");
      if (tab.dataset.tab === "tab-historico") fetchHistory();
      if (tab.dataset.tab === "tab-configuracao") {
        loadEnvironmentStatus();
        loadConfigEnv();
      }
    });
  });

  // Bind config form
  if (configFormEl) {
    configFormEl.addEventListener("submit", saveConfigEnv);
  }
  
  // Bind rebuild environment button
  if (rebuildEnvBtnEl) {
    rebuildEnvBtnEl.addEventListener("click", rebuildEnvironment);
  }
};

window.addEventListener("load", () => {
  bindButtons();
  fetchTargets();
  fetchSamples();
  fetchHistory();
  showPipelineProgress(false);
  setReportContent("");
  updateDBStatus();
});
