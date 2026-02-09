const statusEl = document.getElementById("job-status");
const actionEl = document.getElementById("job-action");
const outputEl = document.getElementById("job-output");
const samplesDatalist = document.getElementById("samples");

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

const runAction = async (action, params = {}) => {
  setStatus("Executando...", action);
  setOutput("");

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
    setOutput(data.output || "");

    if (data.status === "running") {
      setStatus("Executando...", action);
      return;
    }

    if (data.status === "done") {
      setStatus("Concluído", action);
    } else {
      setStatus("Falhou", action);
    }

    clearInterval(interval);
    fetchSamples();
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
};

window.addEventListener("load", () => {
  bindButtons();
  fetchSamples();
});
