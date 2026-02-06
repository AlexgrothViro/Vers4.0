#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUNDLE_DIR="$ROOT/.bundle"
BIN_DIR="$BUNDLE_DIR/bin"
ENV_DIR="$BUNDLE_DIR/env"
MAMBA_ROOT="$BUNDLE_DIR/mamba"
MICRO="$BIN_DIR/micromamba"

export MAMBA_ROOT_PREFIX="$MAMBA_ROOT"

need_cmd() { command -v "$1" >/dev/null 2>&1; }

echo "[INFO] Root: $ROOT"
mkdir -p "$BIN_DIR" "$MAMBA_ROOT"

# Pré-requisitos mínimos (se faltar, tenta instalar via apt com sudo)
if ! need_cmd curl; then
  echo "[WARN] curl não encontrado. Tentando instalar via apt (sudo)..."
  sudo apt-get update && sudo apt-get install -y curl ca-certificates
fi
if ! need_cmd tar; then
  echo "[WARN] tar não encontrado. Tentando instalar via apt (sudo)..."
  sudo apt-get update && sudo apt-get install -y tar
fi
if ! need_cmd bzip2; then
  echo "[WARN] bzip2 não encontrado. Tentando instalar via apt (sudo)..."
  sudo apt-get update && sudo apt-get install -y bzip2
fi
# Baixa micromamba (URL oficial "latest")
if [[ ! -x "$MICRO" ]]; then
  echo "[INFO] Baixando micromamba..."
  tmp="$(mktemp -d)"
  ( cd "$tmp" && curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba )
  install -m 0755 "$tmp/bin/micromamba" "$MICRO"
  rm -rf "$tmp"
else
  echo "[INFO] micromamba já existe: $MICRO"
fi

# Cria/atualiza o ambiente local a partir do environment.yml do projeto
if [[ ! -f "$ROOT/environment.yml" ]]; then
  echo "[ERROR] environment.yml não encontrado em $ROOT"
  exit 1
fi
if [[ ! -d "$ENV_DIR" ]]; then
  echo "[INFO] Criando ambiente em: $ENV_DIR"
  "$MICRO" create -y -p "$ENV_DIR" -f "$ROOT/environment.yml"
else
  echo "[INFO] Atualizando ambiente existente em: $ENV_DIR"
  "$MICRO" install -y -p "$ENV_DIR" -f "$ROOT/environment.yml"
fi

echo "[INFO] Verificando ferramentas dentro do env..."
"$MICRO" run -p "$ENV_DIR" python -V || true
"$MICRO" run -p "$ENV_DIR" blastn -version 2>/dev/null | head -n 1 || true
"$MICRO" run -p "$ENV_DIR" bowtie2 --version 2>/dev/null | head -n 1 || true
"$MICRO" run -p "$ENV_DIR" spades.py --version 2>/dev/null || true
"$MICRO" run -p "$ENV_DIR" mafft --version 2>&1 | head -n 1 || true
"$MICRO" run -p "$ENV_DIR" iqtree2 -version 2>&1 | head -n 1 || true

echo "[OK] Ambiente pronto. Use: bash bundle/run.sh smoke-test"
