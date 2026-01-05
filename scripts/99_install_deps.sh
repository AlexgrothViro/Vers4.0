#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REQUIRED_PACKAGES=(
  velvet             # fornece velveth/velvetg
  ncbi-blast+        # fornece blastn/makeblastdb
  bowtie2            # fornece bowtie2/bowtie2-build
)

OPTIONAL_PACKAGES=(
  edirect            # utilitários esearch/efetch para downloads NCBI
)

APT_GET_BIN="$(command -v apt-get || true)"
if [[ -z "${APT_GET_BIN}" ]]; then
  echo "ERRO: apt-get não encontrado. Instalação automática indisponível neste ambiente." >&2
  exit 1
fi

if [[ "${EUID}" -ne 0 ]]; then
  if command -v sudo >/dev/null 2>&1; then
    SUDO="sudo"
  else
    echo "ERRO: precisa ser root ou ter sudo para instalar pacotes com apt-get." >&2
    exit 1
  fi
else
  SUDO=""
fi

on_error() {
  echo "ERRO: falha durante a instalação com apt-get. Verifique conectividade de rede, proxy e repositórios configurados." >&2
}
trap on_error ERR

install_packages() {
  local label="$1"; shift
  local packages=("$@")
  if (( ${#packages[@]} == 0 )); then
    return 0
  fi
  echo "[INFO] Instalando ${label}: ${packages[*]}"
  DEBIAN_FRONTEND=noninteractive ${SUDO} "${APT_GET_BIN}" update -y
  DEBIAN_FRONTEND=noninteractive ${SUDO} "${APT_GET_BIN}" install -y "${packages[@]}"
}

install_packages "dependências obrigatórias" "${REQUIRED_PACKAGES[@]}"

# Não falha se opcionais não existirem; são úteis para downloads.
set +e
install_packages "dependências opcionais" "${OPTIONAL_PACKAGES[@]}"
set -e

echo
echo "[INFO] Instalação concluída. Reexecute scripts/00_check_env.sh para validar o ambiente."
