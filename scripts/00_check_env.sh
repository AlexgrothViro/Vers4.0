#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  echo "Uso: $0 [--install]" >&2
  echo "  --install    tenta instalar dependências obrigatórias via apt-get" >&2
}

AUTO_INSTALL=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --install)
      AUTO_INSTALL=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      usage
      exit 1
      ;;
  esac
done

# Programas obrigatórios para os testes básicos
REQUIRED_CMDS=(
  gcc
  make
  velveth
  velvetg
  blastn
  makeblastdb
  bowtie2
  bowtie2-build
  python3
  dos2unix
  mafft
  fasttree
  esearch
  efetch
)

IQTREE_CANDIDATES=(iqtree iqtree2)

# Programas opcionais (úteis para download via HTTP)
OPTIONAL_CMDS=(curl)

echo "== Verificando programas necessários no PATH =="

MISSING=0
for cmd in "${REQUIRED_CMDS[@]}"; do
  if command -v "$cmd" >/dev/null 2>&1; then
    path=$(command -v "$cmd")
    printf "  [OK] %s encontrado em %s\n" "$cmd" "$path"
  else
    echo "  [FALTA] $cmd não está no PATH"
    MISSING=1
  fi
done

IQTREE_FOUND=0
for cmd in "${IQTREE_CANDIDATES[@]}"; do
  if command -v "$cmd" >/dev/null 2>&1; then
    path=$(command -v "$cmd")
    printf "  [OK] %s encontrado em %s (iqtree/iqtree2)\n" "$cmd" "$path"
    IQTREE_FOUND=1
    break
  fi
done
if [[ $IQTREE_FOUND -eq 0 ]]; then
  echo "  [FALTA] iqtree ou iqtree2 não encontrado no PATH"
  MISSING=1
fi

for cmd in "${OPTIONAL_CMDS[@]}"; do
  if command -v "$cmd" >/dev/null 2>&1; then
    path=$(command -v "$cmd")
    printf "  [OK] %s encontrado em %s (opcional)\n" "$cmd" "$path"
  else
    echo "  [OPCIONAL] $cmd não está no PATH (usado em downloads NCBI/HTTP)"
  fi
done

echo
PTV_FASTA="data/ptv_db.fa"
if [[ ! -s "$PTV_FASTA" ]]; then
  echo "ATENÇÃO: FASTA de referência ausente em $PTV_FASTA."
  echo "  Sugestão: make ptv-fasta-legacy   # cria data/ptv_db.fa a partir de data/ref/ptv_db.fa"
else
  echo "FASTA de referência encontrado: $PTV_FASTA"
fi

if ! command -v esearch >/dev/null 2>&1 || ! command -v efetch >/dev/null 2>&1; then
  echo "AVISO: EDirect ausente (esearch/efetch) - necessário para baixar FASTA do NCBI."
fi

if [[ -n "${BLAST_DB:-}" ]]; then
  echo "BLAST_DB está definido como: $BLAST_DB"
  missing_db_files=()
  for ext in nhr nin nsq; do
    file="${BLAST_DB}.${ext}"
    [[ -f "$file" ]] || missing_db_files+=("$file")
  done
  if (( ${#missing_db_files[@]} > 0 )); then
    echo "  [AVISO] Índices BLAST ausentes:"
    printf '    - %s\n' "${missing_db_files[@]}"
    echo "  Sugestão: make blastdb"
  else
    echo "  [OK] Índices BLAST encontrados."
  fi
else
  echo "ATENÇÃO: variável de ambiente BLAST_DB não definida."
  echo "  Sugestão: make blastdb   # gera blastdb/ptv e exporte BLAST_DB=blastdb/ptv"
fi

if [[ $MISSING -ne 0 ]]; then
  echo
  if [[ $AUTO_INSTALL -eq 1 ]]; then
    echo "[AÇÃO] Instalando dependências faltantes via ${SCRIPT_DIR}/99_install_deps.sh"
    if "${SCRIPT_DIR}/99_install_deps.sh"; then
      echo
      echo "[INFO] Reavaliando ambiente após instalação..."
      exec "$0"
    else
      echo "Falha na instalação automática. Confira os logs acima." >&2
      exit 1
    fi
  fi
  echo "Alguns programas estão faltando. Veja a seção 'Ambiente e requisitos' no README.md." >&2
  exit 1
fi

echo
echo "Ambiente básico OK."
