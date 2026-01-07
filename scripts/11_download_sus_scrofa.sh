#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

if [[ -f "${REPO_ROOT}/config.env" ]]; then
  source "${REPO_ROOT}/config.env"
fi

source "${SCRIPT_DIR}/lib/common.sh"

DEST_DIR="ref/host"
FA="${DEST_DIR}/sus_scrofa.fa"

# Query para pegar sequências genômicas de Sus scrofa em RefSeq
QUERY='"Sus scrofa"[Organism] AND srcdb_refseq[PROP] AND biomol_genomic[PROP]'

echo "============================================="
echo "  Download do genoma de Sus scrofa (hospedeiro)"
echo "  Fonte: NCBI nuccore via EDirect"
echo "============================================="
echo "Destino: ${FA}"
echo
echo "Query NCBI:"
echo "  ${QUERY}"
echo

mkdir -p "$DEST_DIR"

log_info "Buscando sequências no NCBI..."
esearch -db nucleotide -query "$QUERY" | efetch -format fasta > "$FA"

if [[ ! -s "$FA" ]]; then
  log_error "arquivo FASTA ${FA} está vazio. Verifique a query ou sua conexão."
fi

log_info "Genoma de Sus scrofa salvo em: ${FA}"
echo "Pré-visualização das primeiras linhas:"
head -n 5 "$FA" || true
