#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

if [[ -f "${REPO_ROOT}/config.env" ]]; then
  source "${REPO_ROOT}/config.env"
fi

source "${SCRIPT_DIR}/lib/common.sh"

if [[ ! -f "${REPO_ROOT}/config.env" ]]; then
  log_error "config.env não existe. Crie com: cp config.env.example config.env"
fi

command -v esearch >/dev/null 2>&1 || log_error "EDirect não encontrado (esearch). Instale EDirect."
command -v efetch  >/dev/null 2>&1 || log_error "EDirect não encontrado (efetch). Instale EDirect."
command -v makeblastdb >/dev/null 2>&1 || log_error "makeblastdb não encontrado (blast+)."

DB_ROOT_RESOLVED="$(resolve_path "${DB_ROOT:-data/db}")"
DB_DIR="${DB_ROOT_RESOLVED}/${DB_NAME}"
FASTA="${DB_DIR}/${DB_NAME}.fasta"
BLASTDB="${DB_DIR}/${DB_NAME}"

mkdir -p "$DB_DIR"

QUERY="${NCBI_QUERY:-txid${TARGET_TAXID}[Organism:exp] AND refseq[filter]}"
RETMAX="${EDIRECT_RETMAX:-500}"

log_info "[DB] Query: $QUERY"
log_info "[DB] retmax: $RETMAX"
log_info "[DB] Baixando FASTA..."

esearch -db nucleotide -query "$QUERY" -retmax "$RETMAX" | \
  efetch -format fasta > "$FASTA"

log_info "[DB] Sequências baixadas: $(grep -c '^>' "$FASTA" || true)"
log_info "[DB] Construindo BLAST DB..."
makeblastdb -in "$FASTA" -dbtype nucl -out "$BLASTDB" -parse_seqids

log_info "[DB] OK: $BLASTDB"
