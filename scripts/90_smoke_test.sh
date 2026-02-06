#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

CONFIG_FILE="${REPO_ROOT}/config/picornavirus.env"
LEGACY_CONFIG="${REPO_ROOT}/config.env"
if [[ -f "${CONFIG_FILE}" ]]; then
  source "${CONFIG_FILE}"
elif [[ -f "${LEGACY_CONFIG}" ]]; then
  source "${LEGACY_CONFIG}"
fi

source "${SCRIPT_DIR}/lib/common.sh"

DB="${DB:-ptv}"
DB_QUERY="${DB_QUERY:-\"Teschovirus\"[Organism]}"
REF_FASTA="${REF_FASTA:-data/ref/${DB}.fa}"
BLAST_DB="${BLAST_DB:-blastdb/${DB}}"
BOWTIE2_INDEX="${BOWTIE2_INDEX:-bowtie2/${DB}}"

echo "== Smoke test do pipeline =="

log_info "Preparando DB (DB=${DB})..."
make -C "$REPO_ROOT" db DB="$DB" DB_QUERY="$DB_QUERY" REF_FASTA="$REF_FASTA" \
  BLAST_DB="$BLAST_DB" BOWTIE2_INDEX="$BOWTIE2_INDEX"

check_file "$REF_FASTA"
log_info "FASTA presente: $REF_FASTA"

missing_db=0
for ext in nhr nin nsq; do
  f="${BLAST_DB}.${ext}"
  if [[ ! -f "$f" ]]; then
    log_warn "arquivo BLAST DB faltando: $f"
    missing_db=1
  fi
done
if [[ $missing_db -ne 0 ]]; then
  log_error "Sugestão: make db DB=${DB}"
fi
log_info "Banco BLAST encontrado em prefixo: $BLAST_DB"

missing_bt2=0
for f in "${BOWTIE2_INDEX}".*.bt2*; do
  if [[ ! -e "$f" ]]; then
    missing_bt2=1
  fi
done
if [[ $missing_bt2 -ne 0 ]]; then
  log_error "índice Bowtie2 ausente para prefixo ${BOWTIE2_INDEX}. Sugestão: make db DB=${DB}"
fi
log_info "Índice Bowtie2 encontrado em prefixo: $BOWTIE2_INDEX"

echo "[Teste rápido] blastn (entrada curta via stdin, sem esperar hits)..."
printf \">q1\\nACTGACTGACTG\\n\" | blastn -query - -db "$BLAST_DB" -outfmt 6 >/dev/null
log_info "blastn executou."

log_info "Smoke test concluído com sucesso."
