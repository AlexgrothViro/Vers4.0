#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

if [[ -f "${REPO_ROOT}/config.env" ]]; then
  source "${REPO_ROOT}/config.env"
fi

source "${SCRIPT_DIR}/lib/common.sh"

PTV_FASTA="data/ptv_db.fa"
BLAST_DB="${BLAST_DB:-blastdb/ptv}"
BOWTIE2_INDEX="${BOWTIE2_INDEX:-bowtie2/ptv}"

echo "== Smoke test do pipeline =="

check_file "$PTV_FASTA"
log_info "FASTA presente: $PTV_FASTA"

missing_db=0
for ext in nhr nin nsq; do
  f="${BLAST_DB}.${ext}"
  if [[ ! -f "$f" ]]; then
    log_warn "arquivo BLAST DB faltando: $f"
    missing_db=1
  fi
done
if [[ $missing_db -ne 0 ]]; then
  log_error "Sugestão: make blastdb"
fi
log_info "Banco BLAST encontrado em prefixo: $BLAST_DB"

missing_bt2=0
for f in "${BOWTIE2_INDEX}".*.bt2*; do
  if [[ ! -e "$f" ]]; then
    missing_bt2=1
  fi
done
if [[ $missing_bt2 -ne 0 ]]; then
  log_error "índice Bowtie2 ausente para prefixo ${BOWTIE2_INDEX}. Sugestão: make bowtie2-index"
fi
log_info "Índice Bowtie2 encontrado em prefixo: $BOWTIE2_INDEX"

echo "[Teste rápido] blastn (entrada curta via stdin, sem esperar hits)..."
printf \">q1\\nACTGACTGACTG\\n\" | blastn -query - -db "$BLAST_DB" -outfmt 6 >/dev/null
log_info "blastn executou."

log_info "Smoke test concluído com sucesso."
