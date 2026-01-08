#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

if [[ -f "${REPO_ROOT}/config.env" ]]; then
  source "${REPO_ROOT}/config.env"
fi

source "${SCRIPT_DIR}/lib/common.sh"

if [[ $# -lt 1 ]]; then
  log_error "Uso: $0 SAMPLE [KMER]"
fi

SAMPLE="$1"
KMER="${2:-31}"

ROOT_DIR="$(pwd)"
BLAST_DB="${ROOT_DIR}/db/ptv_teschovirus"
WORK_DIR="${ROOT_DIR}/run_T1/work"
ASSEMBLY_DIR="${ROOT_DIR}/data/assemblies/${SAMPLE}_velvet_k${KMER}"

mkdir -p "$WORK_DIR"

export BLAST_DB

log_info "======================================="
log_info " Pipeline avançado PTV - SAMPLE=${SAMPLE} KMER=${KMER}"
log_info "======================================="

log_info "[1/7] make test (env, filtro hospedeiro, montagem, BLAST básico)..."
make test SAMPLE="${SAMPLE}" KMER="${KMER}"

log_info "[2/7] BLAST de confirmação com qseq/sseq..."
blastn \
  -query "${ASSEMBLY_DIR}/contigs.fa" \
  -db "${BLAST_DB}" \
  -out "${WORK_DIR}/ptv_hits.confirm.tsv" \
  -outfmt "6 qseqid sacc pident length mismatch gapopen evalue bitscore qstart qend sstart send qseq sseq" \
  -max_target_seqs 5 \
  -num_threads 4

log_info "[3/7] Identidade ajustada (adj_identity.py)..."
python3 scripts/adj_identity.py \
  "${WORK_DIR}/ptv_hits.confirm.tsv" \
  "${WORK_DIR}/ptv_hits.adjust.tsv"

log_info "[4/7] Relatório resumido (merge_report.py)..."
python3 scripts/merge_report.py

log_info "[5/7] Plano de extensão de flancos (extend_plan.py)..."
python3 scripts/extend_plan.py

log_info "[6/7] FASTA de flancos (emit_extend_fasta.py)..."
python3 scripts/emit_extend_fasta.py

log_info "[7/7] Simulação de reads (sim_reads_clean.py) e rótulos (label_hits.py)..."
python3 scripts/sim_reads_clean.py
python3 scripts/label_hits.py > "${WORK_DIR}/ptv_hits.labels.tsv"

echo
log_info "Pipeline avançado concluído."
log_info "Arquivos principais em: ${WORK_DIR}"
echo "  - ptv_hits.confirm.tsv"
echo "  - ptv_hits.adjust.tsv"
echo "  - ptv_report.tsv"
echo "  - extend_plan.tsv"
echo "  - extend_regions.fasta"
echo "  - sim_R1.fastq.gz / sim_R2.fastq.gz"
echo "  - ptv_hits.labels.tsv"
