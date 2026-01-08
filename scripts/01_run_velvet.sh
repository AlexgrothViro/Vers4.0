#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

if [[ -f "${REPO_ROOT}/config.env" ]]; then
  source "${REPO_ROOT}/config.env"
fi

source "${SCRIPT_DIR}/lib/common.sh"

if [[ $# -lt 1 ]]; then
  log_error "Uso: $0 NOME_AMOSTRA [KMER]"
fi

SAMPLE="$1"
KMER="${2:-31}"

RAW_DIR="$(resolve_path "${RAW_DIR:-data/raw}")"
ASSEMBLY_DIR="$(resolve_path "${ASSEMBLY_DIR:-data/assemblies}")"

R1="${RAW_DIR}/${SAMPLE}_R1.fastq.gz"
R2="${RAW_DIR}/${SAMPLE}_R2.fastq.gz"
OUTDIR="${ASSEMBLY_DIR}/${SAMPLE}_velvet_k${KMER}"

check_file "$R1"
check_file "$R2"

mkdir -p "$OUTDIR"

log_info "Rodando velveth (k=${KMER})..."
velveth "$OUTDIR" "$KMER" -shortPaired -fastq.gz -separate "$R1" "$R2"

log_info "Rodando velvetg..."
EXTRA_OPTS="${3:-${VELVET_OPTS:-}}"
velvetg "$OUTDIR" -exp_cov auto -cov_cutoff auto $EXTRA_OPTS

if [[ -f "${OUTDIR}/contigs.fa" ]]; then
  log_info "OK: contigs gerados em ${OUTDIR}/contigs.fa"
else
  log_error "ATENÇÃO: contigs.fa não encontrado em ${OUTDIR}"
fi
