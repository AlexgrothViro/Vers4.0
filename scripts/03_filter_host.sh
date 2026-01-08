#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

if [[ -f "${REPO_ROOT}/config.env" ]]; then
  source "${REPO_ROOT}/config.env"
fi

source "${SCRIPT_DIR}/lib/common.sh"

if [[ $# -lt 1 ]]; then
  log_error "Uso: $0 NOME_AMOSTRA"
fi

SAMPLE="$1"

RAW_DIR="$(resolve_path "${RAW_DIR:-data/raw}")"
HOST_REMOVED_DIR="$(resolve_path "${HOST_REMOVED_DIR:-data/host_removed}")"
HOST_INDEX_PREFIX="$(resolve_path "${HOST_INDEX_PREFIX:-ref/host/sus_scrofa_bt2}")"

R1="${RAW_DIR}/${SAMPLE}_R1.fastq.gz"
R2="${RAW_DIR}/${SAMPLE}_R2.fastq.gz"

check_file "$R1"
check_file "$R2"

# Verifica se o índice do Bowtie2 existe
check_file "${HOST_INDEX_PREFIX}.1.bt2"

mkdir -p "$HOST_REMOVED_DIR"

TMP_PREFIX="${HOST_REMOVED_DIR}/${SAMPLE}_host_removed.tmp"

log_info "Filtrando leituras do hospedeiro (Sus scrofa) para amostra ${SAMPLE}..."

bowtie2 \
  -x "$HOST_INDEX_PREFIX" \
  -1 "$R1" \
  -2 "$R2" \
  --very-sensitive \
  -p 4 \
  --un-conc-gz "${TMP_PREFIX}.fastq.gz" \
  -S /dev/null

# Bowtie2 com --un-conc-gz gera dois arquivos:
#   ${TMP_PREFIX}.fastq.1.gz  e  ${TMP_PREFIX}.fastq.2.gz
mv "${TMP_PREFIX}.fastq.1.gz" "${HOST_REMOVED_DIR}/${SAMPLE}_R1.host_removed.fastq.gz"
mv "${TMP_PREFIX}.fastq.2.gz" "${HOST_REMOVED_DIR}/${SAMPLE}_R2.host_removed.fastq.gz"

log_info "Leituras não alinhadas ao hospedeiro salvas em:"
echo "  ${HOST_REMOVED_DIR}/${SAMPLE}_R1.host_removed.fastq.gz"
echo "  ${HOST_REMOVED_DIR}/${SAMPLE}_R2.host_removed.fastq.gz"
