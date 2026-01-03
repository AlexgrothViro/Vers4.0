#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Uso: $0 NOME_AMOSTRA" >&2
  echo "Exemplo: $0 81554_S150" >&2
  exit 1
fi

SAMPLE="$1"

RAW_DIR="data/raw"
HOST_REMOVED_DIR="data/host_removed"
HOST_INDEX_PREFIX="ref/host/sus_scrofa_bt2"

R1="${RAW_DIR}/${SAMPLE}_R1.fastq.gz"
R2="${RAW_DIR}/${SAMPLE}_R2.fastq.gz"

if [[ ! -f "$R1" || ! -f "$R2" ]]; then
  echo "ERRO: FASTQ de entrada não encontrados:" >&2
  echo "  $R1" >&2
  echo "  $R2" >&2
  exit 1
fi

# Verifica se o índice do Bowtie2 existe
if [[ ! -f "${HOST_INDEX_PREFIX}.1.bt2" ]]; then
  echo "ERRO: índice Bowtie2 do hospedeiro não encontrado em ${HOST_INDEX_PREFIX}.*.bt2" >&2
  echo "Crie com:" >&2
  echo "  bowtie2-build ref/host/sus_scrofa.fa ref/host/sus_scrofa_bt2" >&2
  exit 1
fi

mkdir -p "$HOST_REMOVED_DIR"

TMP_PREFIX="${HOST_REMOVED_DIR}/${SAMPLE}_host_removed.tmp"

echo "[$(date)] Filtrando leituras do hospedeiro (Sus scrofa) para amostra ${SAMPLE}..."

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

echo "[$(date)] Leituras não alinhadas ao hospedeiro salvas em:"
echo "  ${HOST_REMOVED_DIR}/${SAMPLE}_R1.host_removed.fastq.gz"
echo "  ${HOST_REMOVED_DIR}/${SAMPLE}_R2.host_removed.fastq.gz"
