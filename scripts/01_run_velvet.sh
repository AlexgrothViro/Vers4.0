#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Uso: $0 NOME_AMOSTRA [KMER]" >&2
  echo "Exemplo: $0 exemplo 31" >&2
  exit 1
fi

SAMPLE="$1"
KMER="${2:-31}"

RAW_DIR="data/raw"
ASSEMBLY_DIR="data/assemblies"

R1="${RAW_DIR}/${SAMPLE}_R1.fastq.gz"
R2="${RAW_DIR}/${SAMPLE}_R2.fastq.gz"
OUTDIR="${ASSEMBLY_DIR}/${SAMPLE}_velvet_k${KMER}"

if [[ ! -f "$R1" || ! -f "$R2" ]]; then
  echo "ERRO: FASTQ de entrada não encontrados:" >&2
  echo "  $R1" >&2
  echo "  $R2" >&2
  echo "Coloque os arquivos de teste em data/raw/ com esse padrão de nome." >&2
  exit 1
fi

mkdir -p "$OUTDIR"

echo "[$(date)] Rodando velveth (k=${KMER})..."
velveth "$OUTDIR" "$KMER" -shortPaired -fastq.gz -separate "$R1" "$R2"

echo "[$(date)] Rodando velvetg..."
velvetg "$OUTDIR" -exp_cov auto -cov_cutoff auto

if [[ -f "${OUTDIR}/contigs.fa" ]]; then
  echo "[$(date)] OK: contigs gerados em ${OUTDIR}/contigs.fa"
else
  echo "[$(date)] ATENÇÃO: contigs.fa não encontrado em ${OUTDIR}" >&2
EXTRA_OPTS="${3:-${VELVET_OPTS:-}}"    
# ...
velvetg "$OUTDIR" $EXTRA_OPTS
  exit 1
fi
