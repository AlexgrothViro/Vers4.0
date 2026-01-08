#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

if [[ -f "${REPO_ROOT}/config.env" ]]; then
  source "${REPO_ROOT}/config.env"
fi

source "${SCRIPT_DIR}/lib/common.sh"

SAMPLE="${1:?SAMPLE obrigatório}"
THREADS="${2:-4}"
SPADES_PARAMS="${3:-}"

RAW_DIR="$(resolve_path "${RAW_DIR:-data/raw}")"
ASSEMBLY_DIR="$(resolve_path "${ASSEMBLY_DIR:-data/assemblies}")"
RAW1="${RAW_DIR}/${SAMPLE}_R1.fastq.gz"
RAW2="${RAW_DIR}/${SAMPLE}_R2.fastq.gz"
OUTDIR="${ASSEMBLY_DIR}/${SAMPLE}_spades"

command -v spades.py >/dev/null 2>&1 || log_error "spades.py não encontrado no PATH"

check_file "$RAW1"
check_file "$RAW2"

mkdir -p "$OUTDIR"
log_info "[SPAdes] sample=$SAMPLE threads=$THREADS params='$SPADES_PARAMS'"
spades.py \
  -1 "$RAW1" -2 "$RAW2" \
  -o "$OUTDIR" \
  -t "$THREADS" \
  $SPADES_PARAMS

check_file "$OUTDIR/contigs.fasta"
log_info "[SPAdes] OK: $OUTDIR/contigs.fasta"
