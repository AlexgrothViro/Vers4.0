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

SAMPLE="${1:?SAMPLE obrigatório}"
THREADS="${2:-4}"
SPADES_PARAMS="${3:-}"

RAW_DIR="$(resolve_path "${RAW_DIR:-data/raw}")"
ASSEMBLY_DIR="$(resolve_path "${ASSEMBLY_DIR:-data/assemblies}")"
OUTDIR="${ASSEMBLY_DIR}/${SAMPLE}_spades"

pick_first() { printf '%s\n' "$@" 2>/dev/null | sed '/^$/d' | sort | head -n1; }
find_read() {
  local sample="$1" raw_dir="$2" which="$3"
  local -a cands=()
  shopt -s nullglob
  if [[ "$which" == "R1" ]]; then
    cands=(
      "$raw_dir/${sample}_R1"*.fastq.gz "$raw_dir/${sample}_R1"*.fq.gz
      "$raw_dir/${sample}_R1"*.fastq    "$raw_dir/${sample}_R1"*.fq
      "$raw_dir/${sample}"*"_R1"*".fastq.gz" "$raw_dir/${sample}"*"_R1"*".fq.gz"
      "$raw_dir/${sample}"*"_R1"*".fastq"    "$raw_dir/${sample}"*"_R1"*".fq"
      "$raw_dir/${sample}_1"*.fastq.gz "$raw_dir/${sample}_1"*.fq.gz
      "$raw_dir/${sample}_1"*.fastq    "$raw_dir/${sample}_1"*.fq
    )
  else
    cands=(
      "$raw_dir/${sample}_R2"*.fastq.gz "$raw_dir/${sample}_R2"*.fq.gz
      "$raw_dir/${sample}_R2"*.fastq    "$raw_dir/${sample}_R2"*.fq
      "$raw_dir/${sample}"*"_R2"*".fastq.gz" "$raw_dir/${sample}"*"_R2"*".fq.gz"
      "$raw_dir/${sample}"*"_R2"*".fastq"    "$raw_dir/${sample}"*"_R2"*".fq"
      "$raw_dir/${sample}_2"*.fastq.gz "$raw_dir/${sample}_2"*.fq.gz
      "$raw_dir/${sample}_2"*.fastq    "$raw_dir/${sample}_2"*.fq
    )
  fi
  shopt -u nullglob
  pick_first "${cands[@]}"
}
# Prioridade: env R1/R2 -> auto-detect em RAW_DIR
RAW1="${R1:-${READ1:-${FASTQ_R1:-}}}"
RAW2="${R2:-${READ2:-${FASTQ_R2:-}}}"

if [[ -z "${RAW1}" || -z "${RAW2}" ]]; then
  RAW1="$(find_read "$SAMPLE" "$RAW_DIR" "R1")"
  RAW2="$(find_read "$SAMPLE" "$RAW_DIR" "R2")"
fi

if [[ -z "${RAW1}" || -z "${RAW2}" ]]; then
  log_error "FASTQs não encontrados para sample='$SAMPLE'.
Procurei em: $RAW_DIR
Exemplos aceitos:
  ${SAMPLE}_R1*.fastq(.gz) / ${SAMPLE}_R2*.fastq(.gz)
  ${SAMPLE}_1*.fq(.gz)     / ${SAMPLE}_2*.fq(.gz)
Alternativa (recomendado):
  R1=/caminho/read1.fastq.gz R2=/caminho/read2.fastq.gz make pipeline SAMPLE=$SAMPLE"
fi

command -v spades.py >/dev/null 2>&1 || log_error "spades.py não encontrado no PATH"

check_file "$RAW1"
check_file "$RAW2"

mkdir -p "$OUTDIR"
log_info "[SPAdes] sample=$SAMPLE"
log_info "[SPAdes] R1=$RAW1"
log_info "[SPAdes] R2=$RAW2"
log_info "[SPAdes] threads=$THREADS params='$SPADES_PARAMS'"
spades.py \
  -1 "$RAW1" -2 "$RAW2" \
  -o "$OUTDIR" \
  -t "$THREADS" \
  $SPADES_PARAMS

check_file "$OUTDIR/contigs.fasta"
log_info "[SPAdes] OK: $OUTDIR/contigs.fasta"
