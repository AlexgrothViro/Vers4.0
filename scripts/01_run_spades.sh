#!/usr/bin/env bash
set -euo pipefail

SAMPLE="${1:?SAMPLE obrigatório}"
THREADS="${2:-4}"
SPADES_PARAMS="${3:-}"

RAW1="data/raw/${SAMPLE}_R1.fastq.gz"
RAW2="data/raw/${SAMPLE}_R2.fastq.gz"
OUTDIR="data/assemblies/${SAMPLE}_spades"

command -v spades.py >/dev/null 2>&1 || { echo "[ERRO] spades.py não encontrado no PATH"; exit 1; }

mkdir -p "$OUTDIR"
echo "[SPAdes] sample=$SAMPLE threads=$THREADS params='$SPADES_PARAMS'"
spades.py \
  -1 "$RAW1" -2 "$RAW2" \
  -o "$OUTDIR" \
  -t "$THREADS" \
  $SPADES_PARAMS

[[ -s "$OUTDIR/contigs.fasta" ]] || { echo "[ERRO] contigs.fasta não gerado"; exit 1; }
echo "[SPAdes] OK: $OUTDIR/contigs.fasta"
