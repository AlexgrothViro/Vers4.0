#!/usr/bin/env bash
set -euo pipefail

# Carrega config
if [[ ! -f config.env ]]; then
  echo "[ERRO] config.env não existe. Crie com: cp config.env.example config.env"
  exit 1
fi
source config.env

RAW1="data/raw/${SAMPLE_NAME}_R1.fastq.gz"
RAW2="data/raw/${SAMPLE_NAME}_R2.fastq.gz"

if [[ ! -s "$RAW1" || ! -s "$RAW2" ]]; then
  echo "[ERRO] FASTQ não encontrados:"
  echo "  $RAW1"
  echo "  $RAW2"
  exit 1
fi

echo "[PIPELINE] Sample: $SAMPLE_NAME"
echo "[PIPELINE] Assembler: $ASSEMBLER"

# Diretório padrão “estável” para o BLAST genérico
STD_ASM_DIR="data/assemblies/${SAMPLE_NAME}_assembly"
mkdir -p "$STD_ASM_DIR"

if [[ "$ASSEMBLER" == "spades" ]]; then
  ./scripts/01_run_spades.sh "$SAMPLE_NAME" "$THREADS" "$SPADES_PARAMS"
  # SPAdes padrão: contigs.fasta
  ln -sf "../${SAMPLE_NAME}_spades/contigs.fasta" "$STD_ASM_DIR/contigs.fa"

elif [[ "$ASSEMBLER" == "velvet" ]]; then
  ./scripts/01_run_velvet.sh "$SAMPLE_NAME" "$VELVET_K" "${VELVET_OPTS:-}"
  ln -sf "../${SAMPLE_NAME}_velvet_k${VELVET_K}/contigs.fa" "$STD_ASM_DIR/contigs.fa"

else
  echo "[ERRO] ASSEMBLER='$ASSEMBLER' inválido. Use 'velvet' ou 'spades'."
  exit 1
fi

echo "[PIPELINE] OK: contigs padronizados em: $STD_ASM_DIR/contigs.fa"
