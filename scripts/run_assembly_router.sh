#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

if [[ -f "${REPO_ROOT}/config.env" ]]; then
  source "${REPO_ROOT}/config.env"
fi

source "${SCRIPT_DIR}/lib/common.sh"

if [[ ! -f "${REPO_ROOT}/config.env" ]]; then
  log_error "config.env não existe. Crie com: cp config.env.example config.env"
fi

RAW_DIR="$(resolve_path "${RAW_DIR:-data/raw}")"
ASSEMBLY_DIR="$(resolve_path "${ASSEMBLY_DIR:-data/assemblies}")"
RAW1="${RAW_DIR}/${SAMPLE_NAME}_R1.fastq.gz"
RAW2="${RAW_DIR}/${SAMPLE_NAME}_R2.fastq.gz"

check_file "$RAW1"
check_file "$RAW2"

log_info "[PIPELINE] Sample: $SAMPLE_NAME"
log_info "[PIPELINE] Assembler: $ASSEMBLER"

# Diretório padrão “estável” para o BLAST genérico
STD_ASM_DIR="${ASSEMBLY_DIR}/${SAMPLE_NAME}_assembly"
mkdir -p "$STD_ASM_DIR"

if [[ "$ASSEMBLER" == "spades" ]]; then
  "${SCRIPT_DIR}/01_run_spades.sh" "$SAMPLE_NAME" "$THREADS" "$SPADES_PARAMS"
  # SPAdes padrão: contigs.fasta
  ln -sf "${ASSEMBLY_DIR}/${SAMPLE_NAME}_spades/contigs.fasta" "$STD_ASM_DIR/contigs.fa"

elif [[ "$ASSEMBLER" == "velvet" ]]; then
  "${SCRIPT_DIR}/01_run_velvet.sh" "$SAMPLE_NAME" "$VELVET_K" "${VELVET_OPTS:-}"
  ln -sf "${ASSEMBLY_DIR}/${SAMPLE_NAME}_velvet_k${VELVET_K}/contigs.fa" "$STD_ASM_DIR/contigs.fa"

else
  log_error "ASSEMBLER='$ASSEMBLER' inválido. Use 'velvet' ou 'spades'."
fi

log_info "[PIPELINE] OK: contigs padronizados em: $STD_ASM_DIR/contigs.fa"
