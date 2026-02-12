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

if [[ ! -f "${CONFIG_FILE}" && ! -f "${LEGACY_CONFIG}" ]]; then
  log_error "config/picornavirus.env não existe. Crie com: cp config/picornavirus.env.example config/picornavirus.env"
fi

RAW_DIR="$(resolve_path "${RAW_DIR:-data/raw}")"
ASSEMBLY_DIR="$(resolve_path "${ASSEMBLY_DIR:-data/assemblies}")"

SAMPLE_KEY="${SAMPLE_ID:-$SAMPLE_NAME}"
if [[ -n "${SAMPLE_SINGLE:-}" ]]; then
  log_error "Montagem requer leituras pareadas. Use SAMPLE_R1/SAMPLE_R2."
fi

if [[ -n "${SAMPLE_R1:-}" ]]; then
  RAW1="$(resolve_path "${SAMPLE_R1}")"
else
  RAW1="${RAW_DIR}/${SAMPLE_KEY}_R1.fastq.gz"
fi

if [[ -n "${SAMPLE_R2:-}" ]]; then
  RAW2="$(resolve_path "${SAMPLE_R2}")"
else
  RAW2="${RAW_DIR}/${SAMPLE_KEY}_R2.fastq.gz"
fi

check_file "$RAW1"
check_file "$RAW2"

log_info "[PIPELINE] Sample: $SAMPLE_NAME"
log_info "[PIPELINE] Assembler: $ASSEMBLER"

# Diretório padrão “estável” para o BLAST genérico
STD_ASM_DIR="${ASSEMBLY_DIR}/${SAMPLE_KEY}_assembly"
mkdir -p "$STD_ASM_DIR"

if [[ "$ASSEMBLER" == "spades" ]]; then
  "${SCRIPT_DIR}/01_run_spades.sh" "$SAMPLE_KEY" "$THREADS" "$SPADES_PARAMS"
  # SPAdes padrão: contigs.fasta
  ln -sf "${ASSEMBLY_DIR}/${SAMPLE_KEY}_spades/contigs.fasta" "$STD_ASM_DIR/contigs.fa"

elif [[ "$ASSEMBLER" == "velvet" ]]; then
  "${SCRIPT_DIR}/01_run_velvet.sh" "$SAMPLE_KEY" "$VELVET_K" "${VELVET_OPTS:-}"
  ln -sf "${ASSEMBLY_DIR}/${SAMPLE_KEY}_velvet_k${VELVET_K}/contigs.fa" "$STD_ASM_DIR/contigs.fa"

else
  log_error "ASSEMBLER='$ASSEMBLER' inválido. Use 'velvet' ou 'spades'."
fi

log_info "[PIPELINE] OK: contigs padronizados em: $STD_ASM_DIR/contigs.fa"
