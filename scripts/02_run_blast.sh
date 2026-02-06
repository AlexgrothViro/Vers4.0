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
KMER="${2:?KMER obrigatório}"

DB="$(resolve_path "${BLAST_DB:-blastdb/ptv}")"
THREADS="${BLAST_THREADS:-4}"

CONTIGS="$(resolve_path "data/assemblies/${SAMPLE}_velvet_k${KMER}/contigs.fa")"
OUTDIR="$(resolve_path "results/blast")"
OUT="${OUTDIR}/${SAMPLE}_k${KMER}_vs_db.tsv"

mkdir -p "$OUTDIR"

check_file "$CONTIGS"
check_file "${DB}.nhr"

log_info "Rodando blastn contra $DB..."
blastn -query "$CONTIGS" -db "$DB" \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
  -max_target_seqs 5 -evalue 1e-5 -num_threads "$THREADS" > "$OUT"

# compat legado (alguns scripts antigos podem ler o nome sem kmer)
ln -sf "$(basename "$OUT")" "${OUTDIR}/${SAMPLE}_vs_db.tsv"

log_info "Resultado salvo em: $OUT"
