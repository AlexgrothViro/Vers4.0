#!/usr/bin/env bash
set -euo pipefail

SAMPLE="${1:-${SAMPLE:-81554_S150}}"
KMER_PTV="${2:-${KMER_PTV:-51}}"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

cd "$ROOT_DIR"

echo "[quickstart] Verificando dependencias..."
./scripts/00_check_env.sh

echo "[quickstart] Preparando diretórios e bancos locais..."
make setup_dirs
make ptv-fasta
make blastdb
make bowtie2-index

echo "[quickstart] Rodando pipeline PTV enriquecido..."
./scripts/05_ptv_enriched_run.sh "$SAMPLE" "$KMER_PTV"

echo "[quickstart] Concluído. Resultados em results/ e data/."
