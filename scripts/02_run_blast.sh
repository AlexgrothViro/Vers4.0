#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Uso: $0 NOME_AMOSTRA KMER" >&2
  echo "Exemplo: $0 81554_S150 31" >&2
  exit 1
fi

SAMPLE="$1"
KMER="$2"

ASSEMBLY_DIR="data/assemblies"
RESULTS_BLAST_DIR="results/blast"

if [[ -z "${BLAST_DB:-}" ]]; then
  echo "ERRO: variável de ambiente BLAST_DB não definida." >&2
  echo "Defina, por exemplo:" >&2
  echo "  export BLAST_DB=/caminho/para/seu_banco/seu_banco" >&2
  echo "    (sem a extensão .nin/.nhr/.nsq)" >&2
  exit 1
fi

CONTIGS="${ASSEMBLY_DIR}/${SAMPLE}_velvet_k${KMER}/contigs.fa"
OUT_TSV="${RESULTS_BLAST_DIR}/${SAMPLE}_vs_db.tsv"

if [[ ! -f "$CONTIGS" ]]; then
  echo "ERRO: contigs de entrada não encontrados em:" >&2
  echo "  $CONTIGS" >&2
  exit 1
fi

mkdir -p "$RESULTS_BLAST_DIR"

echo "[$(date)] Rodando blastn contra ${BLAST_DB}..."
blastn \
  -query "$CONTIGS" \
  -db "$BLAST_DB" \
  -out "$OUT_TSV" \
  -outfmt "6 qseqid sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -max_target_seqs 5 \
  -num_threads 4

echo "[$(date)] Resultado salvo em: $OUT_TSV"
