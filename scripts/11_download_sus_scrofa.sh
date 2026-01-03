#!/usr/bin/env bash
set -euo pipefail

DEST_DIR="ref/host"
FA="${DEST_DIR}/sus_scrofa.fa"

# Query para pegar sequências genômicas de Sus scrofa em RefSeq
QUERY='"Sus scrofa"[Organism] AND srcdb_refseq[PROP] AND biomol_genomic[PROP]'

echo "============================================="
echo "  Download do genoma de Sus scrofa (hospedeiro)"
echo "  Fonte: NCBI nuccore via EDirect"
echo "============================================="
echo "Destino: ${FA}"
echo
echo "Query NCBI:"
echo "  ${QUERY}"
echo

mkdir -p "$DEST_DIR"

echo "[$(date)] Buscando sequências no NCBI..."
esearch -db nucleotide -query "$QUERY" | efetch -format fasta > "$FA"

if [[ ! -s "$FA" ]]; then
  echo "ERRO: arquivo FASTA ${FA} está vazio. Verifique a query ou sua conexão." >&2
  exit 1
fi

echo "[$(date)] Genoma de Sus scrofa salvo em: ${FA}"
echo "Pré-visualização das primeiras linhas:"
head -n 5 "$FA" || true
