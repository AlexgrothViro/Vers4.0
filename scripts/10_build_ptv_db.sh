#!/usr/bin/env bash
set -euo pipefail

DB_DIR="db"
FASTA="${DB_DIR}/ptv_teschovirus.fasta"
DB_PREFIX="${DB_DIR}/ptv_teschovirus"

# Query focada em Teschovirus A (Porcine Teschovirus) com genoma completo
QUERY='"Teschovirus A"[Organism] AND "complete genome"[Title]'

echo "========================================="
echo "  Construção do banco BLAST de PTV"
echo "  (Porcine Teschovirus / Teschovirus A)"
echo "========================================="
echo
echo "Diretório do banco: ${DB_DIR}"
echo "Arquivo FASTA:      ${FASTA}"
echo "Prefixo do DB:      ${DB_PREFIX}"
echo
echo "Query NCBI:"
echo "  ${QUERY}"
echo

mkdir -p "$DB_DIR"

echo "[$(date)] Buscando sequências no NCBI..."
esearch -db nucleotide -query "$QUERY" | efetch -format fasta > "$FASTA"

if [[ ! -s "$FASTA" ]]; then
  echo "ERRO: arquivo FASTA ${FASTA} está vazio. Verifique a query, suas credenciais NCBI ou a conexão." >&2
  exit 1
fi

echo "[$(date)] FASTA salvo em: ${FASTA}"
echo "Pré-visualização das primeiras linhas:"
head -n 5 "$FASTA" || true
echo

echo "[$(date)] Construindo banco BLAST com makeblastdb..."
makeblastdb -in "$FASTA" -dbtype nucl -out "$DB_PREFIX"

echo
echo "[$(date)] Banco BLAST construído com sucesso."
echo "Arquivos gerados em: ${DB_DIR}"
echo
echo "Para usar este banco no pipeline, execute:"
echo "  export BLAST_DB=\"$(pwd)/${DB_PREFIX}\""
