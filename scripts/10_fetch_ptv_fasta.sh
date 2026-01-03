#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Uso: $0 <saida_fasta> (ex.: data/ref/ptv_db.fa)" >&2
  exit 1
fi

OUT_FASTA="$1"
OUT_DIR="$(dirname "$OUT_FASTA")"

if [[ -s "$OUT_FASTA" ]]; then
  echo "[OK] FASTA já existe e não está vazio: $OUT_FASTA"
  exit 0
fi

mkdir -p "$OUT_DIR"

QUERY='"Teschovirus"[Organism]'

have_edirect() {
  command -v esearch >/dev/null 2>&1 && command -v efetch >/dev/null 2>&1
}

if ! have_edirect; then
  cat >&2 <<'EOF'
ERRO: EDirect (esearch/efetch) não encontrado no PATH.
Instale EDirect (https://www.ncbi.nlm.nih.gov/books/NBK179288/)
ou gere o FASTA manualmente e salve em data/ref/ptv_db.fa.
EOF
  exit 1
fi

echo "Baixando sequências de Teschovirus do NCBI (QUERY=${QUERY})..."
esearch -db nucleotide -query "$QUERY" | efetch -format fasta > "$OUT_FASTA"

if [[ ! -s "$OUT_FASTA" ]]; then
  echo "ERRO: download falhou, $OUT_FASTA está vazio." >&2
  exit 1
fi

echo "[OK] FASTA salvo em $OUT_FASTA"
