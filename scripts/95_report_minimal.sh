#!/usr/bin/env bash
set -euo pipefail

SAMPLE=""; CONTIGS=""; BLAST=""; OUT=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="${2:-}"; shift 2 ;;
    --contigs) CONTIGS="${2:-}"; shift 2 ;;
    --blast) BLAST="${2:-}"; shift 2 ;;
    --out) OUT="${2:-}"; shift 2 ;;
    *) echo "[ERRO] arg inválido: $1" >&2; exit 1 ;;
  esac
done
[[ -z "$SAMPLE" || -z "$CONTIGS" || -z "$BLAST" || -z "$OUT" ]] && { echo "[ERRO] uso: --sample --contigs --blast --out" >&2; exit 1; }

mkdir -p "$(dirname "$OUT")"
contig_count="$(grep -c '^>' "$CONTIGS" 2>/dev/null || echo 0)"
max_len="$(awk '
  /^>/ { if (seqlen>max) max=seqlen; seqlen=0; next }
  { seqlen += length($0) }
  END { if (seqlen>max) max=seqlen; print (max+0) }
' "$CONTIGS" 2>/dev/null || echo 0)"
hit_count="$(wc -l < "$BLAST" 2>/dev/null || echo 0)"

{
  echo "# Summary – $SAMPLE"
  echo
  echo "## Assembly"
  echo "- Contigs: ${contig_count}"
  echo "- Maior contig (bp): ${max_len}"
  echo
  echo "## BLAST (vs PTV DB)"
  echo "- Hits (linhas TSV): ${hit_count}"
  echo
  echo "### Top 5 hits (por bitscore)"
  if [[ -s "$BLAST" ]]; then
    echo
    echo "|qseqid|sseqid|pident|length|evalue|bitscore|"
    echo "|---|---:|---:|---:|---:|---:|"
    sort -k12,12gr "$BLAST" | head -n 5 | awk -F'\t' '{printf "|%s|%s|%s|%s|%s|%s|\n",$1,$2,$3,$4,$11,$12}'
  else
    echo "_Sem hits (arquivo BLAST vazio)_"
  fi
} > "$OUT"
echo "[OK] Relatório: $OUT"
