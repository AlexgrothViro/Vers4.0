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

ADJ_TSV="$(mktemp)"
python3 "$(dirname "$0")/adj_identity.py" --blast "$BLAST" --contigs "$CONTIGS" --out "$ADJ_TSV"
top_adj_summary=""
if [[ -s "$ADJ_TSV" ]]; then
  top_adj_summary="$(tail -n +2 "$ADJ_TSV" | sort -t$'\t' -k9,9gr | head -n 1 | awk -F'\t' '{printf "%s vs %s | adj_identity=%s%% | cobertura=%s%%",$1,$2,$9,$8}')"
fi

{
  echo "# Summary – $SAMPLE"
  echo
  echo "## Assembly"
  echo "- Contigs: ${contig_count}"
  echo "- Maior contig (bp): ${max_len}"
  echo
  echo "## BLAST (vs PTV DB)"
  echo "- Hits (linhas TSV): ${hit_count}"
  if [[ -n "$top_adj_summary" ]]; then
    echo "- Melhor identidade ajustada: ${top_adj_summary}"
  fi
  echo
  echo "### Top 5 hits (por adj_identity)"
  if [[ -s "$BLAST" ]]; then
    echo
    echo "|qseqid|sseqid|pident|length|qlen|aln_cov|adj_identity|evalue|bitscore|"
    echo "|---|---:|---:|---:|---:|---:|---:|---:|---:|"
    tail -n +2 "$ADJ_TSV" | sort -t$'\t' -k9,9gr | head -n 5 | awk -F'\t' '{printf "|%s|%s|%s|%s|%s|%s|%s|%s|%s|\n",$1,$2,$3,$4,$7,$8,$9,$5,$6}'
  else
    echo "_Sem hits (arquivo BLAST vazio)_"
  fi
} > "$OUT"

rm -f "$ADJ_TSV"
echo "[OK] Relatório: $OUT"
