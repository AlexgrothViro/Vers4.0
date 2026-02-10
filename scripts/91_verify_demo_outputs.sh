#!/usr/bin/env bash
set -euo pipefail

sample="${1:-DEMO}"
blast_file="results/blast/${sample}_vs_db.tsv"
report_file="results/reports/${sample}_summary.md"

[[ -s "data/assemblies/${sample}_velvet_k31/contigs.fa" ]] || { echo "[ERRO] contigs ausente"; exit 1; }
[[ -s "$blast_file" ]] || { echo "[ERRO] blast ausente: $blast_file"; exit 1; }
[[ -s "$report_file" ]] || { echo "[ERRO] report ausente: $report_file"; exit 1; }

grep -q "adj_identity" "$report_file" || { echo "[ERRO] coluna adj_identity ausente no report"; exit 1; }
grep -q "|qseqid|sseqid|pident|length|qlen|aln_cov|adj_identity|evalue|bitscore|" "$report_file" || {
  echo "[ERRO] tabela esperada ausente no report"
  exit 1
}

echo "[OK] artefatos DEMO validados"
