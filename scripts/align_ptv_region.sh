#!/usr/bin/env bash
set -euo pipefail

IN_FA="${1:-run_T1/work/ptv_region_KX686489_1_7209.fa}"
OUT_FA="${2:-run_T1/work/ptv_region_KX686489_1_7209.aln.fa}"

echo "=== Alinhando região PTV com MAFFT ==="
echo "[INFO] FASTA de entrada: $IN_FA"
echo "[INFO] Saída de alinhamento: $OUT_FA"

if [[ ! -s "$IN_FA" ]]; then
  echo "ERRO: arquivo de entrada não encontrado ou vazio: $IN_FA" >&2
  exit 1
fi

if ! command -v mafft >/dev/null 2>&1; then
  echo "ERRO: 'mafft' não encontrado no PATH. Instale o MAFFT antes de rodar este script." >&2
  exit 1
fi

mafft --auto "$IN_FA" > "$OUT_FA"

echo "[OK] Alinhamento escrito em: $OUT_FA"
