#!/usr/bin/env bash
set -euo pipefail

ALN_FA="run_T1/work/ptv_fragments_plus_ref.aln.fa"
PREFIX="run_T1/work/ptv_fragments_plus_ref.iqtree"

echo "=== Árvore PTV com IQ-TREE ==="

if [[ ! -s "$ALN_FA" ]]; then
  echo "ERRO: Alinhamento não encontrado ou vazio: $ALN_FA" >&2
  exit 1
fi

IQBIN=""

if command -v iqtree2 >/dev/null 2>&1; then
  IQBIN="iqtree2"
elif command -v iqtree >/dev/null 2>&1; then
  IQBIN="iqtree"
else
  echo "ERRO: IQ-TREE não encontrado (nem 'iqtree2' nem 'iqtree' no PATH)." >&2
  echo "Instale, por exemplo, com (Ubuntu):" >&2
  echo "  sudo apt update && sudo apt install iqtree" >&2
  echo "ou baixe o binário oficial em iqtree.org e coloque no PATH." >&2
  exit 1
fi

echo "[INFO] Usando binário: $IQBIN"

"$IQBIN" -s "$ALN_FA" \
  -m GTR+G \
  -nt AUTO \
  -bb 1000 \
  -alrt 1000 \
  -pre "$PREFIX"

echo
echo "[OK] IQ-TREE finalizado."
echo "Arquivos principais em:"
echo "  ${PREFIX}.treefile  (árvore principal)"
echo "  ${PREFIX}.log       (log da análise)"
echo "  ${PREFIX}.iqtree    (resumo do modelo e estatísticas)"
