#!/usr/bin/env bash
set -euo pipefail

ALN_FA="${1:-run_T1/work/ptv_region_KX686489_1_7209.aln.fa}"
PREFIX="${2:-run_T1/work/ptv_region_KX686489_1_7209.iqtree}"

echo "=== Árvore PTV (região KX686489.1:1-7192) com IQ-TREE ==="

if [[ ! -s "$ALN_FA" ]]; then
  echo "ERRO: alinhamento não encontrado ou vazio: $ALN_FA" >&2
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
