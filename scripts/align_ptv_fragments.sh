#!/usr/bin/env bash
set -euo pipefail

# Arquivos de entrada/saída
REF_FA="data/ptv_db.fa"
FRAG_FA="run_T1/work/ptv_hits_fragments.fa"
MERGED_FA="run_T1/work/ptv_fragments_plus_ref.fa"
ALN_FA="run_T1/work/ptv_fragments_plus_ref.aln.fa"

echo "=== Alinhamento PTV: referências + fragmentos ==="

# checagens básicas
if [[ ! -s "$REF_FA" ]]; then
  echo "ERRO: Arquivo de referência não encontrado ou vazio: $REF_FA" >&2
  exit 1
fi

if [[ ! -s "$FRAG_FA" ]]; then
  echo "ERRO: Arquivo de fragmentos não encontrado ou vazio: $FRAG_FA" >&2
  exit 1
fi

if ! command -v mafft >/dev/null 2>&1; then
  echo "ERRO: 'mafft' não encontrado no PATH. Instale o MAFFT antes de rodar este script." >&2
  exit 1
fi

echo "[1/2] Concatenando referências e fragmentos em: $MERGED_FA"
cat "$REF_FA" "$FRAG_FA" > "$MERGED_FA"

echo "[2/2] Rodando MAFFT (modo automático)..."
mafft --auto "$MERGED_FA" > "$ALN_FA"

echo
echo "[OK] Alinhamento pronto em: $ALN_FA"
