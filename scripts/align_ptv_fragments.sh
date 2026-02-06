#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

CONFIG_FILE="${REPO_ROOT}/config/picornavirus.env"
LEGACY_CONFIG="${REPO_ROOT}/config.env"
if [[ -f "${CONFIG_FILE}" ]]; then
  source "${CONFIG_FILE}"
elif [[ -f "${LEGACY_CONFIG}" ]]; then
  source "${LEGACY_CONFIG}"
fi

source "${SCRIPT_DIR}/lib/common.sh"

# Arquivos de entrada/saída
REF_FA="data/ptv_db.fa"
FRAG_FA="run_T1/work/ptv_hits_fragments.fa"
MERGED_FA="run_T1/work/ptv_fragments_plus_ref.fa"
ALN_FA="run_T1/work/ptv_fragments_plus_ref.aln.fa"

echo "=== Alinhamento PTV: referências + fragmentos ==="

# checagens básicas
check_file "$REF_FA"
check_file "$FRAG_FA"

command -v mafft >/dev/null 2>&1 || log_error "'mafft' não encontrado no PATH. Instale o MAFFT antes de rodar este script."

log_info "[1/2] Concatenando referências e fragmentos em: $MERGED_FA"
cat "$REF_FA" "$FRAG_FA" > "$MERGED_FA"

log_info "[2/2] Rodando MAFFT (modo automático)..."
mafft --auto "$MERGED_FA" > "$ALN_FA"

echo
log_info "[OK] Alinhamento pronto em: $ALN_FA"
