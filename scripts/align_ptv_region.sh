#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

if [[ -f "${REPO_ROOT}/config.env" ]]; then
  source "${REPO_ROOT}/config.env"
fi

source "${SCRIPT_DIR}/lib/common.sh"

IN_FA="${1:-run_T1/work/ptv_region_KX686489_1_7209.fa}"
OUT_FA="${2:-run_T1/work/ptv_region_KX686489_1_7209.aln.fa}"

log_info "=== Alinhando região PTV com MAFFT ==="
log_info "FASTA de entrada: $IN_FA"
log_info "Saída de alinhamento: $OUT_FA"

check_file "$IN_FA"

command -v mafft >/dev/null 2>&1 || log_error "'mafft' não encontrado no PATH. Instale o MAFFT antes de rodar este script."

mafft --auto "$IN_FA" > "$OUT_FA"

log_info "Alinhamento escrito em: $OUT_FA"
