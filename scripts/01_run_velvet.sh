#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

usage() {
  cat <<'USAGE' >&2
Uso:
  bash scripts/01_run_velvet.sh SAMPLE [KMER]

Regras:
- SAMPLE passado por argumento SEMPRE vence.
- Procura reads em: data/host_removed -> data/cleaned -> data/raw
USAGE
}

SAMPLE_NAME="${1:-}"
KMER="${2:-${VELVET_K:-31}}"
THREADS="${THREADS:-4}"
MIN_CONTIG_LEN="${MIN_CONTIG_LEN:-100}"

[[ -z "$SAMPLE_NAME" ]] && usage && exit 1

pick_reads() {
  local s="$1"
  for d in "data/host_removed" "data/cleaned" "data/raw"; do
    local r1="${REPO_ROOT}/${d}/${s}_R1.fastq.gz"
    local r2="${REPO_ROOT}/${d}/${s}_R2.fastq.gz"
    if [[ -s "$r1" && -s "$r2" ]]; then
      echo "$r1|$r2|$d"
      return 0
    fi
  done
  return 1
}

if ! picked="$(pick_reads "$SAMPLE_NAME")"; then
  echo "[ERRO] FASTQs não encontrados para SAMPLE='$SAMPLE_NAME'." >&2
  echo "       Procurei em: data/host_removed, data/cleaned, data/raw" >&2
  echo "[DICA] Amostras em data/raw:" >&2
  ls -1 "${REPO_ROOT}/data/raw"/*_R1.fastq.gz 2>/dev/null | sed -E 's#.*/##; s/_R1\.fastq\.gz$//' | sort -u | sed 's/^/  - /' >&2 || true
  exit 1
fi

R1="${picked%%|*}"; rest="${picked#*|}"
R2="${rest%%|*}"; SRC_DIR="${rest#*|}"

OUTDIR="${REPO_ROOT}/data/assemblies/${SAMPLE_NAME}_velvet_k${KMER}"
mkdir -p "$OUTDIR"

echo "[INFO] Velvet: SAMPLE=$SAMPLE_NAME KMER=$KMER"
echo "[INFO] Reads: $SRC_DIR"
echo "[INFO] Out: $OUTDIR"

velveth "$OUTDIR" "$KMER" -shortPaired -fastq.gz -separate "$R1" "$R2"
velvetg "$OUTDIR" -exp_cov auto -cov_cutoff auto -read_trkg yes -min_contig_lgth "$MIN_CONTIG_LEN"

[[ -s "$OUTDIR/contigs.fa" ]] || { echo "[ERRO] Velvet não gerou contigs.fa em $OUTDIR" >&2; exit 1; }
echo "[OK] Contigs: $OUTDIR/contigs.fa"
