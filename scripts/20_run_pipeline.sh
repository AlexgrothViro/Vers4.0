#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

if [[ -f "${REPO_ROOT}/config.env" ]]; then
  source "${REPO_ROOT}/config.env"
  HAS_CONFIG=1
else
  HAS_CONFIG=0
fi

usage() {
  cat <<'USAGE'
Uso: scripts/20_run_pipeline.sh [opções]

Opções:
  --install           instala dependências via apt-get (usa 00_check_env.sh)
  --sample NOME       nome da amostra (padrão: SAMPLE_NAME ou 81554_S150)
  --kmer K            k-mer para Velvet (padrão: VELVET_K ou 31)
  --skip-host-filter  ignora o filtro do hospedeiro
  -h, --help          mostra esta ajuda

Obs.: se existir config.env, ele será usado como base de configuração.
USAGE
}

resolve_path() {
  local path="$1"
  if [[ "$path" == /* ]]; then
    echo "$path"
  else
    echo "${REPO_ROOT}/${path}"
  fi
}

AUTO_INSTALL=0
SKIP_HOST_FILTER=0
SAMPLE_OVERRIDE=""
KMER_OVERRIDE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --install)
      AUTO_INSTALL=1
      shift
      ;;
    --sample)
      SAMPLE_OVERRIDE="${2:-}"
      shift 2
      ;;
    --kmer)
      KMER_OVERRIDE="${2:-}"
      shift 2
      ;;
    --skip-host-filter)
      SKIP_HOST_FILTER=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      usage
      exit 1
      ;;
  esac
done

SAMPLE_NAME="${SAMPLE_OVERRIDE:-${SAMPLE_NAME:-${SAMPLE:-81554_S150}}}"
VELVET_K="${KMER_OVERRIDE:-${VELVET_K:-31}}"
ASSEMBLER="${ASSEMBLER:-velvet}"
THREADS="${THREADS:-4}"
BLAST_DB="$(resolve_path "${BLAST_DB:-blastdb/ptv}")"

log() {
  printf '\n== %s ==\n' "$1"
}

log "[1/6] Verificando ambiente"
if [[ $AUTO_INSTALL -eq 1 ]]; then
  "$SCRIPT_DIR/00_check_env.sh" --install
else
  "$SCRIPT_DIR/00_check_env.sh"
fi

log "[2/6] Preparando diretórios e bancos"
make -C "$REPO_ROOT" setup_dirs
make -C "$REPO_ROOT" ptv-fasta
make -C "$REPO_ROOT" blastdb
make -C "$REPO_ROOT" bowtie2-index

if [[ $SKIP_HOST_FILTER -eq 0 ]]; then
  log "[3/6] Filtrando hospedeiro (opcional)"
  HOST_INDEX_PREFIX="${REPO_ROOT}/ref/host/sus_scrofa_bt2"
  if [[ -f "${HOST_INDEX_PREFIX}.1.bt2" ]]; then
    "$SCRIPT_DIR/03_filter_host.sh" "$SAMPLE_NAME"
  else
    echo "[AVISO] Índice do hospedeiro não encontrado em ${HOST_INDEX_PREFIX}.1.bt2"
    echo "        Se quiser rodar o filtro de hospedeiro, execute scripts/11_download_sus_scrofa.sh"
    echo "        e gere o índice Bowtie2 manualmente."
  fi
fi

log "[4/6] Montagem de contigs"
ASSEMBLY_CONTIGS=""
if [[ $HAS_CONFIG -eq 1 && "$ASSEMBLER" == "spades" ]]; then
  "$SCRIPT_DIR/run_assembly_router.sh"
  ASSEMBLY_CONTIGS="${REPO_ROOT}/data/assemblies/${SAMPLE_NAME}_assembly/contigs.fa"
else
  "$SCRIPT_DIR/01_run_velvet.sh" "$SAMPLE_NAME" "$VELVET_K"
  ASSEMBLY_CONTIGS="${REPO_ROOT}/data/assemblies/${SAMPLE_NAME}_velvet_k${VELVET_K}/contigs.fa"
fi

log "[5/6] BLAST dos contigs"
if [[ ! -f "${ASSEMBLY_CONTIGS}" ]]; then
  echo "[ERRO] Contigs não encontrados em ${ASSEMBLY_CONTIGS}" >&2
  exit 1
fi
if [[ ! -f "${BLAST_DB}.nhr" ]]; then
  echo "[ERRO] Banco BLAST não encontrado em ${BLAST_DB}.nhr" >&2
  exit 1
fi

OUTDIR="${REPO_ROOT}/results/blast"
mkdir -p "$OUTDIR"

if [[ "$ASSEMBLER" == "velvet" ]]; then
  OUT="${OUTDIR}/${SAMPLE_NAME}_k${VELVET_K}_vs_db.tsv"
else
  OUT="${OUTDIR}/${SAMPLE_NAME}_assembly_vs_db.tsv"
fi

blastn -query "$ASSEMBLY_CONTIGS" -db "$BLAST_DB" \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
  -max_target_seqs 5 -evalue 1e-5 -num_threads "$THREADS" > "$OUT"

ln -sf "$(basename "$OUT")" "${OUTDIR}/${SAMPLE_NAME}_vs_db.tsv"

echo "Resultado salvo em: $OUT"

log "[6/6] Pipeline concluído"
echo "Amostra: $SAMPLE_NAME"
echo "Contigs: $ASSEMBLY_CONTIGS"
