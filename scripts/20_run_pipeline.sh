#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

CONFIG_FILE="${REPO_ROOT}/config/picornavirus.env"
LEGACY_CONFIG="${REPO_ROOT}/config.env"
HAS_CONFIG=0
if [[ -f "${CONFIG_FILE}" ]]; then
  source "${CONFIG_FILE}"
  HAS_CONFIG=1
elif [[ -f "${LEGACY_CONFIG}" ]]; then
  source "${LEGACY_CONFIG}"
  HAS_CONFIG=1
fi

usage() {
  cat <<'USAGE'
Uso: scripts/20_run_pipeline.sh [opções]

Opções:
  --install           instala dependências via apt-get (usa 00_check_env.sh)
  --sample NOME       nome da amostra (padrão: SAMPLE_ID/SAMPLE_NAME ou amostra_teste)
  --kmer K            k-mer para Velvet (padrão: VELVET_K ou 31)
  --skip-host-filter  ignora o filtro do hospedeiro
  -h, --help          mostra esta ajuda

Obs.: se existir config/picornavirus.env (ou config.env legado), ele será usado como base de configuração.
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

SAMPLE_ID="${SAMPLE_OVERRIDE:-${SAMPLE_ID:-${SAMPLE_NAME:-amostra_teste}}}"
SAMPLE_NAME="${SAMPLE_ID}"
VELVET_K="${KMER_OVERRIDE:-${VELVET_K:-31}}"
ASSEMBLER="${ASSEMBLER:-velvet}"
THREADS="${THREADS:-4}"
DB="${DB:-ptv}"
BLAST_DB="$(resolve_path "${BLAST_DB:-blastdb/${DB}}")"
RAW_DIR="$(resolve_path "${RAW_DIR:-data/raw}")"

SAMPLE_R1="${SAMPLE_R1:-}"
SAMPLE_R2="${SAMPLE_R2:-}"
SAMPLE_SINGLE="${SAMPLE_SINGLE:-}"

if [[ -n "$SAMPLE_SINGLE" && ( -n "$SAMPLE_R1" || -n "$SAMPLE_R2" ) ]]; then
  echo "[ERRO] Use --single ou --r1/--r2, mas não ambos." >&2
  exit 1
fi

if [[ -n "$SAMPLE_SINGLE" ]]; then
  SAMPLE_SINGLE="$(resolve_path "$SAMPLE_SINGLE")"
else
  if [[ -n "$SAMPLE_R1" ]]; then
    SAMPLE_R1="$(resolve_path "$SAMPLE_R1")"
  else
    SAMPLE_R1="${RAW_DIR}/${SAMPLE_ID}_R1.fastq.gz"
  fi

  if [[ -n "$SAMPLE_R2" ]]; then
    SAMPLE_R2="$(resolve_path "$SAMPLE_R2")"
  else
    SAMPLE_R2="${RAW_DIR}/${SAMPLE_ID}_R2.fastq.gz"
  fi
fi

export SAMPLE_ID SAMPLE_NAME SAMPLE_R1 SAMPLE_R2 SAMPLE_SINGLE RAW_DIR

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
make -C "$REPO_ROOT" db DB="$DB" DB_QUERY="${DB_QUERY:-}"

MISSING_READS=0
if [[ -n "$SAMPLE_SINGLE" ]]; then
  if [[ ! -f "$SAMPLE_SINGLE" ]]; then
    echo "[ERRO] Read single não encontrado: $SAMPLE_SINGLE" >&2
    MISSING_READS=1
  fi
else
  if [[ ! -f "$SAMPLE_R1" ]]; then
    echo "[ERRO] Read R1 não encontrado: $SAMPLE_R1" >&2
    MISSING_READS=1
  fi
  if [[ ! -f "$SAMPLE_R2" ]]; then
    echo "[ERRO] Read R2 não encontrado: $SAMPLE_R2" >&2
    MISSING_READS=1
  fi
fi

if [[ $MISSING_READS -ne 0 ]]; then
  echo "Sugestão: make run ID=${SAMPLE_ID} R1=/caminho/R1.fastq.gz R2=/caminho/R2.fastq.gz DB=${DB}" >&2
  exit 1
fi

if [[ $SKIP_HOST_FILTER -eq 0 ]]; then
  log "[3/6] Filtrando hospedeiro (opcional)"
  if [[ -n "$SAMPLE_SINGLE" ]]; then
    echo "[AVISO] Filtro do hospedeiro ignora amostras single-end. Use --skip-host-filter para silenciar." >&2
  else
  HOST_INDEX_PREFIX="${REPO_ROOT}/ref/host/sus_scrofa_bt2"
  if [[ -f "${HOST_INDEX_PREFIX}.1.bt2" ]]; then
    "$SCRIPT_DIR/03_filter_host.sh" "$SAMPLE_NAME"
  else
    echo "[AVISO] Índice do hospedeiro não encontrado em ${HOST_INDEX_PREFIX}.1.bt2"
    echo "        Se quiser rodar o filtro de hospedeiro, execute scripts/11_download_sus_scrofa.sh"
    echo "        e gere o índice Bowtie2 manualmente."
  fi
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
