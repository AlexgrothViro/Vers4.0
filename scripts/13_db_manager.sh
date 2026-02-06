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

usage() {
  cat <<'USAGE'
Uso: scripts/13_db_manager.sh <comando>

Comandos:
  list   lista perfis básicos suportados
  setup  baixa FASTA + gera BLAST DB + índice Bowtie2

Variáveis (env):
  DB             (padrão: ptv)
  DB_QUERY       (padrão: "\"Teschovirus\"[Organism]")
  REF_FASTA      (padrão: data/ref/<DB>.fa)
  BLAST_DB       (padrão: blastdb/<DB>)
  BOWTIE2_INDEX  (padrão: bowtie2/<DB>)
USAGE
}

CMD="${1:-}"
case "$CMD" in
  list)
    printf "DB\tQuery\n"
    printf "ptv\t\"Teschovirus\"[Organism]\n"
    exit 0
    ;;
  setup)
    ;;
  -h|--help|"")
    usage
    exit 0
    ;;
  *)
    usage
    exit 1
    ;;
esac

DB="${DB:-ptv}"
DB_QUERY="${DB_QUERY:-\"Teschovirus\"[Organism]}"
REF_FASTA="$(resolve_path "${REF_FASTA:-data/ref/${DB}.fa}")"
BLAST_DB="$(resolve_path "${BLAST_DB:-blastdb/${DB}}")"
BOWTIE2_INDEX="$(resolve_path "${BOWTIE2_INDEX:-bowtie2/${DB}}")"

mkdir -p "$(dirname "$REF_FASTA")" "$(dirname "$BLAST_DB")" "$(dirname "$BOWTIE2_INDEX")"

if [[ ! -s "$REF_FASTA" ]]; then
  command -v esearch >/dev/null 2>&1 || log_error "EDirect não encontrado (esearch). Instale EDirect."
  command -v efetch  >/dev/null 2>&1 || log_error "EDirect não encontrado (efetch). Instale EDirect."
  log_info "Baixando FASTA (DB_QUERY=${DB_QUERY})..."
  esearch -db nucleotide -query "$DB_QUERY" | efetch -format fasta > "$REF_FASTA"
  if [[ ! -s "$REF_FASTA" ]]; then
    log_error "Download falhou, FASTA vazio: $REF_FASTA"
  fi
else
  log_info "FASTA já existe: $REF_FASTA"
fi

blast_missing=0
for ext in nhr nin nsq; do
  file="${BLAST_DB}.${ext}"
  if [[ ! -s "$file" ]]; then
    blast_missing=1
  elif [[ "$file" -ot "$REF_FASTA" ]]; then
    blast_missing=1
  fi
done

if [[ $blast_missing -eq 1 ]]; then
  command -v makeblastdb >/dev/null 2>&1 || log_error "makeblastdb não encontrado (blast+)."
  log_info "Gerando BLAST DB em $BLAST_DB"
  makeblastdb -in "$REF_FASTA" -dbtype nucl -out "$BLAST_DB"
else
  log_info "BLAST DB atualizado: $BLAST_DB"
fi

bt2_missing=0
bt2_files=(
  "${BOWTIE2_INDEX}.1.bt2"
  "${BOWTIE2_INDEX}.2.bt2"
  "${BOWTIE2_INDEX}.3.bt2"
  "${BOWTIE2_INDEX}.4.bt2"
  "${BOWTIE2_INDEX}.rev.1.bt2"
  "${BOWTIE2_INDEX}.rev.2.bt2"
)
for file in "${bt2_files[@]}"; do
  if [[ ! -s "$file" ]]; then
    bt2_missing=1
  elif [[ "$file" -ot "$REF_FASTA" ]]; then
    bt2_missing=1
  fi
done

if [[ $bt2_missing -eq 1 ]]; then
  command -v bowtie2-build >/dev/null 2>&1 || log_error "bowtie2-build não encontrado."
  log_info "Gerando índice Bowtie2 em $BOWTIE2_INDEX"
  bowtie2-build "$REF_FASTA" "$BOWTIE2_INDEX"
else
  log_info "Índice Bowtie2 atualizado: $BOWTIE2_INDEX"
fi
