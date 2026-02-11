#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

CONFIG_FILE="${REPO_ROOT}/config/picornavirus.env"
LEGACY_CONFIG="${REPO_ROOT}/config.env"
if [[ -f "${CONFIG_FILE}" ]]; then
  # shellcheck disable=SC1090
  source "${CONFIG_FILE}"
elif [[ -f "${LEGACY_CONFIG}" ]]; then
  # shellcheck disable=SC1090
  source "${LEGACY_CONFIG}"
fi

# shellcheck disable=SC1091
source "${SCRIPT_DIR}/lib/common.sh"

declare -A DB_QUERIES
declare -A DB_DESC

init_profiles() {
   DB_QUERIES[ptv]='"Teschovirus"[Organism]'
   DB_DESC[ptv]='Porcine teschovirus (Teschovirus)'

   DB_QUERIES[evg]='"Enterovirus G"[Organism]'
   DB_DESC[evg]='Enterovirus G (suínos)'

   DB_QUERIES[psv]='"Sapelovirus A"[Organism]'
   DB_DESC[psv]='Sapelovirus A (porcine sapelovirus)'

   DB_QUERIES[svv]='"Senecavirus A"[Organism]'
   DB_DESC[svv]='Senecavirus A'

   DB_QUERIES[fmdv]='"Foot-and-mouth disease virus"[Organism]'
   DB_DESC[fmdv]='FMDV (aftosa)'

   # ---- Picornaviridae (perfis do UX) ----
   DB_QUERIES[picornaviridae_refseq]='"Picornaviridae"[Organism] AND refseq[filter]'
   DB_DESC[picornaviridae_refseq]='Picornaviridae (RefSeq) [recomendado]'

   DB_QUERIES[picornaviridae_complete]='"Picornaviridae"[Organism] AND ("complete genome"[Title] OR "complete cds"[Title])'
   DB_DESC[picornaviridae_complete]='Picornaviridae (complete genome/cds)'

   DB_QUERIES[picornaviridae_all]='"Picornaviridae"[Organism]'
   DB_DESC[picornaviridae_all]='Picornaviridae (ALL) [gigante]'

   # (opcional) mantém o alias antigo:
   DB_QUERIES[picornaviridae]='"Picornaviridae"[Organism]'
   DB_DESC[picornaviridae]='TODOS Picornaviridae (alias antigo)'
}

usage() {
  cat <<'USAGE'
Uso: scripts/13_db_manager.sh <comando>

Comandos:
  list   lista perfis básicos suportados (DB -> query NCBI)
  setup  baixa FASTA + gera BLAST DB + índice Bowtie2

Variáveis (env):
  DB             (padrão: ptv)
  DB_QUERY       (se definido, sobrescreve a query padrão do perfil)
  NCBI_DB        (padrão: nucleotide)

  REF_FASTA      (padrão: data/ref/<DB>.fa)
  BLAST_DB       (padrão: blastdb/<DB>)
  BOWTIE2_INDEX  (padrão: bowtie2/<DB>)

Exemplos:
  scripts/13_db_manager.sh list
  DB=evg scripts/13_db_manager.sh setup
  DB=custom DB_QUERY='"Porcine kobuvirus"[All Fields]' scripts/13_db_manager.sh setup
USAGE
}

init_profiles

CMD="${1:-}"
case "$CMD" in
  list)
      # Ordem estável (não depende da ordem do associative array)
      DB_ORDER=(
        ptv evg psv svv fmdv
        picornaviridae_refseq picornaviridae_complete picornaviridae_all
      )

      if [[ "${2:-}" == "--json" ]]; then
        printf '[
'
        first=1
        for id in "${DB_ORDER[@]}"; do
          q="${DB_QUERIES[$id]:-}"
          d="${DB_DESC[$id]:-}"
          [[ -n "$q" ]] || continue

          # escape básico pra JSON
          q_esc="${q//\\/\\\\}"
          q_esc="${q_esc//\"/\\\"}"
          d_esc="${d//\\/\\\\}"
          d_esc="${d_esc//\"/\\\"}"

          if [[ $first -eq 0 ]]; then printf ',
'; fi
          first=0
          printf '  {"id":"%s","label":"%s","query":"%s"}' "$id" "$d_esc" "$q_esc"
        done
        printf '
]
'
        exit 0
      fi

      printf "DB	Desc	Query
"
      for id in "${DB_ORDER[@]}"; do
        q="${DB_QUERIES[$id]:-}"
        d="${DB_DESC[$id]:-}"
        [[ -n "$q" ]] || continue
        printf "%s	%s	%s
" "$id" "$d" "$q"
      done
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
NCBI_DB="${NCBI_DB:-nucleotide}"

DEFAULT_QUERY="${DB_QUERIES[$DB]:-}"
DB_QUERY="${DB_QUERY:-$DEFAULT_QUERY}"

if [[ -z "${DB_QUERY}" ]]; then
  log_error "DB_QUERY vazio e DB '${DB}' não tem perfil conhecido. Use DB_QUERY=... ou DB=ptv/evg/psv/svv/fmdv/picornaviridae_refseq/picornaviridae_complete/picornaviridae_all."
fi

REF_FASTA="$(resolve_path "${REF_FASTA:-data/ref/${DB}.fa}")"
BLAST_DB="$(resolve_path "${BLAST_DB:-blastdb/${DB}}")"
BOWTIE2_INDEX="$(resolve_path "${BOWTIE2_INDEX:-bowtie2/${DB}}")"

mkdir -p "$(dirname "$REF_FASTA")" "$(dirname "$BLAST_DB")" "$(dirname "$BOWTIE2_INDEX")"

if [[ ! -s "$REF_FASTA" ]]; then
  command -v esearch >/dev/null 2>&1 || log_error "EDirect não encontrado (esearch). Instale EDirect."
  command -v efetch  >/dev/null 2>&1 || log_error "EDirect não encontrado (efetch). Instale EDirect."
  log_info "Baixando FASTA (NCBI_DB=${NCBI_DB}; DB_QUERY=${DB_QUERY})..."
  esearch -db "$NCBI_DB" -query "$DB_QUERY" | efetch -format fasta > "$REF_FASTA"
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

log_info "OK (DB=${DB})"
