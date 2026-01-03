#!/usr/bin/env bash
set -euo pipefail

# Programas obrigatórios para os testes básicos
REQUIRED_CMDS=(
  velveth
  velvetg
  blastn
  makeblastdb
  bowtie2
  bowtie2-build
  python3
)

# Programas opcionais (úteis para download via NCBI/HTTP)
OPTIONAL_CMDS=(esearch efetch curl)

echo "== Verificando programas necessários no PATH =="

MISSING=0
for cmd in "${REQUIRED_CMDS[@]}"; do
  if command -v "$cmd" >/dev/null 2>&1; then
    path=$(command -v "$cmd")
    printf "  [OK] %s encontrado em %s\n" "$cmd" "$path"
  else
    echo "  [FALTA] $cmd não está no PATH"
    MISSING=1
  fi
done

for cmd in "${OPTIONAL_CMDS[@]}"; do
  if command -v "$cmd" >/dev/null 2>&1; then
    path=$(command -v "$cmd")
    printf "  [OK] %s encontrado em %s (opcional)\n" "$cmd" "$path"
  else
    echo "  [OPCIONAL] $cmd não está no PATH (usado em downloads NCBI/HTTP)"
  fi
done

echo
PTV_FASTA="data/ptv_db.fa"
if [[ ! -s "$PTV_FASTA" ]]; then
  echo "ATENÇÃO: FASTA de referência ausente em $PTV_FASTA."
  echo "  Sugestão: make ptv-fasta-legacy   # cria data/ptv_db.fa a partir de data/ref/ptv_db.fa"
else
  echo "FASTA de referência encontrado: $PTV_FASTA"
fi

if ! command -v esearch >/dev/null 2>&1 || ! command -v efetch >/dev/null 2>&1; then
  echo "AVISO: EDirect ausente (esearch/efetch) - necessário para baixar FASTA do NCBI."
fi

if [[ -n "${BLAST_DB:-}" ]]; then
  echo "BLAST_DB está definido como: $BLAST_DB"
  missing_db_files=()
  for ext in nhr nin nsq; do
    file="${BLAST_DB}.${ext}"
    [[ -f "$file" ]] || missing_db_files+=("$file")
  done
  if (( ${#missing_db_files[@]} > 0 )); then
    echo "  [AVISO] Índices BLAST ausentes:"
    printf '    - %s\n' "${missing_db_files[@]}"
    echo "  Sugestão: make blastdb"
  else
    echo "  [OK] Índices BLAST encontrados."
  fi
else
  echo "ATENÇÃO: variável de ambiente BLAST_DB não definida."
  echo "  Sugestão: make blastdb   # gera blastdb/ptv e exporte BLAST_DB=blastdb/ptv"
fi

if [[ $MISSING -ne 0 ]]; then
  echo
  echo "Alguns programas estão faltando. Veja a seção 'Ambiente e requisitos' no README.md." >&2
  exit 1
fi

echo
echo "Ambiente básico OK."
