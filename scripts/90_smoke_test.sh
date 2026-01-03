#!/usr/bin/env bash
set -euo pipefail

PTV_FASTA="data/ptv_db.fa"
BLAST_DB="${BLAST_DB:-blastdb/ptv}"
BOWTIE2_INDEX="${BOWTIE2_INDEX:-bowtie2/ptv}"

echo "== Smoke test do pipeline =="

if [[ ! -s "$PTV_FASTA" ]]; then
  echo "ERRO: FASTA de referência ausente ou vazio: $PTV_FASTA" >&2
  exit 1
fi
echo "[OK] FASTA presente: $PTV_FASTA"

missing_db=0
for ext in nhr nin nsq; do
  f="${BLAST_DB}.${ext}"
  if [[ ! -f "$f" ]]; then
    echo "ERRO: arquivo BLAST DB faltando: $f" >&2
    missing_db=1
  fi
done
if [[ $missing_db -ne 0 ]]; then
  echo "Sugestão: make blastdb" >&2
  exit 1
fi
echo "[OK] Banco BLAST encontrado em prefixo: $BLAST_DB"

missing_bt2=0
for f in "${BOWTIE2_INDEX}".*.bt2*; do
  if [[ ! -e "$f" ]]; then
    missing_bt2=1
  fi
done
if [[ $missing_bt2 -ne 0 ]]; then
  echo "ERRO: índice Bowtie2 ausente para prefixo ${BOWTIE2_INDEX}" >&2
  echo "Sugestão: make bowtie2-index" >&2
  exit 1
fi
echo "[OK] Índice Bowtie2 encontrado em prefixo: $BOWTIE2_INDEX"

echo "[Teste rápido] blastn (entrada curta via stdin, sem esperar hits)..."
printf \">q1\\nACTGACTGACTG\\n\" | blastn -query - -db "$BLAST_DB" -outfmt 6 >/dev/null
echo "[OK] blastn executou."

echo "Smoke test concluído com sucesso."
