#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Uso:
  bash scripts/00_import_sample.sh --sample NOME --r1 CAMINHO --r2 CAMINHO [--copy]

Faz symlink (padrão) ou cópia para:
  data/raw/NOME_R1.fastq.gz
  data/raw/NOME_R2.fastq.gz
EOF
}

SAMPLE=""
R1=""
R2=""
MODE="link"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="${2:-}"; shift 2;;
    --r1) R1="${2:-}"; shift 2;;
    --r2) R2="${2:-}"; shift 2;;
    --copy) MODE="copy"; shift 1;;
    -h|--help) usage; exit 0;;
    *) echo "[ERRO] argumento inválido: $1"; usage; exit 2;;
  esac
done

[[ -n "$SAMPLE" && -n "$R1" && -n "$R2" ]] || { echo "[ERRO] faltou --sample/--r1/--r2"; usage; exit 2; }

# Windows -> WSL se necessário
if command -v wslpath >/dev/null 2>&1; then
  [[ "$R1" =~ ^[A-Za-z]:\\ ]] && R1="$(wslpath -u "$R1")"
  [[ "$R2" =~ ^[A-Za-z]:\\ ]] && R2="$(wslpath -u "$R2")"
fi

mkdir -p data/raw

out1="data/raw/${SAMPLE}_R1.fastq.gz"
out2="data/raw/${SAMPLE}_R2.fastq.gz"

for f in "$R1" "$R2"; do
  [[ -s "$f" ]] || { echo "[ERRO] arquivo não existe ou vazio: $f"; exit 1; }
done

do_place() {
  local src="$1" dst="$2"
  if [[ "$MODE" == "copy" ]]; then
    cp -f "$src" "$dst"
  else
    ln -sf "$(realpath "$src")" "$dst"
  fi
}

do_place "$R1" "$out1"
do_place "$R2" "$out2"

echo "[OK] Amostra importada:"
echo "  $out1 -> $R1"
echo "  $out2 -> $R2"
