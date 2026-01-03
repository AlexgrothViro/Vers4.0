#!/usr/bin/env bash
set -euo pipefail

# Programas que o pipeline precisa pelo menos para os testes básicos
REQUIRED_CMDS=(velveth velvetg blastn)

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

echo
if [[ -n "${BLAST_DB:-}" ]]; then
  echo "BLAST_DB está definido como: $BLAST_DB"
  echo "Certifique-se de que os arquivos de índice do BLAST existem nesse caminho."
else
  echo "ATENÇÃO: variável de ambiente BLAST_DB não definida."
  echo "O alvo 'test-blast' do Makefile vai precisar dela."
fi

if [[ $MISSING -ne 0 ]]; then
  echo
  echo "Alguns programas estão faltando. Veja a seção 'Ambiente e requisitos' no README.md." >&2
  exit 1
fi

echo
echo "Ambiente básico OK."
