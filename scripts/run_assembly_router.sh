#!/usr/bin/env bash
set -euo pipefail

# Carrega as configurações
source config.env

echo "[PIPELINE] Iniciando montagem para: $SAMPLE_NAME"
echo "[PIPELINE] Montador escolhido: $ASSEMBLER"

if [[ "$ASSEMBLER" == "spades" ]]; then
    # Chama o SPAdes usando os parâmetros do config
    scripts/01_run_spades.sh "$SAMPLE_NAME" "$THREADS" "$SPADES_PARAMS"

elif [[ "$ASSEMBLER" == "velvet" ]]; then
    # Chama o Velvet usando os parâmetros do config
    scripts/01_run_velvet.sh "$SAMPLE_NAME" "$THREADS" "$VELVET_K" "$VELVET_OPTS"

else
    echo "ERRO: Montador '$ASSEMBLER' não reconhecido no config.env"
    exit 1
fi
