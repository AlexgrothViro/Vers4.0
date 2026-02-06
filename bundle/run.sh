#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MICRO="$ROOT/.bundle/bin/micromamba"
ENV_DIR="$ROOT/.bundle/env"

if [[ ! -x "$MICRO" || ! -d "$ENV_DIR" ]]; then
  echo "[INFO] Ambiente não instalado. Rodando installer..."
  bash "$ROOT/bundle/install_wsl.sh"
fi

TARGET="${1:-help}"
shift || true

# Roda Makefile dentro do env, garantindo PATH correto
exec "$MICRO" run -p "$ENV_DIR" make -C "$ROOT" "$TARGET" "$@"
