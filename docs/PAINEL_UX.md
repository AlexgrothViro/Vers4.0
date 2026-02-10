# Painel de Uso (UX) — Picornavirus-quali2026

## Como iniciar

No WSL:

```bash
make ux
```

No Windows: `http://localhost:8000`.

Ou use diretamente no Windows:

```bat
bundle\start_platform.bat
```

## Funcionalidades atuais

- **Check de ambiente** (`make test-env`).
- **DEMO reprodutível** (`make demo`).
- **DB viral genérico (NCBI)**:
  - dropdown de alvo (`config/targets.json`),
  - query customizada opcional,
  - taxid opcional,
  - execução do script `scripts/10_build_viral_db.sh`.
- **Importação de amostras**:
  - por caminho colado (R1/R2) com normalização Windows→WSL,
  - por upload de R1/R2,
  - por upload de `.zip` (auto-detecção R1/R2),
  - padronização em `data/raw/<sample>_R1.fastq.gz` e `_R2.fastq.gz`.
- **Execução do pipeline por sample_id** (sem digitar path completo).
- **Histórico e logs** com artefatos de execução.

## Validações e mensagens amigáveis

- `R2 não encontrado`
- `arquivo vazio`
- `NCBI query retornou 0 sequências`
- validações de extensão: `.fastq` / `.fastq.gz`

## Logs

- Dashboard gera logs timestampados em `logs/ux_dashboard_*.log`.
- Execuções de pipeline são persistidas em `results/runs/<timestamp>_<sample>/run.json`.
