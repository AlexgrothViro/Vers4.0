# Painel de uso (UX) — Picornavirus-quali2026

O painel UX é um dashboard simples para usuários não técnicos rodarem o pipeline no WSL/Ubuntu
com acesso via navegador do Windows.

## 1. Como iniciar

No WSL, dentro do repositório:

```bash
make dashboard
```

O painel ficará acessível em:

- `http://localhost:8787` (Windows, mesmo PC)
- ou `http://<IP_DO_WSL>:8787` (se o `localhost` não funcionar)

Para descobrir o IP do WSL:

```bash
hostname -I
```

> O terminal deve permanecer aberto enquanto o painel estiver em uso.

## 2. O que o painel faz

O painel oferece quatro ações principais:

1. **Verificar ambiente** — roda `make test-env` para checar dependências.
2. **Gerar DEMO** — roda `make demo` e cria FASTQs em `data/raw`.
3. **Importar amostra** — executa `scripts/00_import_sample.sh` para copiar/registrar R1/R2.
4. **Rodar pipeline** — executa `scripts/20_run_pipeline.sh` com assembler configurável.

## 3. Observações importantes

- **SPAdes** só é usado quando `config.env` estiver presente e `ASSEMBLER=spades`.
- Logs ficam em `logs/ux_dashboard_*.log` (ignorados pelo git).
- Caso precise parar o painel, use `Ctrl+C` no terminal.
