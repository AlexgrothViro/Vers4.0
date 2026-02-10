# Painel de uso (UX) — Picornavirus-quali2026

O painel UX é um dashboard simples para usuários não técnicos rodarem o pipeline no WSL/Ubuntu
com acesso via navegador do Windows.

## 1. Como iniciar

No WSL, dentro do repositório:

```bash
make fix-wsl
make deps
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

## 1.1. Checklist rápido de instalação

1. **WSL/Ubuntu instalado** (Windows).
2. **Repositório clonado** no WSL.
3. **Dependências instaladas** via `make deps`.
4. **Painel rodando** via `make dashboard`.

Se alguma etapa falhar, veja “Observações importantes” ao final deste documento.

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

## 4. Dicas de uso no Windows (WSL)

### 4.1. Caminhos de arquivos FASTQ

Para arquivos que estão no Windows, use caminhos do WSL, por exemplo:

```
/mnt/c/Users/SEU_USUARIO/Downloads/minha_amostra_R1.fastq.gz
```

### 4.2. Se o `localhost` não abrir

Use o IP retornado por:

```bash
hostname -I
```

E acesse no navegador: `http://<IP_DO_WSL>:8787`.

### 4.3. Se `make deps` falhar

Geralmente é bloqueio de rede/proxy. Ajuste o proxy do `apt-get` ou instale
as dependências manualmente e rode `make test-env` para validar.
