# Painel de Uso (UX) — Picornavirus-quali2026

Este painel é uma interface web simples para usuários não técnicos executarem as principais
etapas do pipeline no WSL/Ubuntu e acessar pelo navegador do Windows.

## Como iniciar

Dentro do WSL, na raiz do repositório:

```bash
make ux
```

Por padrão o painel usa a porta **8000**. Para acessar no Windows, abra o navegador e visite:

```
http://localhost:8000
```

> **Importante:** mantenha o terminal do WSL aberto enquanto o painel estiver em execução.

## Funcionalidades

- **Verificar ambiente**: roda `make test-env` e confirma dependências instaladas.
- **Gerar DEMO**: cria FASTQs de demonstração em `data/raw`.
- **Importar amostras**: cria links ou cópias de R1/R2 para `data/raw`.
- **Rodar pipeline**: executa a pipeline com Velvet (padrão) ou SPAdes (se `config.env` existir).
- **Status visual de término**: mostra SUCESSO/FALHA com link para report ou tail do log (30 linhas).
- **Histórico**: aba “Histórico” lista runs persistidos em `results/runs/*/run.json` com botões para abrir report/log/blast e reexecutar.

## Observações importantes

- O painel executa comandos locais; não envia dados para a internet.
- Logs e resultados permanecem nos diretórios padrão (`data/`, `results/`, `logs/`).
- Cada execução cria `logs/ux_dashboard_*.log` (log completo) e, para pipeline, snapshot em `results/runs/<timestamp>_<sample>/`.
- Para SPAdes, assegure que `config.env` está configurado conforme `config.env.example`.
- Caminhos de arquivos do Windows podem ser informados usando o formato `C:\\...`.

## Encerrar

Use `Ctrl + C` no terminal onde o painel está rodando para interromper o servidor.
