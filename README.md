# Pipeline de montagem e anotação preliminar de Porcine Teschovirus (quali-2025)

Neste repositório estou organizando o pipeline que estou desenvolvendo para a qualificação do meu mestrado em Virologia.  
O foco é montar e anotar **Porcine Teschovirus (Teschovirus A)** a partir de dados de sequenciamento (FASTQ pareado), com:

- **remoção de leituras do hospedeiro** (Sus scrofa),
- **montagem de contigs** com Velvet/SPAdes,
- **BLAST** dos contigs contra um **banco específico de Teschovirus**.

A ideia é ter algo **simples de rodar em qualquer máquina Linux/WSL2**, usando `make` e alguns scripts bash, para poder replicar o fluxo em outros computadores e em outras amostras.

---

## 1. Estrutura do repositório

A organização básica é:

```text
quali-2025/
├── Makefile
├── README.md
├── data/
│   ├── raw/             # FASTQ brutos (.fastq.gz)
│   ├── cleaned/         # reservado para futuros pré-processamentos
│   ├── host_removed/    # FASTQ após remoção de hospedeiro
│   └── assemblies/      # contigs e stats das montagens
├── db/
│   └── .gitkeep         # databases locais (PTV) são gerados por script
├── docs/                # documentação, textos, figuras (futuro)
├── results/
│   ├── qc/              # reservado para QC
│   ├── blast/           # resultados de BLAST (TSV)
│   ├── phylogeny/       # arquivos de filogenia (futuro)
│   └── reports/         # relatórios, tabelas finais (futuro)
└── scripts/
    ├── 00_check_env.sh          # checa se os programas principais estão instalados
    ├── 01_run_velvet.sh         # monta contigs com Velvet
    ├── 02_run_blast.sh          # roda BLAST dos contigs contra o banco de PTV
    ├── 03_filter_host.sh        # filtra leituras do hospedeiro (Sus scrofa) com Bowtie2
    ├── 10_build_ptv_db.sh       # baixa sequências de Teschovirus A (NCBI) e gera banco BLAST
    └── 11_download_sus_scrofa.sh# baixa o genoma do hospedeiro (Sus scrofa) via NCBI/EDirect

---

## 2. Instalando dependências automaticamente

Se estiver em um ambiente Debian/Ubuntu (inclui WSL), basta rodar:

```bash
make deps
```

Esse alvo usa `apt-get` para instalar os requisitos do pipeline (build-essential, SPAdes, Velvet, BLAST+, Bowtie2, MAFFT, FastTree, IQ-TREE/IQ-TREE2, dos2unix e EDirect via `ncbi-entrez-direct`), e depois revalida o ambiente com `scripts/00_check_env.sh`. Caso `apt-get` não esteja disponível, o script aborta e mantém as mensagens de requisitos para instalação manual.

> Dica: se a instalação falhar por restrições de rede/proxy, ajuste o proxy de `apt-get` conforme o ambiente antes de rodar `make deps`.

## 3. WSL/Windows (CRLF)

Ambientes Windows/WSL podem clonar o repositório com finais de linha em CRLF ou perder o bit de execução dos scripts. Para evitar erros do tipo `^M` ou `Permission denied`, execute antes de rodar o pipeline:

```bash
make fix-wsl
./scripts/00_check_env.sh --install
```

O alvo `fix-wsl` normaliza os arquivos rastreados para LF e restaura permissões de execução nos scripts. Em seguida, `00_check_env.sh --install` garante que as dependências (incluindo ferramentas de filogenia) estejam disponíveis. Depois disso, siga o fluxo normal (`make ptv-fasta`, `make blastdb`, etc.).

---

## 4. Quickstart: rodar amostra

O fluxo recomendado para rodar uma amostra é via `make run`, que faz:
1. staging da amostra (symlink por padrão, com gzip automático);
2. preparação do banco (FASTA + BLAST + Bowtie2);
3. pipeline completo.

Exemplo (pareado):

```bash
make run ID=demo \
  R1=/caminho/reads_R1.fastq.gz \
  R2=/caminho/reads_R2.fastq.gz \
  DB=ptv
```

Para forçar cópia em vez de symlink, use `COPY=1` no `make sample-add` (ou configure staging manualmente).

Também é possível usar diretamente o script principal:

```bash
./scripts/20_run_pipeline.sh --install --sample demo --kmer 31
```

Se existir `config/picornavirus.env` (ou `config.env` legado), ele é usado como base (ex.: `ASSEMBLER=spades`, `VELVET_K`, `THREADS`).

## 5. DBs suportados

O DB manager oferece perfis básicos e suporta queries personalizadas:

```bash
make db-list
make db DB=ptv
```

Para alterar a query NCBI:

```bash
make db DB=ptv DB_QUERY='"Teschovirus"[Organism] AND "complete genome"[Title]'
```

Se existir `config.env`, ele é usado como base (ex.: `ASSEMBLER=spades`, `VELVET_K`, `THREADS`). Também dá para chamar via `make pipeline`.

---

## 5. Painel de uso (UX)

Para facilitar a execução por usuários não técnicos, existe um painel web simples. Dentro do WSL:

```bash
make ux
```

Depois, no Windows, acesse: `http://localhost:8000`.

O painel permite checar o ambiente, gerar DEMO, importar amostras e rodar o pipeline (Velvet por padrão, SPAdes opcional). Consulte detalhes em `docs/PAINEL_UX.md`.

---

## 6. Banco viral genérico (dashboard + script)

Agora o projeto suporta DB parametrizável (sem hardcode de PTV) via:

```bash
bash scripts/10_build_viral_db.sh --target teschovirus_a
bash scripts/10_build_viral_db.sh --target enterovirus_g --query '"Enterovirus G"[Organism]'
bash scripts/10_build_viral_db.sh --target sapelovirus_a --taxid 12115
```

Catálogo versionado em `config/targets.json`.

Estrutura gerada:

```text
db/<target>/
  ├─ sequences.fasta
  ├─ blastdb/
  └─ metadata.json
```

`metadata.json` registra query/taxid, data, número de sequências e fonte/versão.

## 7. Início rápido no Windows

Use o launcher:

```bat
bundle\start_platform.bat
```

Ele valida `wsl.exe`, valida distro (default: `Ubuntu`), converte paths Windows→WSL e inicia o dashboard em background, abrindo o navegador em `http://localhost:8000`.
