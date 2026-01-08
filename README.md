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

## 4. Comando único (quick start)

Para rodar toda a verificação e o pipeline básico em um único comando, use:

```bash
./scripts/20_run_pipeline.sh --install --sample 81554_S150 --kmer 31
```

Esse script faz:
1. checagem/instalação de dependências;
2. preparação de diretórios e bancos (FASTA + BLAST + Bowtie2);
3. montagem dos contigs (Velvet por padrão);
4. BLAST dos contigs contra o banco de PTV.

Se existir `config.env`, ele é usado como base (ex.: `ASSEMBLER=spades`, `VELVET_K`, `THREADS`). Também dá para chamar via `make pipeline`.
