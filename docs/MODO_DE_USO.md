# Modo de uso do pipeline (visão rápida)

Este documento resume como rodar o pipeline de PTV em uma amostra, incluindo a parte avançada
(identidade ajustada, planejamento de flancos, FASTA de flancos e simulação de reads).

## 1. Pré-requisitos

- Linux ou WSL2 com Ubuntu
- Programas instalados:
  - velveth, velvetg
  - blastn, makeblastdb
  - bowtie2
  - python3
  - curl + EDirect (para baixar sequências do NCBI)

O repositório deve estar clonado e os diretórios básicos criados:

    git clone https://github.com/AlexgrothViro/quali-2025.git
    cd quali-2025
    make setup_dirs

## 2. Preparar bancos de referência (uma vez por máquina)

### 2.1. Banco de Porcine Teschovirus (PTV) para BLAST

    cd quali-2025
    ./scripts/10_build_ptv_db.sh

    # definir variável BLAST_DB (se ainda não estiver definida)
    export BLAST_DB="\$PWD/db/ptv_teschovirus"
### 2.2. Genoma de Sus scrofa + índice Bowtie2

    cd quali-2025
    ./scripts/11_download_sus_scrofa.sh
    bowtie2-build ref/host/sus_scrofa.fa ref/host/sus_scrofa_bt2

## 3. Colocar os FASTQs da amostra

Copiar os arquivos FASTQ.gz pareados para `data/raw` com o padrão:

- data/raw/<SAMPLE>_R1.fastq.gz  
- data/raw/<SAMPLE>_R2.fastq.gz  

Exemplo real:

- data/raw/81554_S150_R1.fastq.gz  
- data/raw/81554_S150_R2.fastq.gz  

## 4. Rodar o pipeline avançado para uma amostra

O script abaixo executa:

1. make test (checagem de ambiente, filtro de hospedeiro, montagem com Velvet, BLAST básico)
2. BLAST de confirmação com qseq/sseq
3. identidade ajustada
4. relatório resumido de hits
5. plano de flancos na referência PTV
6. FASTA das regiões de flanco
7. simulação de reads de flanco + rótulos dos hits
Exemplo de uso:

    cd quali-2025

    # garantir variável BLAST_DB
    export BLAST_DB="\$PWD/db/ptv_teschovirus"

    ./scripts/run_ptv_advanced.sh 81554_S150 31

As principais saídas ficam em:

- run_T1/work/ptv_hits.confirm.tsv  
- run_T1/work/ptv_hits.adjust.tsv  
- run_T1/work/ptv_report.tsv  
- run_T1/work/extend_plan.tsv  
- run_T1/work/extend_regions.fasta  
- run_T1/work/sim_R1.fastq.gz  
- run_T1/work/sim_R2.fastq.gz  
- run_T1/work/ptv_hits.labels.tsv  

Esses arquivos são a base para análise, geração de tabelas, gráficos e etapas de filogenia.
