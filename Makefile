SHELL := /bin/bash
.ONESHELL:
.SHELLFLAGS := -euo pipefail -c

SCRIPTS_DIR := scripts

# Parâmetros padrão (podem ser sobrescritos: SAMPLE=... KMER=...)
SAMPLE ?= 81554_S150
KMER   ?= 31
ID     ?= amostra_teste
R1     ?=
R2     ?=
SINGLE ?=
DB     ?= ptv
DB_QUERY ?=
KMERS ?=

# Caminhos padrão (podem ser sobrescritos via ambiente ou linha de comando)
REF_DIR       ?= data/ref
REF_FASTA     ?= $(REF_DIR)/ptv_db.fa
BLAST_DB      ?= blastdb/ptv
BOWTIE2_INDEX ?= bowtie2/ptv

.PHONY: help setup_dirs deps test-env filter-host test-velvet test-blast \
	ptv-fasta ptv-fasta-legacy blastdb bowtie2-index pipeline test clean fix-wsl \
	db db-list sample-add run

-include config/picornavirus.env
-include config.env

.PHONY: cfg-all cfg-db cfg-assembly cfg-blast

cfg-all: cfg-db cfg-assembly cfg-blast

cfg-db:
	bash scripts/10_build_custom_db.sh

cfg-assembly:
	bash scripts/run_assembly_router.sh

cfg-blast:
	bash scripts/02_run_blast.sh

help:
	@echo "Alvos disponíveis:"
	@echo "  make deps                  # instala dependências (apt-get) e roda check de ambiente"
	@echo "  make setup_dirs             # cria estrutura básica (data/, results/, docs/)"
	@echo "  make ptv-fasta              # baixa/gera FASTA de PTV em $(REF_FASTA)"
	@echo "  make ptv-fasta-legacy       # cria symlink data/ptv_db.fa -> $(REF_FASTA)"
	@echo "  make blastdb                # gera banco BLAST em $(BLAST_DB) (usa $(REF_FASTA))"
	@echo "  make bowtie2-index          # gera índice Bowtie2 em $(BOWTIE2_INDEX) (usa $(REF_FASTA))"
	@echo "  make db                     # prepara FASTA + BLAST DB + Bowtie2 (via db_manager)"
	@echo "  make db-list                # lista perfis básicos de DB"
	@echo "  make sample-add             # faz staging de amostra em data/raw"
	@echo "  make run                    # sample-add + db + pipeline"
	@echo "  make pipeline               # roda verificação + pipeline completo (scripts/20_run_pipeline.sh)"
	@echo "  make test-env               # verifica dependências básicas"
	@echo "  make smoke-test             # roda smoke test (prep + 90_smoke_test.sh)"
	@echo "  make test                   # roda smoke test (prep + 90_smoke_test.sh)"
	@echo "  make filter-host/test-velvet/test-blast # alvos individuais legados"
	@echo "  make clean                  # remove artefatos gerados (blastdb, bowtie2, run_T1, logs/tmp/results)"
	@echo
	@echo "Variáveis úteis:"
	@echo "  REF_FASTA=$(REF_FASTA)"
	@echo "  BLAST_DB=$(BLAST_DB)"
	@echo "  BOWTIE2_INDEX=$(BOWTIE2_INDEX)"
	@echo "  SAMPLE=$(SAMPLE) KMER=$(KMER)"
	@echo "  ID=$(ID) R1=$(R1) R2=$(R2) SINGLE=$(SINGLE)"
	@echo "  DB=$(DB) DB_QUERY=$(DB_QUERY)"
	@echo "  THREADS=$(THREADS) KMERS=$(KMERS)"

setup_dirs:
	mkdir -p data/raw data/cleaned data/host_removed data/assemblies
	mkdir -p $(REF_DIR) results/qc results/blast results/phylogeny results/reports docs

deps:
	$(SCRIPTS_DIR)/00_check_env.sh --install

fix-wsl:
	bash scripts/00_fix_wsl.sh

test-env:
	BLAST_DB="$(BLAST_DB)" BOWTIE2_INDEX="$(BOWTIE2_INDEX)" \
	$(SCRIPTS_DIR)/00_check_env.sh

ptv-fasta: setup_dirs
	$(SCRIPTS_DIR)/10_fetch_ptv_fasta.sh "$(REF_FASTA)"

ptv-fasta-legacy: ptv-fasta
	mkdir -p data
	ln -sf "$(abspath $(REF_FASTA))" data/ptv_db.fa

blastdb: ptv-fasta
	mkdir -p $(dir $(BLAST_DB))
	if [[ -s "$(BLAST_DB).nhr" && -s "$(BLAST_DB).nin" && -s "$(BLAST_DB).nsq" && "$(BLAST_DB).nhr" -nt "$(REF_FASTA)" ]]; then
		echo "[INFO] BLAST DB já existe e está atualizado: $(BLAST_DB)"
	else
		echo "[INFO] Gerando BLAST DB em $(BLAST_DB) a partir de $(REF_FASTA)"
		makeblastdb -in "$(REF_FASTA)" -dbtype nucl -out "$(BLAST_DB)"
	fi

bowtie2-index: ptv-fasta
	mkdir -p $(dir $(BOWTIE2_INDEX))
	if [[ -s "$(BOWTIE2_INDEX).1.bt2" && -s "$(BOWTIE2_INDEX).2.bt2" && -s "$(BOWTIE2_INDEX).3.bt2" && -s "$(BOWTIE2_INDEX).4.bt2" \
	      && -s "$(BOWTIE2_INDEX).rev.1.bt2" && -s "$(BOWTIE2_INDEX).rev.2.bt2" \
	      && "$(BOWTIE2_INDEX).1.bt2" -nt "$(REF_FASTA)" ]]; then
		echo "[INFO] Índice Bowtie2 já existe e está atualizado: $(BOWTIE2_INDEX)"
	else
		echo "[INFO] Gerando índice Bowtie2 em $(BOWTIE2_INDEX) a partir de $(REF_FASTA)"
		bowtie2-build "$(REF_FASTA)" "$(BOWTIE2_INDEX)"
	fi

pipeline:
	$(SCRIPTS_DIR)/20_run_pipeline.sh

db:
	DB="$(DB)" DB_QUERY="$(DB_QUERY)" REF_FASTA="$(REF_FASTA)" \
	BLAST_DB="$(BLAST_DB)" BOWTIE2_INDEX="$(BOWTIE2_INDEX)" \
	$(SCRIPTS_DIR)/13_db_manager.sh setup

db-list:
	$(SCRIPTS_DIR)/13_db_manager.sh list

sample-add:
	ARGS=(--id "$(ID)")
	if [[ -n "$(R1)" ]]; then ARGS+=(--r1 "$(R1)"); fi
	if [[ -n "$(R2)" ]]; then ARGS+=(--r2 "$(R2)"); fi
	if [[ -n "$(SINGLE)" ]]; then ARGS+=(--single "$(SINGLE)"); fi
	if [[ -n "$(COPY)" ]]; then ARGS+=(--copy); fi
	$(SCRIPTS_DIR)/12_stage_sample.sh "$${ARGS[@]}"

run: sample-add db
	RUN_R1=""
	RUN_R2=""
	RUN_SINGLE=""
	if [[ -n "$(R1)" ]]; then RUN_R1="data/raw/$(ID)_R1.fastq.gz"; fi
	if [[ -n "$(R2)" ]]; then RUN_R2="data/raw/$(ID)_R2.fastq.gz"; fi
	if [[ -n "$(SINGLE)" ]]; then RUN_SINGLE="data/raw/$(ID).fastq.gz"; fi
	SAMPLE_ID="$(ID)" SAMPLE_R1="$$RUN_R1" SAMPLE_R2="$$RUN_R2" SAMPLE_SINGLE="$$RUN_SINGLE" \
	DB="$(DB)" DB_QUERY="$(DB_QUERY)" THREADS="$(THREADS)" KMERS="$(KMERS)" \
	$(SCRIPTS_DIR)/20_run_pipeline.sh

filter-host:
	$(SCRIPTS_DIR)/03_filter_host.sh $(SAMPLE)

test-velvet:
	$(SCRIPTS_DIR)/01_run_velvet.sh $(SAMPLE) $(KMER)

test-blast:
	$(SCRIPTS_DIR)/02_run_blast.sh $(SAMPLE) $(KMER)

smoke-test: test-env ptv-fasta-legacy blastdb bowtie2-index
	BLAST_DB="$(BLAST_DB)" BOWTIE2_INDEX="$(BOWTIE2_INDEX)" \
	$(SCRIPTS_DIR)/90_smoke_test.sh

test: smoke-test

clean:
	rm -rf run_T1 blastdb bowtie2 results logs tmp
# ---- PTV-enriched (Bowtie2 -> Velvet -> BLAST -> hits -> MAFFT addfragments) ----
KMER_PTV ?= 51

.PHONY: ptv-enriched-run
ptv-enriched-run: bowtie2-index blastdb
	STOP_AFTER_ALIGN=$(STOP_AFTER_ALIGN) ./scripts/05_ptv_enriched_run.sh "$(SAMPLE)" "$(KMER_PTV)"

# ---- Bundle WSL (micromamba) ----
BUNDLE_TAG ?= dev

.PHONY: bundle-wsl test-bundle-wsl test-all

bundle-wsl:
	bash scripts/98_build_bundle_wsl.sh "$(BUNDLE_TAG)"

test-bundle-wsl:
	bash scripts/99_test_bundle_wsl.sh "$(BUNDLE_TAG)"

test-all: test test-bundle-wsl
