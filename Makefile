SHELL := /usr/bin/env bash

SCRIPTS_DIR := scripts

# Parâmetros padrão (podem ser sobrescritos: SAMPLE=... KMER=...)
SAMPLE ?= 81554_S150
KMER   ?= 31

.PHONY: help setup_dirs test-env filter-host test-velvet test-blast test

help:
	@echo "Alvos disponíveis:"
	@echo "  make help                 # mostra esta ajuda"
	@echo "  make setup_dirs           # cria a estrutura básica de pastas (data/, results/, docs/)"
	@echo "  make test-env             # verifica se os programas básicos estão instalados"
	@echo "  make filter-host          # remove leituras alinhadas ao genoma de Sus scrofa"
	@echo "  make test-velvet          # roda montagem de teste com Velvet usando SAMPLE e KMER"
	@echo "  make test-blast           # roda BLAST dos contigs contra o banco definido em BLAST_DB"
	@echo "  make test                 # executa test-env, filter-host, test-velvet e test-blast em sequência"

setup_dirs:
	mkdir -p data/raw data/cleaned data/host_removed data/assemblies
	mkdir -p results/qc results/blast results/phylogeny results/reports
	mkdir -p docs

test-env:
	$(SCRIPTS_DIR)/00_check_env.sh

filter-host:
	$(SCRIPTS_DIR)/03_filter_host.sh $(SAMPLE)

test-velvet:
	$(SCRIPTS_DIR)/01_run_velvet.sh $(SAMPLE) $(KMER)

test-blast:
	$(SCRIPTS_DIR)/02_run_blast.sh $(SAMPLE) $(KMER)

test: test-env filter-host test-velvet test-blast
