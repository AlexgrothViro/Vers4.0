# DEMO (execução reprodutível)

Rodar do zero:
1) `make deps`
2) `make demo`
3) `BLAST_DB=blastdb/ptv BOWTIE2_INDEX=bowtie2/ptv ASSEMBLER=velvet make pipeline SAMPLE=DEMO`
4) `make test-demo SAMPLE=DEMO`

Saídas esperadas:
- `data/assemblies/DEMO_velvet_k31/contigs.fa`
- `results/blast/DEMO_k31_vs_db.tsv` (e link `results/blast/DEMO_vs_db.tsv`)
- `results/reports/DEMO_summary.md` (com colunas `aln_cov` e `adj_identity`)
- `results/runs/<timestamp>_DEMO/run.json` + cópias de log/report/blast quando executado via dashboard UX
