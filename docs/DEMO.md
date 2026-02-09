# DEMO (execução reprodutível)

Rodar:
1) make deps
2) make demo
3) BLAST_DB=blastdb/ptv BOWTIE2_INDEX=bowtie2/ptv ASSEMBLER=velvet make pipeline SAMPLE=DEMO

Saídas esperadas:
- data/assemblies/DEMO_velvet_k31/contigs.fa
- results/blast/DEMO_k31_vs_db.tsv
- results/reports/DEMO_summary.md
