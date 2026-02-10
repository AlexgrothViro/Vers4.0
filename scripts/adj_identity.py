#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def parse_fasta_lengths(path: Path):
    lengths = {}
    current = None
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                lengths[current] = 0
            elif current is not None:
                lengths[current] += len(line)
    return lengths


def compute_adjusted(blast_path: Path, contigs_path: Path, out_path: Path):
    qlens = parse_fasta_lengths(contigs_path)

    with blast_path.open("r", encoding="utf-8", errors="replace") as blast_in, out_path.open(
        "w", encoding="utf-8", newline=""
    ) as out:
        reader = csv.reader(blast_in, delimiter="\t")
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "evalue",
            "bitscore",
            "qlen",
            "aln_cov",
            "adj_identity",
        ])

        for row in reader:
            if len(row) < 12:
                continue
            qseqid, sseqid = row[0], row[1]
            pident = float(row[2])
            aln_len = float(row[3])
            evalue, bitscore = row[10], row[11]
            qlen = float(qlens.get(qseqid, 0))
            aln_cov = (aln_len / qlen) if qlen > 0 else 0.0
            adj_identity = pident * aln_cov
            writer.writerow(
                [
                    qseqid,
                    sseqid,
                    f"{pident:.3f}",
                    f"{int(aln_len)}",
                    evalue,
                    bitscore,
                    f"{int(qlen)}",
                    f"{aln_cov:.5f}",
                    f"{adj_identity:.3f}",
                ]
            )


def main():
    parser = argparse.ArgumentParser(description="Calcula identidade ajustada por cobertura de alinhamento")
    parser.add_argument("--blast", required=True, help="TSV BLAST outfmt 6 (12 colunas)")
    parser.add_argument("--contigs", required=True, help="FASTA de contigs")
    parser.add_argument("--out", required=True, help="TSV de saída com aln_cov e adj_identity")
    args = parser.parse_args()

    compute_adjusted(Path(args.blast), Path(args.contigs), Path(args.out))


if __name__ == "__main__":
    main()
