#!/usr/bin/env python3
import argparse, gzip, random
from pathlib import Path

def read_first_fasta_seq(path: Path) -> str:
    seq = []
    with path.open("rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq:
                    break
                continue
            seq.append(line)
    s = "".join(seq).upper()
    if not s:
        raise SystemExit(f"[ERRO] FASTA sem sequência útil: {path}")
    return s

def revcomp(s: str) -> str:
    comp = str.maketrans("ACGTN", "TGCAN")
    return s.translate(comp)[::-1]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sample", default="DEMO")
    ap.add_argument("--pairs", type=int, default=2000)
    ap.add_argument("--len", dest="read_len", type=int, default=150)
    ap.add_argument("--insert", type=int, default=300)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    ref = Path(args.ref)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    seq = read_first_fasta_seq(ref)
    if len(seq) < args.insert:
        raise SystemExit("[ERRO] Referência curta demais para o insert atual.")
    if args.read_len > args.insert:
        raise SystemExit("[ERRO] read_len maior que insert; ajuste --len ou --insert.")

    random.seed(args.seed)
    r1_path = outdir / f"{args.sample}_R1.fastq.gz"
    r2_path = outdir / f"{args.sample}_R2.fastq.gz"
    qual = "I" * args.read_len

    max_start = len(seq) - args.insert
    with gzip.open(r1_path, "wt", encoding="utf-8") as r1, gzip.open(r2_path, "wt", encoding="utf-8") as r2:
        for i in range(1, args.pairs + 1):
            start = random.randint(0, max_start)
            frag = seq[start:start + args.insert]
            read1 = frag[:args.read_len]
            read2 = revcomp(frag[-args.read_len:])
            rid = f"{args.sample}_{i}"
            r1.write(f"@{rid}/1\n{read1}\n+\n{qual}\n")
            r2.write(f"@{rid}/2\n{read2}\n+\n{qual}\n")

    print(f"[OK] DEMO gerada: {r1_path} e {r2_path}")

if __name__ == "__main__":
    main()
