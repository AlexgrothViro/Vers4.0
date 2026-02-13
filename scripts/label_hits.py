#!/usr/bin/env python3
import argparse
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("adj_path", nargs="?", default="run_T1/work/ptv_hits.adjust.tsv")
    args = ap.parse_args()

    path = Path(args.adj_path)
    if not path.exists():
        raise SystemExit(f"[ERRO] Arquivo não encontrado: {path}")

    with path.open("r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            # Espera colunas geradas pelo adj_identity.py (pelo menos 13)
            if len(cols) < 13:
                continue
            try:
                qseqid = cols[0]
                alen   = int(float(cols[3]))
                evalue = cols[4]
                bits   = cols[5]
                consider = int(cols[10])
                match    = int(cols[11])
                adj      = float(cols[12])
            except Exception:
                # ignora linhas mal-formadas
                continue

            label = "REVIEW"
            if alen >= 60 and adj >= 90.0:
                label = "PASS"
            elif 40 <= alen < 60 and adj >= 95.0:
                label = "SHORT_PASS"

            print("\t".join([
                qseqid, label, f"{adj:.2f}", str(alen), evalue, bits, str(consider), str(match)
            ]))

if __name__ == "__main__":
    main()
