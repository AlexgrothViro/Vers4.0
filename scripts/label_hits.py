#!/usr/bin/env python3
import sys

ADJ_PATH = sys.argv[1] if len(sys.argv) > 1 else "run_T1/work/ptv_hits.adjust.tsv"

with open(ADJ_PATH, "r") as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        cols = line.split("\t")
        # Espera as mesmas colunas geradas pelo adj_identity.py
        if len(cols) < 15:
            continue

        qseqid = cols[0]
        alen   = int(float(cols[3]))
        evalue = cols[4]
        bits   = cols[5]
        consider = int(cols[10])
        match    = int(cols[11])
        adj      = float(cols[12])

        label = "REVIEW"
        if alen >= 60 and adj >= 90.0:
            label = "PASS"
        elif 40 <= alen < 60 and adj >= 95.0:
            label = "SHORT_PASS"

        print("\t".join([
            qseqid, label, f"{adj:.2f}", str(alen), evalue, bits, str(consider), str(match)
        ]))
