#!/usr/bin/env python3
import sys, csv

PLAN   = "run_T1/work/extend_plan.tsv"
REF_FA = "data/ptv_db.fa"
OUTFA  = "run_T1/work/extend_regions.fasta"

def read_fasta(path):
    seqs = {}
    cur_id = None
    cur = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(cur)
                cur_id = line[1:].split()[0]
                cur = []
            else:
                cur.append(line.replace(" ", "").replace("\t", ""))
        if cur_id is not None:
            seqs[cur_id] = "".join(cur)
    return seqs

ref = read_fasta(REF_FA)

def emit_record(out, name, seq, width=70):
    out.write(f">{name}\n")
    for i in range(0, len(seq), width):
        out.write(seq[i:i+width] + "\n")

with open(PLAN, newline="") as fh, open(OUTFA, "w") as out:
    r = csv.DictReader(fh, delimiter="\t")
    for d in r:
        sid = d["sseqid"]
        if sid not in ref:
            continue
        genome = ref[sid]
        Lref = len(genome)

        qid    = d["qseqid"]
        strand = d["strand"]
        sstart = int(d["sstart"])
        send   = int(d["send"])
        left_w  = int(d["left_w"])
        right_w = int(d["right_w"])

        sL = min(sstart, send)
        sR = max(sstart, send)

        # fatias de extensão, respeitando limites
        if left_w > 0:
            start = max(1, sL - left_w)  # 1-based
            end   = sL - 1
            if end >= start:
                subseq = genome[start-1:end]  # python 0-based
                emit_record(out, f"{qid}|{sid}|left|{start}-{end}|strand={strand}", subseq)

        if right_w > 0:
            start = sR + 1
            end   = min(Lref, sR + right_w)
            if end >= start:
                subseq = genome[start-1:end]
                emit_record(out, f"{qid}|{sid}|right|{start}-{end}|strand={strand}", subseq)

print("[OK] extend_regions.fasta")
