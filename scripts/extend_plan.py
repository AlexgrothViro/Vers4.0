#!/usr/bin/env python3
import sys, csv

REPORT = "run_T1/work/ptv_report.tsv"
REF_FA = "data/ptv_db.fa"
OUT_TSV = "run_T1/work/extend_plan.tsv"

# parâmetros de janela de extensão (pode ajustar depois)
FLANK_MIN = 60     # tamanho mínimo desejado por lado
FLANK_MAX = 150    # tamanho máximo desejado por lado

# Carrega referências em memória
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
ref_len = {k: len(v) for k,v in ref.items()}

# Le o report e calcula plano
rows = []
with open(REPORT, newline="") as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for d in r:
        label = d.get("label","")
        if label not in ("PASS","SHORT_PASS"):
            # só planeja extensão para aprovados (ajuste se quiser incluir REVIEW)
            continue

        qid = d["qseqid"]
        sid = d["sseqid"]
        L   = int(float(d["length"]))
        qstart = int(d["qstart"])
        qend   = int(d["qend"])
        sstart = int(d["sstart"])
        send   = int(d["send"])

        if sid not in ref_len:
            # referência não encontrada (inconsistência)
            strand = "?"
            Lref = 0
            left_gap = right_gap = 0
            rows.append([qid, sid, strand, L, Lref, sstart, send, left_gap, right_gap, 0,0, "REF_NOT_FOUND"])
            continue

        Lref = ref_len[sid]
        strand = "+" if send >= sstart else "-"  # BLAST coords

        # normaliza intervalos no sentido 5'->3' da referência
        sL = min(sstart, send)
        sR = max(sstart, send)

        # lacunas até as extremidades da referência
        left_gap  = sL - 1            # bases disponíveis à esquerda (antes de sL)
        right_gap = Lref - sR         # bases disponíveis à direita (depois de sR)

        # tamanhos desejados de flanco
        left_w  = min(max(FLANK_MIN, 0), min(FLANK_MAX, left_gap))
        right_w = min(max(FLANK_MIN, 0), min(FLANK_MAX, right_gap))

        # status (se tem material para estender)
        status = []
        status.append("HAS_LEFT" if left_w  > 0 else "NO_LEFT")
        status.append("HAS_RIGHT" if right_w > 0 else "NO_RIGHT")

        rows.append([qid, sid, strand, L, Lref, sstart, send,
                     left_gap, right_gap, left_w, right_w, ",".join(status)])

# escreve plano
with open(OUT_TSV, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["qseqid","sseqid","strand","aln_len","ref_len","sstart","send",
                "left_gap","right_gap","left_w","right_w","status"])
    for r in rows:
        w.writerow(r)

print("[OK] extend_plan.tsv")
