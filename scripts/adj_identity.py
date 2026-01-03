import sys
inp = sys.argv[1]; out = sys.argv[2]
def is_base(c): return c in "ACGTacgt"
with open(inp) as fi, open(out,"w") as fo:
    for line in fi:
        line=line.rstrip("\n")
        if not line: continue
        c=line.split("\t")
        if len(c) < 14: continue
        qseqid,sseqid=c[0],c[1]; pid,alen=c[2],c[3]
        evalue,bits=c[6],c[7]; qstart,qend=c[8],c[9]; sstart,send=c[10],c[11]
        qseq,sseq=c[12],c[13]
        consider=match=0
        for qc,sc in zip(qseq,sseq):
            if qc in "Nn": continue
            if not is_base(sc): continue
            consider+=1
            if qc.upper()==sc.upper(): match+=1
        adj = (100.0*match/consider) if consider>0 else 0.0
        fo.write("\t".join([qseqid,sseqid,pid,alen,evalue,bits,qstart,qend,sstart,send,
                            str(consider),str(match),f"{adj:.2f}",qseq,sseq])+"\n")
