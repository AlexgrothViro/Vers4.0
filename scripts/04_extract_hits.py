#!/usr/bin/env python3
import argparse, os, sys, statistics as st
from collections import defaultdict

def read_fasta(path):
    seqs = {}
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs
def pick_input_paths(sample, kmer, tsv_arg=None, contigs_arg=None):
    tsv_candidates = []
    if tsv_arg:
        tsv_candidates.append(tsv_arg)
    tsv_candidates += [
        f"results/blast/{sample}_k{kmer}_vs_db.tsv",
        f"results/blast/{sample}_vs_db.tsv",  # fallback legado
    ]
    tsv = next((p for p in tsv_candidates if os.path.exists(p)), None)
    if not tsv:
        raise FileNotFoundError("Nenhum TSV encontrado. Tente gerar com make test-blast e/ou informe --tsv.")

    contig_candidates = []
    if contigs_arg:
        contig_candidates.append(contigs_arg)
    contig_candidates += [
        f"data/assemblies/{sample}_velvet_k{kmer}/contigs.fa",
        f"data/assemblies/{sample}_velvet_k{kmer}/contigs.fa",  # redundante, intencional
    ]
    contigs = next((p for p in contig_candidates if os.path.exists(p)), None)
    if not contigs:
        raise FileNotFoundError("contigs.fa não encontrado para esse sample/kmer. Rode make test-velvet e/ou informe --contigs.")

    return tsv, contigs

def fasta_len_stats(lengths):
    if not lengths:
        return "nenhuma sequência escrita"
    s = sorted(lengths)
    def pct(p): return s[int((p/100)*(len(s)-1))]
    return (
        f"min={min(s)} p50={pct(50)} p90={pct(90)} max={max(s)} mean={round(st.mean(s),2)} | "
        f">=50={sum(x>=50 for x in s)} >=80={sum(x>=80 for x in s)} >=120={sum(x>=120 for x in s)} >=200={sum(x>=200 for x in s)}"
    )

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--kmer", required=True, type=int)
    ap.add_argument("--ref", default="data/ref/ptv_db.fa")
    ap.add_argument("--min-pident", type=float, default=90.0)
    ap.add_argument("--min-aln-len", type=int, default=50)
    ap.add_argument("--tsv", default=None, help="opcional: caminho explícito do TSV BLAST")
    ap.add_argument("--contigs", default=None, help="opcional: caminho explícito do contigs.fa")
    args = ap.parse_args()

    sample = args.sample
    k = args.kmer

    tsv, contigs_fa = pick_input_paths(sample, k, args.tsv, args.contigs)

    outdir = f"results/phylogeny/{sample}_k{k}"
    os.makedirs(outdir, exist_ok=True)

    print("== extract_hits ==")
    print("sample:", sample, "kmer:", k)
    print("TSV:", tsv)
    print("contigs.fa:", contigs_fa)
    print("ref.fa:", args.ref)
    print(f"filtro: pident>={args.min_pident} aln_len>={args.min_aln_len}")

    contigs = read_fasta(contigs_fa)
    refs_all = read_fasta(args.ref)

    best = {}  # qseqid -> (bits, alen, pident, sseqid, sstart, send, evalue)
    refs_needed = set()

    with open(tsv) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 12:
                continue
            qseqid, sseqid = cols[0], cols[1]
            pident = float(cols[2])
            alen = int(cols[3])
            sstart = cols[8]
            send   = cols[9]
            evalue = cols[10]
            bits   = float(cols[11])

            if pident < args.min_pident or alen < args.min_aln_len:
                continue

            cur = best.get(qseqid)
            if (cur is None) or (bits > cur[0]):
                best[qseqid] = (bits, alen, pident, sseqid, sstart, send, evalue)
                refs_needed.add(sseqid)
    hits_fa = os.path.join(outdir, "hits.fa")
    summary_tsv = os.path.join(outdir, "hits_summary.tsv")
    refs_fa_out = os.path.join(outdir, "refs.fa")

    missing = []
    written_lengths = []

    with open(hits_fa, "w") as out_fa, open(summary_tsv, "w") as out_tsv:
        out_tsv.write("contig\tcontig_len\tbest_ref\tpident\taln_len\tbitscore\tevalue\tsstart\tsend\n")
        for qseqid, (bits, alen, pident, sseqid, sstart, send, evalue) in sorted(best.items(), key=lambda x: x[1][0], reverse=True):
            seq = contigs.get(qseqid)
            if not seq:
                missing.append(qseqid)
                continue
            out_fa.write(f">{qseqid}\n{seq}\n")
            out_tsv.write(f"{qseqid}\t{len(seq)}\t{sseqid}\t{pident}\t{alen}\t{bits}\t{evalue}\t{sstart}\t{send}\n")
            written_lengths.append(len(seq))

    miss_refs = []
    with open(refs_fa_out, "w") as out:
        for acc in sorted(refs_needed):
            seq = refs_all.get(acc)
            if not seq:
                miss_refs.append(acc)
                continue
            out.write(f">{acc}\n{seq}\n")

    print("OUTDIR:", outdir)
    print("Contigs_com_hit (no TSV):", len(best))
    print("Contigs_escritos (encontrados no contigs.fa):", len(written_lengths))
    print("Refs:", len(refs_needed))
    print("hits.fa:", hits_fa)
    print("hits_summary.tsv:", summary_tsv)
    print("refs.fa:", refs_fa_out)
    print("Tamanho_seqs_hits:", fasta_len_stats(written_lengths))

    if missing:
        print(f"ATENÇÃO: {len(missing)} contigs do TSV não foram encontrados no contigs.fa (provável TSV/assembly desencontrado).", file=sys.stderr)
        print("Exemplos (até 10):", ", ".join(missing[:10]), file=sys.stderr)

    if miss_refs:
        print("ATENÇÃO: refs não encontradas no ref.fa (até 10):", ", ".join(miss_refs[:10]), file=sys.stderr)

if __name__ == "__main__":
    main()
