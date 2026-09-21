#!/usr/bin/env python3
"""build_gene_matrix.py — gene-level count + TPM matrices across all samples.

Collapses each sample's transcript-level quant.sf to gene level using the GTF
(transcript_id -> gene_id) and attaches gene descriptions from the RefSeq
transcriptome FASTA headers. The summed gene counts equal what tximport
produces with countsFromAbundance="no".

Outputs into <out>/:
  gene_counts_matrix.tsv   gene_id, description, <one integer column per sample>
  gene_tpm_matrix.tsv      gene_id, description, <one TPM column per sample>

Standalone helper, not part of the Snakemake DAG: run it by hand when someone
wants a spreadsheet. 02_aggregate.R is what the pipeline itself uses.
"""
import argparse, gzip, os, re, sys


def clean_desc(d, species=""):
    """Tidy a RefSeq FASTA description.

    `species` is stripped when it is present as a leading binomial. It used to
    be hardcoded to one organism, which silently left the species name in place
    on every other genome.
    """
    d = re.sub(r'^PREDICTED:\s*', '', d)
    if species:
        d = re.sub(r'^' + re.escape(species) + r'\s+', '', d)
    d = re.sub(r',\s*(transcript variant [^,]+,\s*)?(partial )?mRNA\.?\s*$', '', d)
    d = re.sub(r'\s*\([A-Za-z0-9_]+\)\s*$', '', d)
    return d.strip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quant", required=True, help="quant/ dir with per-sample subdirs")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--rna", required=True, help="transcriptome FASTA (.fna.gz) for descriptions")
    ap.add_argument("--out", required=True)
    ap.add_argument("--species", default="",
                    help="scientific name to strip from descriptions, e.g. 'Vitis vinifera'")
    ap.add_argument("--samples", default=None,
                    help="sample sheet pinning which samples to include; without it "
                         "every directory under --quant holding a quant.sf is used, "
                         "including any left from an earlier run")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    # transcript -> gene
    tx2gene = {}
    for line in open(a.gtf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        g = re.search(r'gene_id "([^"]+)"', f[8])
        t = re.search(r'transcript_id "([^"]+)"', f[8])
        if g and t:
            tx2gene[t.group(1)] = g.group(1)

    # transcript -> description
    tx_desc = {}
    with gzip.open(a.rna, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                acc = line[1:].split()[0]
                tx_desc[acc] = clean_desc(
                    line[1:].rstrip("\n")[len(acc):].strip(), a.species)

    on_disk = sorted(s for s in os.listdir(a.quant)
                     if os.path.isfile(os.path.join(a.quant, s, "quant.sf")))
    if a.samples:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from qc_report import read_sample_ids
        wanted = read_sample_ids(a.samples)
        samples = [s for s in wanted if s in on_disk]
        absent = [s for s in wanted if s not in on_disk]
        if absent:
            print(f"WARNING {len(absent)} sample(s) in the sheet have no quant.sf "
                  f"and are omitted: {', '.join(absent)}")
        for s in sorted(set(on_disk) - set(wanted)):
            print(f"ignoring {s}: not in {a.samples}")
    else:
        samples = on_disk
        print(f"no --samples given; including all {len(samples)} directory/ies "
              f"under {a.quant}")
    if not samples:
        raise SystemExit(f"build_gene_matrix: no samples with quant.sf under {a.quant}")
    counts = {}
    tpm = {}
    gene_desc = {}
    for s in samples:
        with open(os.path.join(a.quant, s, "quant.sf")) as q:
            next(q)
            for line in q:
                name, length, efflen, t, n = line.rstrip("\n").split("\t")
                g = tx2gene.get(name, name)
                n = float(n); t = float(t)
                counts.setdefault(g, {}); counts[g][s] = counts[g].get(s, 0.0) + n
                tpm.setdefault(g, {}); tpm[g][s] = tpm[g].get(s, 0.0) + t
                d = tx_desc.get(name, "")
                if g not in gene_desc or n > gene_desc[g][0]:
                    gene_desc[g] = (n, d)

    all_genes = sorted(counts.keys())

    def write(mat, fn, integer):
        with open(os.path.join(a.out, fn), "w") as o:
            o.write("gene_id\tdescription\t" + "\t".join(samples) + "\n")
            for g in all_genes:
                vals = [mat[g].get(s, 0.0) for s in samples]
                vs = ([str(int(round(v))) for v in vals] if integer
                      else ["%.3f" % v for v in vals])
                o.write(f"{g}\t{gene_desc.get(g, (0,''))[1]}\t" + "\t".join(vs) + "\n")

    write(counts, "gene_counts_matrix.tsv", True)
    write(tpm, "gene_tpm_matrix.tsv", False)

    libs = {s: sum(counts[g].get(s, 0.0) for g in all_genes) for s in samples}
    det = {s: sum(1 for g in all_genes if counts[g].get(s, 0.0) > 0) for s in samples}
    print(f"genes: {len(all_genes)}   samples: {len(samples)}")
    print("%-6s %14s %10s" % ("sample", "library_size", "genes>0"))
    for s in samples:
        print("%-6s %14d %10d" % (s, round(libs[s]), det[s]))


if __name__ == "__main__":
    main()
