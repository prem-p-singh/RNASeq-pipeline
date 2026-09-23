#!/usr/bin/env python3
"""Deterministic raw-read fixture: real reference downloads, fastp and Salmon.

One unique transcript per gene, with planted up/down regulation and independently
known fragment counts. This is technical qualification, not a real-data benchmark.
"""
import csv
import gzip
from contextlib import ExitStack
import json
from pathlib import Path
import random
import sys

import yaml


def build(dest, repo, layout):
    rng = random.Random(9122026)
    source = dest / "inputs"
    source.mkdir(parents=True)
    (dest / "config").mkdir()
    seqs = {f"G{i:03d}": "".join(rng.choices("ACGT", k=2000)) for i in range(120)}
    (source / "transcripts.fa").write_text("".join(f">{g}_tx\n{s}\n" for g, s in seqs.items()))
    (source / "genes.gtf").write_text("".join(
        f'chr1\tfixture\texon\t{i*3000+1}\t{i*3000+2000}\t.\t+\t.\tgene_id "{g}"; transcript_id "{g}_tx";\n'
        for i, g in enumerate(seqs)))
    expected = {}
    sheet = []
    for sample in range(6):
        # Exercise numeric-looking IDs through the entire single-end workflow.
        sid = f"{sample+1:03d}" if layout == "single" else f"S{sample+1}"
        condition = "control" if sample < 3 else "treated"
        r1, r2 = source / f"{sid}_R1.fq.gz", source / f"{sid}_R2.fq.gz"
        expected[sid] = {}
        with gzip.open(r1, "wt") as a, gzip.open(r2, "wt") as b:
            for i, (g, seq) in enumerate(seqs.items()):
                effect = (4 if i < 8 else 0.25 if i < 16 else 1) if sample >= 3 else 1
                n = round(600 * effect * rng.uniform(0.9, 1.1))
                expected[sid][g] = n
                for j in range(n):
                    length = rng.randint(250, 400)
                    start = rng.randrange(len(seq) - length)
                    frag = seq[start:start+length]
                    left = frag[:100]
                    right = frag[-100:].translate(str.maketrans("ACGT", "TGCA"))[::-1]
                    a.write(f"@{sid}_{g}_{j}/1\n{left}\n+\n{'I'*100}\n")
                    b.write(f"@{sid}_{g}_{j}/2\n{right}\n+\n{'I'*100}\n")
        sheet.append([sid, str(r1), str(r2) if layout == "paired" else "", condition])
    with open(dest / "config/samples.tsv", "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample_id", "fastq_url", "fastq_url_r2", "treatment"])
        w.writerows(sheet)
    cfg = {
        "project": {"name": f"raw_{layout}", "output_dir": "results/"},
        "organism": {"scientific_name": "Fixtura synthetica", "tax_id": 1,
                     "kegg_code": None, "orgdb_package": "org.Fixture.eg.db"},
        "reference": {"accession": "FIXTURE_RAW_1", "cache_dir": None,
                      "transcriptome_fasta_url": (source / "transcripts.fa").as_uri(),
                      "gtf_url": (source / "genes.gtf").as_uri(),
                      "annotation_tsv": {"path": None}},
        "samples": {"seq_type": f"rnaseq_{layout}", "sheet": "config/samples.tsv"},
        "model": {"fixed_effects": "~ treatment", "random_effects": None,
                  "primary_factor": "treatment"},
        "contrasts": [{"id": "treatment", "type": "pairwise", "factor": "treatment", "reverse": True}],
        "downstream": {"run_go": False, "run_kegg": False, "run_wgcna": False},
        "orgdb": {"strategy": "skip"},
        "hpc": {"delete_fastq_after_quant": True},
        "thresholds": {"sample_qc": {"min_reads_on_genes_rnaseq": 1}},
    }
    (dest / "config/config.yaml").write_text(yaml.safe_dump(cfg))
    (dest / "config/thresholds.yaml").write_text((repo / "config/thresholds.yaml").read_text())
    (dest / "expected.json").write_text(json.dumps(expected))


def canonicalize(dest, repo):
    """Use two sequencing runs per library, each reusing lane number 1."""
    sys.path.insert(0, str(repo / "scripts"))
    import metadata
    cfg = yaml.safe_load((dest / "config/config.yaml").read_text())
    tables, _ = metadata.from_legacy(metadata.read_tsv(dest / "config/samples.tsv"), cfg["samples"]["seq_type"])
    reads = []
    for record in tables["reads"]:
        source = Path(record["uri"])
        with ExitStack() as stack:
            outputs = [source.with_name(f"run{run}_{source.name}") for run in (1, 2)]
            handles = [stack.enter_context(gzip.open(p, "wt")) for p in outputs]
            with gzip.open(source, "rt") as handle:
                i = 0
                while header := handle.readline():
                    handles[i % 2].write(header + ''.join(handle.readline() for _ in range(3)))
                    i += 1
        for run, path in enumerate(outputs, 1):
            reads.append(dict(record, read_unit_id=f"run{run}_{record['read_unit_id']}",
                              uri=str(path), run=str(run), lane="1"))
    tables["reads"] = reads
    for lib in tables["libraries"]:
        lib.update(umi="no", strandedness="unknown")
    metadata.write_tables(tables, dest / "metadata")
    cfg["samples"].update(metadata_dir="metadata", sheet="metadata/samples.tsv")
    (dest / "config/config.yaml").write_text(yaml.safe_dump(cfg))


def verify(dest):
    expected = json.loads((dest / "expected.json").read_text())
    for sample, counts in expected.items():
        path = dest / "results/quant" / sample
        metrics = json.loads((path / "metrics.json").read_text())
        assert metrics["mapping_rate"] > 0.98, metrics
        assert metrics["num_reads_processed"] == sum(counts.values())
        with open(path / "quant.sf") as f:
            got = {r["Name"].removesuffix("_tx"): float(r["NumReads"]) for r in csv.DictReader(f, delimiter="\t")}
        assert set(got) == set(counts)
        assert max(abs(got[g]-counts[g]) / counts[g] for g in counts) < 0.03
        assert (dest / "inputs" / f"{sample}_R1.fq.gz").is_file(), "source read deleted"
        assert not list(path.glob("*.trim.fastq.gz")), "owned intermediates retained"
    for rel in ("counts.tsv", "de_manifest.tsv", "enrichment_status.tsv", "wgcna_manifest.tsv",
                "qc_report/qc_charts.html", "qc_report/multiqc/multiqc_report.html"):
        assert (dest / "results" / rel).stat().st_size > 0, rel
    with open(dest / "results/counts.tsv") as f:
        assert next(csv.reader(f, delimiter="\t"))[1:] == list(expected), "sample IDs changed in counts"
    with open(dest / "results/de_manifest.tsv") as f:
        manifest = list(csv.DictReader(f, delimiter="\t"))
    with open(dest / "results" / manifest[0]["analysis_table"]) as f:
        de = {r["gene_id"]: r for r in csv.DictReader(f, delimiter="\t")}
    for i in range(16):
        row = de[f"G{i:03d}"]
        assert float(row["logFC"]) * (1 if i < 8 else -1) > 1, row
        assert float(row["adj.P.Val"]) < 0.05, row
    print("Raw-read fixture: reference, quantification, DE signs, reports and cleanup passed")


if __name__ == "__main__":
    mode, dest = sys.argv[1], Path(sys.argv[2]).resolve()
    if mode == "build":
        build(dest, Path(sys.argv[3]).resolve(), sys.argv[4])
    elif mode == "canonicalize":
        canonicalize(dest, Path(sys.argv[3]).resolve())
    else:
        verify(dest)
