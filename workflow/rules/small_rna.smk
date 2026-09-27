# Small-RNA animal route: Bowtie 1 to species mature miRNAs, unique hits counted.
# Tools come from the srna module (submit.sh adds it for small_rna_* routes).
rule mature_mirna_reference:
    input:
        helper = REPO_DIR / "scripts/small_rna.py",
    output:
        fa = REF / "mature_species.fa",
    params:
        url = config["reference"]["mature_mirna_url"],
        prefix = config["organism"]["mirbase_prefix"],
        runtime_identity = RUNTIME_ID,
    resources:
        mem_mb = 2000,
        runtime = 30,
    log:
        "logs/mature_mirna_reference.log",
    shell:
        """
        set -euo pipefail
        bash {REPO_DIR:q}/workflow/scripts/fetch_reference_file.sh --url {params.url:q} --out {output.fa:q}.all --kind fasta > {log:q} 2>&1
        python3 {input.helper:q} reference {output.fa:q}.all {params.prefix:q} {output.fa:q} >> {log:q} 2>&1
        bowtie-build {output.fa:q} {output.fa:q} >> {log:q} 2>&1
        """

# -v 0: exact matches; -a --best --strata: report every equally best hit so
# multi-mapping reads are seen and excluded; --norc: mature miRNA reads are sense.
rule small_rna_counts:
    input:
        ref = REF / "mature_species.fa",
        trims = lambda wc: [p.format(sample=wc.sample) for p in _trim_patterns],
        preprocessing = QUANT / "{sample}/preprocessing.json",
        helper = REPO_DIR / "scripts/small_rna.py",
    output:
        counts = QUANT / "{sample}/small_rna/counts.tsv",
        metrics = QUANT / "{sample}/small_rna/metrics.json",
    params:
        runtime_identity = RUNTIME_ID,
    threads: 4
    resources:
        mem_mb = 4000,
        runtime = 60,
    log:
        "logs/small_rna/{sample}.log",
    shell:
        """
        set -euo pipefail
        tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
        bowtie -p {threads} -v 0 -a --best --strata --norc -S -x {input.ref:q} \
            {input.trims:q} "$tmp/a.sam" 2> {log:q}
        python3 {input.helper:q} count "$tmp/a.sam" {output.counts:q} {output.metrics:q} {wildcards.sample:q} >> {log:q} 2>&1
        """
