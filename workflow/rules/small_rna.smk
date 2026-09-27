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

# Small-RNA plant route: ShortStack clusters loci jointly, so it runs once for the cohort.
# Reads are fastp-trimmed already, so shortstack_command never passes --autotrim.
if PLANT_SRNA_ROUTE:
    ruleorder: small_rna_plant_metrics > small_rna_counts

    rule shortstack_run:
        input:
            genome = REF / "genome.fa",
            known = [REF / "mature_species.fa"] if config["reference"].get("mature_mirna_url") else [],
            trims = expand(str(QUANT / "{sample}/{sample}.trim.fastq.gz"), sample=SAMPLES),
            helper = REPO_DIR / "scripts/small_rna.py",
        output:
            counts = OUT / "small_rna_counts.tsv",
            loci = OUT / "small_rna_loci.tsv",
        params:
            known = str(REF / "mature_species.fa") if config["reference"].get("mature_mirna_url") else "",
            runtime_identity = RUNTIME_ID,
        threads: 8
        resources:
            mem_mb = 32000,
            runtime = 360,
        log:
            "logs/shortstack.log",
        shell:
            """
            set -euo pipefail
            tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
            mkdir "$tmp/reads"
            # {sample}.fastq.gz symlinks so ShortStack's per-file count columns are the sample IDs.
            for f in {input.trims:q}; do ln -s "$(readlink -f "$f")" "$tmp/reads/$(basename "$f" .trim.fastq.gz).fastq.gz"; done
            python3 -c 'import subprocess, sys; sys.path.insert(0, sys.argv[1]); from small_rna import shortstack_command; subprocess.run(shortstack_command(sys.argv[2], sys.argv[6:], sys.argv[3], sys.argv[4], sys.argv[5]), check=True)' \
                {REPO_DIR:q}/scripts {input.genome:q} {params.known:q} "$tmp/run" {threads} "$tmp"/reads/*.fastq.gz > {log:q} 2>&1
            python3 {input.helper:q} shortstack_matrix "$tmp/run/Counts.txt" "$tmp/run/Results.txt" {output.counts:q} {output.loci:q} >> {log:q} 2>&1
            """

    rule small_rna_plant_metrics:
        input:
            counts = OUT / "small_rna_counts.tsv",
            preprocessing = QUANT / "{sample}/preprocessing.json",
            helper = REPO_DIR / "scripts/small_rna.py",
        output:
            metrics = QUANT / "{sample}/small_rna/metrics.json",
        params:
            runtime_identity = RUNTIME_ID,
        resources:
            mem_mb = 2000,
            runtime = 10,
        shell:
            "python3 {input.helper:q} plant_metrics {input.counts:q} {input.preprocessing:q} {wildcards.sample:q} {output.metrics:q}"
