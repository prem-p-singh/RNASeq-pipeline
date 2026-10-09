# Genomic alignment route. Barcode/UMI libraries never enter these rules.
rule star_index:
    input:
        genome = REF / "genome.fa",
        gtf = REF / "annotation.gtf",
        helper = REPO_DIR / "scripts/star_counts.py",
    output:
        index = directory(REF / "star_idx"),
        lock = REF / "star_idx/reference.lock.json",
    params:
        overhang = config["reference"]["star_overhang"],
        runtime_identity = RUNTIME_ID,
    threads: 4
    resources:
        mem_mb = 64000,
        runtime = 180,
    log:
        "logs/star_index.log",
    shell:
        "python3 {input.helper:q} index --genome {input.genome:q} --gtf {input.gtf:q} --out {output.index:q} --threads {threads} --overhang {params.overhang} > {log:q} 2>&1"


def star_strand(wc):
    declared = (CANONICAL_INPUTS[wc.sample]["expected_libtype"] if CANONICAL_INPUTS is not None
                else config["samples"].get("expected_libtype"))
    return {"U": 0, "IU": 0, "SF": 1, "ISF": 1, "SR": 2, "ISR": 2}[declared]


rule star_quant:
    input:
        index = REF / "star_idx/reference.lock.json",
        gtf = REF / "annotation.gtf",
        trims = lambda wc: [p.format(sample=wc.sample) for p in _trim_patterns],
        preprocessing = QUANT / "{sample}/preprocessing.json",
        helper = REPO_DIR / "scripts/star_counts.py",
    output:
        folder = directory(QUANT / "{sample}/star"),
        counts = QUANT / "{sample}/star/gene_counts.tsv",
        lengths = QUANT / "{sample}/star/gene_lengths.tsv",
        bam = QUANT / "{sample}/star/Aligned.sortedByCoord.out.bam",
        bai = QUANT / "{sample}/star/Aligned.sortedByCoord.out.bam.bai",
        metrics = QUANT / "{sample}/star/metrics.json",
        alignment_summary = QUANT / "{sample}/star/Log.final.out",
        counting_summary = QUANT / "{sample}/star/featureCounts.tsv.summary",
    params:
        index = str(REF / "star_idx"),
        r1 = lambda wc: _trim_patterns[0].format(sample=wc.sample),
        r2 = lambda wc: _trim_patterns[1].format(sample=wc.sample) if len(_trim_patterns) == 2 else "",
        strand = star_strand,
        runtime_identity = RUNTIME_ID,
    threads: 4
    resources:
        mem_mb = 40000,
        runtime = 180,
    log:
        "logs/star/{sample}.log",
    shell:
        """
        python3 {input.helper:q} quant --r1 {params.r1:q} --r2={params.r2:q} \
            --index {params.index:q} --gtf {input.gtf:q} --out {output.folder:q} \
            --sample {wildcards.sample:q} --strand {params.strand} --threads {threads} \
            --replace-owned > {log:q} 2>&1
        """
