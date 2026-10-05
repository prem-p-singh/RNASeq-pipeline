# Host and pathogen in one library: combined reference, read assignment and split.
# Included only when dual_organism.enabled; single-organism runs never parse it.
DUAL_REF = REF / "dual"
DUAL_SPLIT = [str(QUANT / ("{sample}/dual/%s%s.fastq.gz" % (role, mate))) for role in ("host", "pathogen")
              for mate in (("_R1", "_R2") if len(_trim_patterns) == 2 else ("",))]

rule dual_fetch:
    output:
        fa = DUAL_REF / "{role}.genome.fa",
    params:
        url = lambda wc: DUAL[wc.role]["genome_fasta_url"],
    wildcard_constraints:
        role = "host|pathogen",
    shell:
        "bash {REPO_DIR:q}/workflow/scripts/fetch_reference_file.sh --url {params.url:q} --out {output.fa:q} --kind fasta"

rule dual_fetch_gtf:
    output:
        gtf = DUAL_REF / "{role}.annotation.gtf",
    params:
        url = lambda wc: DUAL[wc.role]["gtf_url"],
    wildcard_constraints:
        role = "host|pathogen",
    shell:
        "bash {REPO_DIR:q}/workflow/scripts/fetch_reference_file.sh --url {params.url:q} --out {output.gtf:q} --kind gtf"

# Contig names get the organism label and "|" as prefix.
rule dual_combine:
    input:
        host = DUAL_REF / "host.genome.fa",
        pathogen = DUAL_REF / "pathogen.genome.fa",
        host_gtf = DUAL_REF / "host.annotation.gtf",
        pathogen_gtf = DUAL_REF / "pathogen.annotation.gtf",
        helper = REPO_DIR / "scripts/dual_organism.py",
    output:
        fa = DUAL_REF / "combined.fa",
        gtf = DUAL_REF / "combined.gtf",
    params:
        host_label = DUAL["host"]["label"],
        pathogen_label = DUAL["pathogen"]["label"],
        runtime_identity = RUNTIME_ID,
    resources:
        mem_mb = 4000,
        runtime = 60,
    log:
        "logs/dual_combine.log",
    shell:
        "python3 {input.helper:q} combine --host {input.host:q} --pathogen {input.pathogen:q} "
        "--host-gtf {input.host_gtf:q} --pathogen-gtf {input.pathogen_gtf:q} "
        "--host-label {params.host_label:q} --pathogen-label {params.pathogen_label:q} "
        "--out-fa {output.fa:q} --out-gtf {output.gtf:q} > {log:q} 2>&1"

rule dual_index:
    input:
        genome = DUAL_REF / "combined.fa",
        gtf = DUAL_REF / "combined.gtf",
        helper = REPO_DIR / "scripts/star_counts.py",
    output:
        index = directory(DUAL_REF / "star_idx"),
        lock = DUAL_REF / "star_idx/reference.lock.json",
    params:
        overhang = config["reference"]["star_overhang"],
        runtime_identity = RUNTIME_ID,
    threads: 4
    resources:
        mem_mb = 64000,
        runtime = 240,
    log:
        "logs/dual_index.log",
    shell:
        "python3 {input.helper:q} index --genome {input.genome:q} --gtf {input.gtf:q} --out {output.index:q} --threads {threads} --overhang {params.overhang} > {log:q} 2>&1"

rule dual_assign:
    input:
        lock = DUAL_REF / "star_idx/reference.lock.json",
        trims = lambda wc: [p.format(sample=wc.sample) for p in _trim_patterns],
        preprocessing = QUANT / "{sample}/preprocessing.json",
        helper = REPO_DIR / "scripts/dual_organism.py",
    output:
        record = QUANT / "{sample}/dual/assignment.json",
        split = DUAL_SPLIT,
    params:
        index = str(DUAL_REF / "star_idx"),
        outdir = lambda wc: str(QUANT / wc.sample / "dual"),
        host_label = DUAL["host"]["label"],
        pathogen_label = DUAL["pathogen"]["label"],
        limit = DUAL["multimap_limit"],
        runtime_identity = RUNTIME_ID,
    threads: 4
    resources:
        mem_mb = 40000,
        runtime = 180,
    log:
        "logs/dual/{sample}.log",
    shell:
        "python3 {input.helper:q} assign --index {params.index:q} --reads {input.trims:q} "
        "--outdir {params.outdir:q} --sample {wildcards.sample:q} "
        "--host-label {params.host_label:q} --pathogen-label {params.pathogen_label:q} "
        "--multimap-limit {params.limit} --threads {threads} > {log:q} 2>&1"

rule dual_table:
    input:
        records = expand(str(QUANT / "{sample}/dual/assignment.json"), sample=SAMPLES),
        helper = REPO_DIR / "scripts/dual_organism.py",
    output:
        table = OUT / "assignment.tsv",
    params:
        warn = DUAL["ambiguous_warn_fraction"],
    shell:
        "python3 {input.helper:q} table --warn {params.warn} --out {output.table:q} {input.records:q}"
