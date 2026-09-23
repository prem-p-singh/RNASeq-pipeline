# Current architecture

The implemented release is a Linux bulk RNA-seq workflow. Broader RNA assay families remain planned; [README.md](README.md) and [release qualification](docs/RELEASE.md) define the actual supported boundary.

1. **Environment:** install the exact Linux package lock into a dedicated prefix; verify package metadata, versions, imports and runtime paths before analysis. Keep custom OrgDb packages outside that prefix.
2. **Preflight:** resolve configuration defaults, validate sample identity and requested capability, test model/contrast estimability, and write issues and recommendation reasons. Reject unsupported requested methods.
3. **Resources:** measure local inputs and available filesystem space; estimate unknown inputs, cache construction, retained outputs and concurrent working sets. Insufficient or unknown capacity only warns. Storage never blocks launch or reduces concurrency. Actual write failures remain errors.
4. **Reference:** download through temporary files, validate format, check FASTA-to-GTF transcript/gene compatibility, build Salmon under a filesystem lock, and atomically publish the index. Record source checksums, construction settings and the actual builder version.
5. **Quantification:** run fastp and Salmon per library. Keep seekable trimmed reads during quantification, validate quantitative outputs, record QC and input hashes, then clean only owned intermediates according to policy.
6. **Aggregation:** apply the sample-disposition policy, verify the retained design, import counts with the route's length policy, and save counts, effective lengths, tximport data and provenance.
7. **Inference:** use limma-voom for fixed effects or dream for declared random effects. Record contrast identities and result paths. Optional enrichment uses signed statistics; WGCNA uses explicit settings and memory-aware blocking.
8. **Reporting/restart:** produce comparative QC and MultiQC. Validate result manifests before reusing completion markers; let Snakemake rerun changed inputs/configuration or missing outputs.

The source checkout is separate from each project's output tree. The environment and any shared reference cache must be visible on every SLURM worker. Node-local temporary directories are suitable for isolated release tests, but cannot substitute for shared project paths in a distributed run.

The cache key includes URLs and index parameters, not a remote content lookup. Source checksums detect local corruption; immutable source URLs remain necessary. Reference download/index publication is safe per artifact, but publication of the whole bundle is not a single transaction. See the release limitations before enabling concurrent fresh builds against a shared cache.

The three-table metadata validators are a foundation for future library/lane execution. The current DAG consumes the legacy single sample sheet; validation support must not be advertised as executable multi-lane or multi-assay support.
