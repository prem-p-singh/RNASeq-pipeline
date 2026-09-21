# Sample-count strategy (auto-selected)

The pipeline picks one of three strategies **automatically** based on how many rows are in `config/samples.tsv`. You don't choose — the submit script counts and decides, then writes the choice to `gates/decisions.log`.

---

## The three tiers

| Tier | # samples | How it runs | Why |
|---|---|---|---|
| **Small** | **≤ 20** | One SLURM job loops through samples serially (~30 min/sample, ~10 hr total worst case) | Job-array overhead isn't worth it; easier to debug; single log file |
| **Medium** | **21 – 200** | SLURM job array, **up to 20 samples running at once** | Sweet spot — parallel speedup without flooding the queue |
| **Large** | **201+** | Chunked array of **50 at a time**, throttled (`--jobs 20`, `--max-jobs-per-second 1`) | Being a good queue citizen; checkpoints after each chunk so interruptions don't cost much |

All three run the **same Snakemake rules** — the only difference is the SLURM profile they submit through.

---

## What actually changes between tiers

| Setting | Small | Medium | Large |
|---|---|---|---|
| SLURM submission | single job | job array | chunked job array |
| Max concurrent samples | 1 | 20 | 20 (but 50 queued) |
| `--jobs` (Snakemake) | 1 | 20 | 20 |
| `--max-jobs-per-second` | — | — | 1 |
| Peak disk usage | 1 FASTQ (~0.5 GB) | 20 FASTQ (~10 GB) | 20 FASTQ (~10 GB) |
| Wall time (rough) | hours | minutes–hours | hours–day |
| Checkpoint granularity | per-sample | per-sample | per-chunk |

The 20 GB cap is the ceiling — medium/large tiers are already close to it at peak (10 GB of FASTQ in flight + 1 GB index + ~3 GB outputs). If a medium project has unusually big files, the submit script drops you to a stricter `--jobs` value automatically.

---

## How the auto-selection works

When you run `./submit.sh`:

1. Counts rows in `config/samples.tsv` → `N`
2. Estimates per-sample FASTQ size from `samples.seq_type` in `config.yaml` (TAGseq ≈ 0.5 GB, RNA-Seq ≈ 2 GB)
3. Computes **max safe concurrency** = `floor((20 GB - 4 GB reserved) / avg_fastq_size)`
4. Picks tier:
   - `N ≤ 20` → small
   - `N ≤ 200` → medium (concurrency capped by step 3)
   - `N > 200` → large (concurrency capped by step 3, chunked)
5. Writes decision to `gates/decisions.log`:
   ```
   [2026-04-18 14:22] STRATEGY: medium (N=87, avg_fastq=0.5GB, max_concurrent=20)
   ```
6. Launches `snakemake --profile profiles/<tier>/ ...`

---

## Manual override

If the auto-choice is wrong (e.g., you want to stress-test with 5 samples on the medium profile), pass `--profile` directly:

```bash
snakemake --profile profiles/medium/ --configfile config/config.yaml
```

This skips `submit.sh` entirely.

---

## Edge cases the tiers handle

| Situation | What happens |
|---|---|
| Sample FASTQ larger than remaining budget | `submit.sh` reduces `--jobs` until safe; logs warning |
| Queue rejects submission | Snakemake retries with exponential backoff (built-in) |
| Mid-run crash | `snakemake` picks up from last completed `quant.sf` — nothing re-runs |
| One sample fails validation (mapping < 60%) | That sample is flagged in `gates/decisions.log`; downstream R scripts exclude it |
| All samples finish before aggregate step | Aggregate stage (tximport + DE + GO + WGCNA) runs on login-node-friendly resources, not array |
