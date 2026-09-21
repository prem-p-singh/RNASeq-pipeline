# =============================================================================
# Shared helpers — utilities reused by retrieve / per_sample / aggregate rules.
# =============================================================================

import json
from datetime import datetime


def log_decision(gate: str, message: str):
    """Append an auto-decision to gates/decisions.log with timestamp."""
    line = f"[{datetime.now().isoformat(timespec='seconds')}] {gate}: {message}\n"
    with open(GATES / "decisions.log", "a") as f:
        f.write(line)


def write_metrics(stage: str, payload: dict):
    """Persist a stage's metrics as JSON for downstream gates to read."""
    path = METRICS / f"{stage}.json"
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)
    return path


def read_metrics(stage: str) -> dict:
    path = METRICS / f"{stage}.json"
    if not path.exists():
        return {}
    with open(path) as f:
        return json.load(f)


def get_fastq_url(wildcards):
    """Look up the (R1 / single-end) FASTQ URL for a given sample from the sheet."""
    return samples.loc[wildcards.sample, "fastq_url"]


def get_fastq_url_r2(wildcards):
    """Look up the R2 FASTQ URL for paired-end samples.

    Returns "" when there is no fastq_url_r2 column (single-end / TAGseq) or
    when the cell is blank, so the per-sample script can stay single-end.
    """
    if "fastq_url_r2" not in samples.columns:
        return ""
    val = samples.loc[wildcards.sample, "fastq_url_r2"]
    return "" if pd.isna(val) else val


def sample_metadata_col(sample: str, col: str):
    """Fetch a metadata field for a sample (used by aggregate scripts)."""
    return samples.loc[sample, col]
