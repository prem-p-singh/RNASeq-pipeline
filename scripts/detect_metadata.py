#!/usr/bin/env python3
"""
detect_metadata.py <metadata_path>

Inspect a metadata CSV or XLSX and emit a JSON summary:
  - n_samples (rows)
  - columns
  - sample_id_candidates (by header-name heuristic)
  - factor_candidates  (columns with 2-10 unique values, good DE factors)

Used by new_project.sh to auto-fill setup_inputs.yaml with sensible guesses.
"""
from __future__ import annotations
import json
import re
import sys
from pathlib import Path

import pandas as pd


SAMPLE_ID_PATTERNS = [
    r"^sample[_\s]?id$", r"^sample$", r"^sampleid$",
    r"^specimen[_\s]?id$", r"^specimen$",
    r"^id$", r"^name$",
]


def classify(df: pd.DataFrame) -> dict:
    cols = [str(c) for c in df.columns]

    # Sample ID candidates — header regex match, case insensitive
    sample_id_candidates = []
    for c in cols:
        for pat in SAMPLE_ID_PATTERNS:
            if re.match(pat, c.strip(), flags=re.IGNORECASE):
                sample_id_candidates.append(c)
                break

    # Factor candidates — 2 to 10 unique values, not a sample-id column
    factor_candidates = []
    for c in cols:
        if c in sample_id_candidates:
            continue
        try:
            n_unique = int(df[c].nunique(dropna=True))
        except Exception:
            continue
        if 2 <= n_unique <= 10:
            levels = [str(v) for v in df[c].dropna().unique().tolist()]
            factor_candidates.append({
                "name": c,
                "n_unique": n_unique,
                "levels": levels,
            })

    return {
        "n_samples": len(df),
        "columns": cols,
        "sample_id_candidates": sample_id_candidates,
        "factor_candidates": factor_candidates,
    }


def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: detect_metadata.py <metadata_path>")
    path = Path(sys.argv[1]).expanduser()
    if not path.exists():
        sys.exit(f"File not found: {path}")

    if path.suffix.lower() in {".xlsx", ".xls"}:
        df = pd.read_excel(path)
    elif path.suffix.lower() in {".csv"}:
        df = pd.read_csv(path)
    elif path.suffix.lower() in {".tsv", ".txt"}:
        df = pd.read_csv(path, sep="\t")
    else:
        sys.exit(f"Unknown file type: {path.suffix}")

    print(json.dumps(classify(df), indent=2, default=str))


if __name__ == "__main__":
    main()
