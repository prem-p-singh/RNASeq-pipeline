#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from recommend import recommend

cfg = {"samples": {"seq_type": "rnaseq_paired"},
       "model": {"fixed_effects": "~ condition", "random_effects": None}}
a, b = recommend(cfg, 6), recommend(cfg, 300)
assert a["route"] == b["route"] == "salmon_bulk"
assert a["backend"] == b["backend"] == "limma_voom"
assert a["handbook_ids"] and a["policy_sha256"] and a["rule_ids"]
cfg["model"]["random_effects"] = "(1|subject)"
assert recommend(cfg)["backend"] == "dream"
for analysis in ({"objectives": ["splicing"]}, {"umi": True}, {"backend": "deseq2"}):
    cfg["analysis"] = analysis
    assert recommend(cfg)["status"] == "blocked"
cfg.pop("analysis")
for seq in ("small_rna", "single_cell", "spatial", "long_read", "unknown"):
    cfg["samples"]["seq_type"] = seq
    assert recommend(cfg)["status"] == "blocked"
print("check_recommend.py: scope, dependence, traceability and sample-count invariance passed")
