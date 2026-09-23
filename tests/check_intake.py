#!/usr/bin/env python3
"""Exercise the shipped workbook through setup, without live reference lookups."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

import openpyxl
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from intake import read_intake, setup_inputs

spec = importlib.util.spec_from_file_location("project_setup", ROOT / "scripts/setup.py")
setup = importlib.util.module_from_spec(spec)
spec.loader.exec_module(setup)

with tempfile.TemporaryDirectory(prefix="intake test ") as tmp:
    base = Path(tmp).resolve()
    book = base / "filled intake.xlsx"
    w = openpyxl.load_workbook(ROOT / "intake_template.xlsx")
    s = w["Bulk RNA-seq"]
    answers = {
        "project_name": "intake_test", "contact_email": "test@example.org",
        "biological_question": "Treatment response", "organism": "human",
        "tissue_type": "blood", "platform": "NextSeq", "primary_factor": "treatment",
        "min_samples_per_group": 3, "repeated_measures": "no", "batch_effect_column": "none",
        "fastq_path": "reads", "metadata_file": "samples.tsv", "sample_id_column": "sample_id",
        "run_enrichment": "no", "orgdb_strategy": "auto", "run_wgcna": "no",
        "library_type": "polyA mRNA-seq", "read_config": "paired-end", "strandedness": "reverse",
        "input_stage": "raw FASTQ", "umi": "no",
    }
    q = yaml.safe_load((ROOT / "config/intake_questions.yaml").read_text())
    # The shipped workbook is a generated interface to the same schema.
    assert w["README"]["B2"].value == q["schema_version"]
    for tab in q["tabs"]:
        expected = q["common_questions"] + tab["extra_questions"]
        actual = list(w[tab["name"]].iter_rows(min_row=5, max_row=4 + len(expected), max_col=5))
        for question, row in zip(expected, actual):
            assert row[4].value == question["id"], (tab["name"], question["id"])
            assert row[1].value == question["text"]
            assert row[3].value == question["help"]
            assert (row[2].value or None) == (question.get("default") or None)
    questions = q["common_questions"] + q["tabs"][0]["extra_questions"]
    labels = {item["text"]: item["id"] for item in questions}
    cells = {}
    for row in s.iter_rows(min_col=1, max_col=3):
        key = labels.get(row[1].value)
        if key:
            cells[key] = row[2].coordinate
            if key in answers:
                row[2].value = answers[key]
    w.save(book)
    record = read_intake(book)
    assert record["schema_version"] == 3
    inputs = setup_inputs(record)
    assert inputs["tax_id"] == 9606 and inputs["seq_type"] == "rnaseq_paired"
    assert inputs["expected_libtype"] == "ISR" and inputs["model_random_effects"] is None
    assert inputs["metadata_file"] == str(base / "samples.tsv")
    assert inputs["orgdb_strategy"] == "skip" and not inputs["run_enrichment"]
    assert "platform" in record["context_only"]

    # Moving the question rows must not change interpretation.
    original = [(s.cell(5, c).value, s.cell(6, c).value) for c in range(1, 6)]
    for c, (a, b) in enumerate(original, 1):
        s.cell(5, c, b); s.cell(6, c, a)
    w.save(base / "reordered.xlsx")
    assert setup_inputs(read_intake(base / "reordered.xlsx")) == inputs
    for c, (a, b) in enumerate(original, 1):
        s.cell(5, c, a); s.cell(6, c, b)

    def rejected(key, value, message, convert=False):
        old = s[cells[key]].value
        s[cells[key]] = value
        w.save(base / "bad.xlsx")
        try:
            r = read_intake(base / "bad.xlsx")
            if convert:
                setup_inputs(r)
        except ValueError as exc:
            assert message in str(exc), str(exc)
        else:
            raise AssertionError(f"accepted invalid {key}")
        s[cells[key]] = old

    rejected("project_name", "=1+1", "formula")
    rejected("read_config", "mystery", "choose one")
    rejected("metadata_file", None, "required")
    rejected("min_samples_per_group", 2.5, "whole number")
    rejected("expected_libtype", "ISF", "conflicts", True)
    rejected("repeated_measures", "yes", "subject/model", True)
    rejected("umi", "yes", "UMI processing", True)
    rejected("umi", "don't know", "unknown status", True)
    rejected("input_stage", "counts or processed objects", "Only raw FASTQ", True)
    rejected("umi", None, "required")

    w["README"]["B2"] = 99
    w.save(base / "future.xlsx")
    try:
        read_intake(base / "future.xlsx")
    except ValueError as exc:
        assert "unsupported intake schema" in str(exc)
    else:
        raise AssertionError("accepted unknown schema version")
    w["README"]["B2"] = 2

    # Stable IDs survive relabeling, including the project-selection question.
    project_row = s[cells["project_name"]].row
    old_label = s.cell(project_row, 2).value
    s.cell(project_row, 2, "Study identifier")
    w.save(base / "relabelled.xlsx")
    assert setup_inputs(read_intake(base / "relabelled.xlsx")) == inputs
    s.cell(project_row, 2, old_label)

    reads = base / "reads"
    reads.mkdir()
    rows = ["sample_id\ttreatment"]
    for i in range(1, 7):
        sid = "NA" if i == 6 else f"00{i}"
        rows.append(f"{sid}\t{'control' if i <= 3 else 'treated'}")
        for mate in (1, 2):
            (reads / f"{sid}_R{mate}.fastq.gz").write_bytes(b"fixture")
    (base / "samples.tsv").write_text("\n".join(rows) + "\n")
    project = base / "project"
    org = dict(genus="Homo", species="sapiens", tax_id=9606, kegg_code="hsa", orgdb_package="org.Hs.eg.db")
    ref = dict(accession="fixture", assembly_name="fixture",
               genome_fasta_url="file:///fixture/genome.fa", transcriptome_fasta_url="file:///fixture/tx.fa",
               gtf_url="file:///fixture/genes.gtf", gene_info_url=None)
    argv = ["setup.py", "--intake", str(book), "--project-dir", str(project)]
    with patch.object(sys, "argv", argv), patch.object(setup, "resolve_organism", return_value=org), \
         patch.object(setup, "resolve_reference_urls", return_value=ref), \
         patch.object(setup, "fetch_annotation_info", return_value=False):
        setup.main()
    cfg = yaml.safe_load((project / "config/config.yaml").read_text())
    assert cfg["samples"]["expected_libtype"] == "ISR"
    assert not any(cfg["downstream"].values())
    assert cfg["model"]["fixed_effects"] == "~ treatment" and cfg["model"]["random_effects"] is None
    assert "001\t" in (project / "config/samples.tsv").read_text(), "leading zeros lost"
    assert "NA\t" in (project / "config/samples.tsv").read_text(), "literal NA identifier lost"
    assert cfg["contrasts"][0]["factor"] == "treatment"
    saved = next((project / "intake").glob("*/source_map.json"))
    assert json.loads(saved.read_text())["sources"]["read_config"].startswith("Bulk RNA-seq!C")
    before = (project / "config/config.yaml").read_bytes()
    with patch.object(sys, "argv", argv):
        try:
            setup.main()
        except SystemExit as exc:
            assert "will not overwrite" in str(exc)
        else:
            raise AssertionError("existing configuration overwritten")
    assert before == (project / "config/config.yaml").read_bytes()

    # The same setup must consume workbook records without scanning filenames.
    s[cells["metadata_source"]] = "workbook tables"
    s[cells["metadata_file"]] = None
    s[cells["fastq_path"]] = None
    for rownum, line in enumerate(rows[1:], 5):
        sid, condition = line.split("\t")
        w["Samples"].cell(rownum, 1, sid)
        w["Samples"].cell(rownum, 2, sid)
        w["Samples"].cell(rownum, 3, condition)
        for col, value in enumerate([f"lib_{sid}", sid, "bulk", "paired", "reverse", "no"], 1):
            w["Libraries"].cell(rownum, col, value)
        for mate in (1, 2):
            readrow = 5 + (rownum - 5) * 2 + mate - 1
            for col, value in enumerate([f"read_{sid}_{mate}", f"lib_{sid}", f"R{mate}",
                                         f"reads/{sid}_R{mate}.fastq.gz", "run1", "1"], 1):
                w["Reads"].cell(readrow, col, value)
    w.save(base / "inline.xlsx")
    inline_project = base / "inline project"
    argv = ["setup.py", "--intake", str(base / "inline.xlsx"), "--project-dir", str(inline_project)]
    with patch.object(sys, "argv", argv), patch.object(setup, "resolve_organism", return_value=org), \
         patch.object(setup, "resolve_reference_urls", return_value=ref), \
         patch.object(setup, "fetch_annotation_info", return_value=False), \
         patch.object(setup, "build_samples_tsv", side_effect=AssertionError("Inline setup scanned files")):
        setup.main()
    import metadata
    cfg = yaml.safe_load((inline_project / "config/config.yaml").read_text())
    assert cfg["samples"]["sheet"] == "metadata/samples.tsv"
    compiled = metadata.execution_inputs(inline_project, cfg)
    assert set(compiled) == {"001", "002", "003", "004", "005", "NA"}
    assert compiled["001"]["units"][0]["r1"]["uri"] == str(reads / "001_R1.fastq.gz")
    w["Samples"]["A5"] = 1
    w.save(base / "numeric-id.xlsx")
    try:
        read_intake(base / "numeric-id.xlsx")
    except ValueError as exc:
        assert "identifiers must be stored as Text" in str(exc)
    else:
        raise AssertionError("Numeric sample identifier accepted")
    w["Samples"]["A5"] = "001"
    w["Samples"]["B6"] = "001"
    w.save(base / "repeat-unit.xlsx")
    try:
        setup_inputs(read_intake(base / "repeat-unit.xlsx"))
    except ValueError as exc:
        assert "Repeated biological units" in str(exc)
    else:
        raise AssertionError("Repeated biological units treated as independent")
    w["Samples"]["B6"] = "002"
    w["Reads"]["B5"] = "missing_library"
    w.save(base / "orphan.xlsx")
    try:
        read_intake(base / "orphan.xlsx")
    except ValueError as exc:
        assert "no matching" in str(exc)
    else:
        raise AssertionError("Orphan read accepted")

print("check_intake.py: all assertions passed")
