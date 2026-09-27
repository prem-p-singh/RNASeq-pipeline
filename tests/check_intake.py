#!/usr/bin/env python3
"""Exercise the shipped workbook through setup, without live reference lookups."""
import importlib.util
import json
import os
import subprocess
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
    assert record["schema_version"] == 5
    inputs = setup_inputs(record)
    assert inputs["tax_id"] == 9606 and inputs["seq_type"] == "rnaseq_paired"
    assert inputs["expected_libtype"] == "ISR" and inputs["model_random_effects"] is None
    assert inputs["metadata_file"] == str(base / "samples.tsv")
    assert inputs["orgdb_strategy"] == "skip" and not inputs["run_enrichment"]
    assert "platform" not in record["context_only"] and "library_kit" in record["context_only"]
    legacy=openpyxl.load_workbook(book)
    legacy['README']['B2']=3
    legacy['Bulk RNA-seq'].delete_rows(30,10)
    legacy.save(base/'schema3.xlsx')
    assert setup_inputs(read_intake(base/'schema3.xlsx'))['analysis']=={'quantifier':'auto','backend':'auto'}

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
    rejected("repeated_measures", "yes", "subject column", True)
    rejected("subject_column", "donor", "only used with repeated measures", True)
    rejected("model_formula", "~ genotype", "must include", True)
    rejected("model_formula", "treatment + x", "use ~ followed", True)
    s[cells["analysis_goal"]] = "QC only"; w.save(base / "qc.xlsx")
    assert setup_inputs(read_intake(base / "qc.xlsx"))["objectives"] == ["qc"]
    s[cells["run_wgcna"]] = "yes"; w.save(base / "qc.xlsx"); s[cells["run_wgcna"]] = "no"
    try:
        setup_inputs(read_intake(base / "qc.xlsx"))
    except ValueError as exc:
        assert "WGCNA needs expression" in str(exc)
    else:
        raise AssertionError("QC-only goal accepted WGCNA")
    s[cells["analysis_goal"]] = "differential expression"
    rejected("umi", "yes", "UMI processing", True)
    rejected("umi", "don't know", "unknown status", True)
    rejected("input_stage", "counts or processed objects", "Only raw FASTQ", True)
    rejected("umi", None, "required")
    rejected('analysis_backend', 'magic', 'choose one')
    rejected('preprocess_minimum_length', 0, 'whole number')
    rejected('preprocess_minimum_length', 10001, '1..10000', True)
    rejected('preprocess_adapter_r1', 'INVALID', 'A/C/G/T', True)
    rejected('preprocess_adapter_r2', 'ACGT', 'adapter_r1', True)
    rejected('reference_genome', 'custom.fa', 'Incomplete custom reference', True)

    # Explicit methods/references/settings survive the Excel -> setup contract.
    advanced = dict(analysis_quantifier='star', analysis_backend='edger_ql',
                    reference_genome='genome.fa', reference_gtf='genes.gtf',
                    preprocess_poly_g='off', preprocess_adapter_r1='AGATCGGA',
                    preprocess_minimum_length=25)
    saved_answers = {key:s[cells[key]].value for key in advanced}
    for key,value in advanced.items(): s[cells[key]]=value
    w.save(base/'advanced.xlsx')
    advanced_inputs = setup_inputs(read_intake(base/'advanced.xlsx'))
    assert advanced_inputs['analysis'] == {'quantifier':'star','backend':'edger_ql'}
    assert advanced_inputs['reference_overrides']['genome_fasta_url'] == (base/'genome.fa').as_uri()
    assert advanced_inputs['reference_overrides']['transcriptome_fasta_url'] is None
    assert advanced_inputs['preprocessing']['minimum_length'] == 25
    assert advanced_inputs['preprocessing']['poly_g'] == 'off'
    rejected('strandedness', "don't know", 'STAR requires known', True)
    for key,value in saved_answers.items(): s[cells[key]]=value

    # Schema 5: goal, subject model, formula, kit profile, platform poly-G and screening.
    schema5 = dict(analysis_goal='expression only', repeated_measures='yes', subject_column='donor',
                   subject_model='fixed block', model_formula='~ donor + treatment + age',
                   library_protocol='TruSeq Stranded mRNA', screening='yes',
                   screening_extra='Botrytis=refs/botrytis.fa', screening_fragments=5000)
    saved_answers = {key:s[cells[key]].value for key in schema5}
    for key,value in schema5.items(): s[cells[key]]=value
    w.save(base/'schema5.xlsx')
    five = setup_inputs(read_intake(base/'schema5.xlsx'))
    assert five['objectives'] == ['gene_expression'], five['objectives']
    assert five['model_fixed_effects'] == '~ donor + treatment + age' and five['model_random_effects'] is None
    assert five['expected_libtype'] == 'ISR' and five['preprocessing']['poly_g'] == 'on'
    assert five['preprocessing']['adapter_r1'] == 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
    assert five['preprocessing']['adapter_r2'] == 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
    refs = five['screening']['references']
    assert [r['name'] for r in refs] == ['Host', 'PhiX', 'Botrytis'] and refs[0] == dict(name='Host', role='expected', fasta='transcriptome')
    assert refs[2]['fasta'] == str(base/'refs/botrytis.fa') and five['screening']['fragments'] == 5000
    s[cells['subject_model']] = 'random'
    s[cells['model_formula']] = None
    w.save(base/'schema5r.xlsx')
    random5 = setup_inputs(read_intake(base/'schema5r.xlsx'))
    assert random5['model_fixed_effects'] == '~ treatment' and random5['model_random_effects'] == '(1|donor)'
    rejected('strandedness', 'forward', 'conflicts with TruSeq Stranded mRNA', True)
    rejected('screening_extra', 'Botrytis', 'name=FASTA', True)
    rejected('run_enrichment', 'yes', 'differential expression goal', True)
    for key,value in saved_answers.items(): s[cells[key]]=value

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
    assert cfg['reference']['decoys'] == 'genome'
    assert cfg["samples"]["expected_libtype"] == "ISR"
    assert not any(cfg["downstream"].values())
    assert cfg["model"]["fixed_effects"] == "~ treatment" and cfg["model"]["random_effects"] is None
    assert "001\t" in (project / "config/samples.tsv").read_text(), "leading zeros lost"
    assert "NA\t" in (project / "config/samples.tsv").read_text(), "literal NA identifier lost"
    assert cfg["contrasts"][0]["factor"] == "treatment"
    assert cfg["analysis"]["backend"] == "deseq2", "handbook default: 3-12 units per group"
    assert cfg["analysis"]["objectives"] == ["gene_expression", "differential_expression"]
    saved = next((project / "intake").glob("*/source_map.json"))
    assert json.loads(saved.read_text())["sources"]["read_config"].startswith("Bulk RNA-seq!C")
    advanced_project = base/'advanced project'
    with patch.object(sys, 'argv', ['setup.py','--intake',str(base/'advanced.xlsx'),
                                   '--project-dir',str(advanced_project),'--plan-only']), \
         patch.object(setup, 'resolve_organism', return_value=org), \
         patch.object(setup, 'resolve_reference_urls', side_effect=AssertionError('Custom reference replaced by lookup')):
        setup.main()
    advanced_cfg=yaml.safe_load((advanced_project/'config/config.yaml').read_text())
    assert advanced_cfg['analysis']['quantifier']=='star' and advanced_cfg['analysis']['backend']=='edger_ql'
    assert advanced_cfg['reference']['genome_fasta_url']==(base/'genome.fa').as_uri()
    assert advanced_cfg['reference']['transcriptome_fasta_url'] is None
    assert advanced_cfg['preprocessing']['minimum_length']==25
    before = (project / "config/config.yaml").read_bytes()
    with patch.object(sys, "argv", argv):
        try:
            setup.main()
        except SystemExit as exc:
            assert "will not overwrite" in str(exc)
        else:
            raise AssertionError("existing configuration overwritten")
    assert before == (project / "config/config.yaml").read_bytes()

    # One launcher command resumes exactly the imported workbook and keeps user config.
    env = dict(os.environ, PATH=str(Path(sys.executable).parent)+os.pathsep+os.environ['PATH'])
    command = ['bash', str(ROOT/'submit.sh'), '--plan-only', '--executor', 'local',
               '--intake', str(book), '-d', str(project)]
    launched = subprocess.run(command, env=env, text=True, capture_output=True)
    assert launched.returncode == 0, launched.stdout+launched.stderr
    assert 'executor:      local' in launched.stdout
    assert before == (project / 'config/config.yaml').read_bytes()
    rejected_resume = subprocess.run(command+['--intake-sheet','TAGseq'], env=env, text=True, capture_output=True)
    assert rejected_resume.returncode != 0 and 'worksheet differs' in rejected_resume.stderr

    # A changed workbook reports its impact; only --accept-revision applies it.
    s[cells["preprocess_minimum_length"]] = 20
    w.save(base / "revised.xlsx")
    revised_cmd = command[:command.index(str(book))] + [str(base / "revised.xlsx")] + command[command.index(str(book)) + 1:]
    report = subprocess.run(revised_cmd, env=env, text=True, capture_output=True)
    assert report.returncode != 0 and "preprocessing: preprocessing, screening" in report.stdout, report.stdout + report.stderr
    assert "--accept-revision" in report.stderr and before == (project / "config/config.yaml").read_bytes()
    applied = subprocess.run(revised_cmd + ["--accept-revision"], env=env, text=True, capture_output=True)
    assert applied.returncode == 0, applied.stdout + applied.stderr
    revised_cfg = yaml.safe_load((project / "config/config.yaml").read_text())
    assert revised_cfg["preprocessing"]["minimum_length"] == 20
    old_digest = json.loads(saved.read_text())["workbook_sha256"]
    new_digest = (project / "intake/current").read_text().strip()
    assert new_digest != old_digest and (project / "intake" / old_digest / "config.yaml").read_bytes() == before
    impact = json.loads((project / "intake" / new_digest / "revision_impact.json").read_text())
    assert impact["changes"] == [{"field": "preprocessing", "stages": "preprocessing, screening, quantification, counts, statistics, report"}]
    stale = subprocess.run(command, env=env, text=True, capture_output=True)
    assert stale.returncode != 0 and "Workbook revision impact" in stale.stdout, "old workbook resumed over revised config"
    s[cells["preprocess_minimum_length"]] = 15

    # Optional Contrasts tab replaces generated contrasts.
    for col, value in enumerate(["trt", "treatment", None, "treated", "control"], 1):
        w["Contrasts"].cell(5, col, value)
    w["Contrasts"].cell(6, 1, "coef"); w["Contrasts"].cell(6, 6, "treatmenttreated=1")
    w.save(base / "contrasts.xlsx")
    specs = setup_inputs(read_intake(base / "contrasts.xlsx"))["contrasts"]
    assert specs == [dict(id="trt", type="pairwise", factor="treatment", by=None, numerator="treated", denominator="control"),
                     dict(id="coef", type="linear", weights={"treatmenttreated": 1.0})], specs
    w["Contrasts"].cell(6, 2, "treatment"); w.save(base / "contrasts.xlsx")
    try:
        setup_inputs(read_intake(base / "contrasts.xlsx"))
    except ValueError as exc:
        assert "use weights alone" in str(exc)
    else:
        raise AssertionError("mixed linear/pairwise contrast accepted")
    for row in (5, 6):
        for col in range(1, 7):
            w["Contrasts"].cell(row, col).value = None

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
