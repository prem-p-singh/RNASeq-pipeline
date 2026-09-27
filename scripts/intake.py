"""Read the existing questionnaire workbook into setup inputs (M1 first slice).

No Excel formulas or macros are executed. The legacy template is matched by
question label, never row number; future templates may carry field IDs in E.
"""
from pathlib import Path
import hashlib
import math
import re

import yaml
import metadata

ROOT = Path(__file__).resolve().parent.parent
# NCBI nuccore FASTA for the Illumina PhiX control; identity is hashed at download.
PHIX_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta&retmode=text"


def read_intake(path, sheet_name=None):
    import openpyxl

    path = Path(path).expanduser().resolve()
    definition = yaml.safe_load((ROOT / "config/intake_questions.yaml").read_text())
    tabs = {tab["name"]: tab for tab in definition["tabs"]}
    workbook = openpyxl.load_workbook(path, data_only=False, keep_links=False)
    try:
        version = workbook["README"]["B2"].value if "README" in workbook.sheetnames else None
        version = 1 if version in (None, "") else version
        if type(version) is not int or version not in (1, 2, 3, 4, 5):
            raise ValueError("README!B2: unsupported intake schema version; expected 1 to 5")
        candidates = []
        for name in tabs:
            if name in workbook.sheetnames:
                # Defaults exist on every tab. A project name selects a tab.
                if any((row[4].value == "project_name" or row[1].value == "Project name")
                       and row[2].value not in (None, "")
                       for row in workbook[name].iter_rows(min_col=1, max_col=5)):
                    candidates.append(name)
        if sheet_name is None:
            if len(candidates) != 1:
                raise ValueError("Fill Project name on exactly one assay tab, or select --intake-sheet explicitly.")
            sheet_name = candidates[0]
        if sheet_name not in tabs or sheet_name not in workbook.sheetnames:
            raise ValueError(f"Unknown intake tab: {sheet_name!r}")
        questions = definition["common_questions"] + tabs[sheet_name]["extra_questions"]
        by_label = {q["text"]: q for q in questions}
        by_id = {q["id"]: q for q in questions}
        answers, sources, errors = {}, {}, []
        for row in workbook[sheet_name].iter_rows(min_col=1, max_col=5):
            label, cell, field = row[1].value, row[2], row[4].value
            # Accept prior site-specific labels without binding to a host name.
            if isinstance(label, str):
                label = re.sub(r"^FASTQ folder path \(on [^)]+\)$", "FASTQ folder path (on analysis host)", label)
                label = re.sub(r"^(Metadata file \(CSV or Excel\) — full path on ).+$", r"\1analysis host", label)
            question = by_id.get(field) if field else by_label.get(label)
            location = f"{sheet_name}!{cell.coordinate}"
            if question is None:
                if cell.value not in (None, "", "Your answer"):
                    errors.append(f"{location}: unrecognized question/field {field or label!r}")
                continue
            key = question["id"]
            if key in answers:
                errors.append(f"{location}: duplicate field {key}")
                continue
            value = cell.value
            if cell.data_type == "f":
                errors.append(f"{location}: enter a literal answer, not an Excel formula")
                continue
            if isinstance(value, str):
                value = value.strip()
            if value in (None, ""):
                if question.get("required") and key not in ("fastq_path", "metadata_file", "sample_id_column"):
                    errors.append(f"{location}: {question['text']} is required")
                value = None
            elif question["type"] == "dropdown" and value not in question["options"]:
                errors.append(f"{location}: choose one of {question['options']}")
            elif question["type"] == "number":
                if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 1 or int(value) != value:
                    errors.append(f"{location}: expected a positive whole number")
            elif not isinstance(value, (str, int, float)) or isinstance(value, bool):
                errors.append(f"{location}: expected text or a number; dates are not identifiers")
            answers[key] = value
            sources[key] = location
        for q in questions:
            if q.get("required") and q.get("introduced_in", 1) <= version and q["id"] not in answers:
                errors.append(f"{sheet_name}: missing required question {q['id']}")
        inline = answers.get("metadata_source") == "workbook tables"
        table_records = {}
        for name in list(metadata.SCHEMAS) + ["contrasts"]:
            title = name.title()
            if title not in workbook.sheetnames:
                if inline and name != "contrasts":
                    errors.append(f"Missing worksheet {title}")
                continue
            ws = workbook[title]
            headers = [c.value for c in ws[4]]
            while headers and headers[-1] is None:
                headers.pop()
            if not headers or any(not isinstance(c, str) or not c.strip() for c in headers) or len(set(headers)) != len(headers):
                errors.append(f"{title}!4: expected unique nonblank column names")
                continue
            rows = []
            for row in ws.iter_rows(min_row=5):
                if all(c.value in (None, "") for c in row):
                    continue
                values = {}
                for i, cell in enumerate(row):
                    if i >= len(headers):
                        if cell.value not in (None, ""):
                            errors.append(f"{title}!{cell.coordinate}: value has no column heading")
                        continue
                    key, value = headers[i], cell.value
                    location = f"{title}!{cell.coordinate}"
                    if cell.data_type == "f" or isinstance(value, bool) or (value is not None and not isinstance(value, (str, int, float))):
                        errors.append(f"{location}: enter literal text or a number, not formulas or dates")
                    if key.endswith("_id") and value is not None and not isinstance(value, str):
                        errors.append(f"{location}: identifiers must be stored as Text before entry")
                    values[key] = str(value).strip() if value is not None else ""
                    sources[f"{name}.{cell.row}.{key}"] = location
                rows.append(values)
            table_records[name] = rows
        contrast_rows = table_records.pop("contrasts", [])
        if inline:
            if answers.get("fastq_path") or answers.get("metadata_file"):
                errors.append("Workbook tables selected: clear FASTQ folder and external metadata file to avoid conflicting sources")
            errors.extend(f"{code}: {detail}" for code, detail in metadata.validate(table_records))
        else:
            if any(table_records.values()):
                errors.append("Workbook tables contain records: select workbook tables as metadata source or clear those records")
            for key in ("fastq_path", "metadata_file", "sample_id_column"):
                if not answers.get(key):
                    errors.append(f"{sources.get(key, sheet_name)}: {key} is required with external files")
        if errors:
            raise ValueError("Intake needs correction:\n  " + "\n  ".join(errors))
        return {
            "schema_version": version, "workbook_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "source": str(path), "sheet": sheet_name, "answers": answers, "sources": sources,
            "metadata_tables": table_records if inline else None,
            "contrast_rows": contrast_rows,
        }
    finally:
        workbook.close()


def setup_inputs(record):
    """Convert only executable questionnaire choices; never invent a model."""
    a, sources = record["answers"], record["sources"]

    def fail(key, detail):
        raise ValueError(f"{sources.get(key, record['sheet'])}: {detail}")

    if record["sheet"] != "Bulk RNA-seq":
        fail("project_name", "Direct intake execution currently supports Bulk RNA-seq only. "
             "TAG-seq needs kit qualification; small RNA needs its own processing route. "
             "No bulk fallback was selected.")
    if a.get("input_stage", "raw FASTQ") != "raw FASTQ":
        fail("input_stage", "Only raw FASTQ input is implemented; no counts/object or instrument-data import is available.")
    if a.get("umi", "no") != "no":
        fail("umi", "Confirm a non-UMI library. UMI processing is not implemented; unknown status must be resolved.")
    if not re.fullmatch(r"[A-Za-z0-9_]+", str(a["project_name"])):
        fail("project_name", "use letters, digits and underscores")
    layout = {"paired-end": "rnaseq_paired", "single-end": "rnaseq_single"}.get(a["read_config"])
    if layout is None:
        fail("read_config", "read configuration must be resolved before setup")
    factor = str(a["primary_factor"])
    batch = a["batch_effect_column"]
    subject = a.get("subject_column")
    subject_model = a.get("subject_model") or "random"
    if a["repeated_measures"] == "yes" and not subject:
        fail("subject_column", "Repeated measures need the subject column that identifies each biological unit")
    if subject and a["repeated_measures"] != "yes":
        fail("subject_column", "A subject column is only used with repeated measures = yes")
    terms = ([str(subject)] if subject and subject_model == "fixed block" else []) + \
            ([str(batch)] if batch and batch != "none" else []) + [factor]
    if any(not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", term) for term in terms + ([str(subject)] if subject else [])):
        fail("primary_factor", "factor/batch/subject column names must be simple identifiers")
    fixed = "~ " + " + ".join(dict.fromkeys(terms))
    if a.get("model_formula"):
        fixed = str(a["model_formula"]).strip()
        if not re.fullmatch(r"~\s*[A-Za-z][A-Za-z0-9_]*(\s*[+*:]\s*[A-Za-z][A-Za-z0-9_]*)*\s*", fixed):
            fail("model_formula", "use ~ followed by column names joined with +, * or :")
        missing = [term for term in terms if term not in re.findall(r"[A-Za-z][A-Za-z0-9_]*", fixed)]
        if missing:
            fail("model_formula", "formula must include the declared primary factor, batch and fixed subject block: " + ", ".join(missing))
    random = f"(1|{subject})" if subject and subject_model == "random" else None
    goal = a.get("analysis_goal") or "differential expression"
    objectives = {"differential expression": ["gene_expression", "differential_expression"],
                  "expression only": ["gene_expression"], "QC only": ["qc"]}[goal]
    organism = str(a["organism"])
    if organism.isdigit():
        tax_id = int(organism)
    elif isinstance(a["organism"], (int, float)) and math.isfinite(a["organism"]) and int(a["organism"]) == a["organism"]:
        tax_id = int(a["organism"])
    else:
        matches = []
        for p in (ROOT / "scripts/presets").glob("*.yaml"):
            preset = yaml.safe_load(p.read_text())
            if organism.casefold() in {str(preset.get(k, "")).casefold() for k in ("scientific_name", "common_name")}:
                matches.append(preset["tax_id"])
        if len(set(matches)) != 1:
            fail("organism", "supply an NCBI taxID; this name has no unambiguous local preset")
        tax_id = matches[0]
    if tax_id <= 0:
        fail("organism", "NCBI taxID must be positive")

    def absolute(value):
        if "://" in value:
            return value
        p = Path(value).expanduser()
        return str(p if p.is_absolute() else Path(record["source"]).parent / p)

    kits = yaml.safe_load((ROOT / "config/kit_profiles.yaml").read_text())
    protocol = a.get("library_protocol") or "explicit"
    kit = kits["profiles"].get(protocol, {})
    if protocol != "explicit":
        for key, field in (("strandedness", "strandedness"), ("adapter_r1", "preprocess_adapter_r1"), ("adapter_r2", "preprocess_adapter_r2")):
            given = a.get(field)
            if key != "strandedness" and given:
                given = str(given).upper()
            if given not in (None, "", "don't know") and given != kit[key]:
                fail(field, f"conflicts with {protocol} profile ({kit[key]})")
    strand = kit.get("strandedness") or a["strandedness"]
    expected = {"reverse": "ISR" if layout == "rnaseq_paired" else "SR",
                "forward": "ISF" if layout == "rnaseq_paired" else "SF",
                "unstranded": "IU" if layout == "rnaseq_paired" else "U"}.get(strand)
    supplied = a.get("expected_libtype")
    if supplied and expected and supplied.upper() != expected:
        fail("expected_libtype", f"conflicts with declared strand/layout ({expected})")
    if supplied:
        supplied = supplied.upper()
        allowed = {"IU", "ISF", "ISR", "OU", "OSF", "OSR", "MU", "MSF", "MSR"} if layout == "rnaseq_paired" else {"U", "SF", "SR"}
        if supplied not in allowed:
            fail("expected_libtype", f"invalid library type for {a['read_config']}")
    if a["orgdb_strategy"] == "use_existing":
        fail("orgdb_strategy", "This workbook has no OrgDb package field. Use YAML with an explicit package.")
    quantifier = a.get('analysis_quantifier') or 'auto'
    backend = a.get('analysis_backend') or 'auto'
    decoys = None if a.get('reference_decoys') == 'transcriptome only' else 'genome'
    custom_ref = {key: a.get(field) for key, field in (
        ('genome_fasta_url', 'reference_genome'), ('transcriptome_fasta_url', 'reference_transcriptome'),
        ('gtf_url', 'reference_gtf'))}
    if any(custom_ref.values()):
        required_refs = ['genome_fasta_url', 'gtf_url'] if quantifier == 'star' else ['transcriptome_fasta_url', 'gtf_url'] + (['genome_fasta_url'] if decoys else [])
        missing = [key for key in required_refs if not custom_ref[key]]
        if missing:
            fail('reference_gtf', 'Incomplete custom reference set: missing ' + ', '.join(missing))
        custom_ref = {key: (value if '://' in str(value) else Path(absolute(str(value))).as_uri())
                      if value else None for key, value in custom_ref.items()}
        custom_ref.update(accession='custom', assembly_name=None, gene_info_url=None)
    else:
        custom_ref = None
    if quantifier == 'star':
        tables = record.get('metadata_tables')
        if tables:
            if any(row['strandedness'] not in ('unstranded', 'forward', 'reverse') for row in tables['libraries']):
                fail('analysis_quantifier', 'STAR requires known strandedness for every library in Libraries')
        elif (supplied or expected) not in ({'IU', 'ISF', 'ISR'} if layout == 'rnaseq_paired' else {'U', 'SF', 'SR'}):
            fail('strandedness', 'STAR requires known strandedness and inward paired orientation')
    from preprocessing import policy_args
    poly_g = a.get('preprocess_poly_g') or 'auto'
    if poly_g == 'auto':
        poly_g = kits['platform_poly_g'].get(a.get('platform'), 'auto')
    paired = layout == 'rnaseq_paired'
    policy = dict(poly_g=poly_g,
                  adapter_r1=kit.get('adapter_r1') or str(a.get('preprocess_adapter_r1') or '').upper(),
                  adapter_r2=(kit.get('adapter_r2') if paired else '') or str(a.get('preprocess_adapter_r2') or '').upper(),
                  minimum_length=int(a.get('preprocess_minimum_length') or 15))
    try:
        policy, _ = policy_args(policy, layout == 'rnaseq_paired')
    except ValueError as exc:
        key = next((key for key in policy if key in str(exc)), 'minimum_length')
        fail('preprocess_' + key, str(exc))
    enrichment = a["run_enrichment"] == "yes"
    if enrichment and a["orgdb_strategy"] == "skip":
        fail("orgdb_strategy", "skip conflicts with requested GO + KEGG enrichment; select no enrichment or a strategy")
    if enrichment and goal != "differential expression":
        fail("run_enrichment", "enrichment needs the differential expression goal")
    if a["run_wgcna"] == "yes" and goal == "QC only":
        fail("run_wgcna", "WGCNA needs expression; QC only produces no counts")
    objectives += (["enrichment"] if enrichment else []) + (["coexpression"] if a["run_wgcna"] == "yes" else [])
    contrasts = []
    for row in record.get("contrast_rows") or []:
        cid = row.get("contrast_id", "")
        if not cid:
            fail("contrasts", "every Contrasts row needs a contrast_id")
        if row.get("weights"):
            weights = {}
            for part in row["weights"].split(";"):
                name, _, value = part.partition("=")
                try:
                    weights[name.strip()] = float(value)
                except ValueError:
                    fail("contrasts", f"Contrasts row {cid}: weights need coefficient=number pairs")
            if any(row.get(k) for k in ("factor", "by", "numerator", "denominator")):
                fail("contrasts", f"Contrasts row {cid}: use weights alone for a linear contrast")
            contrasts.append(dict(id=cid, type="linear", weights=weights))
        else:
            spec = dict(id=cid, type="pairwise", factor=row.get("factor", ""), by=row.get("by") or None)
            if row.get("numerator") or row.get("denominator"):
                spec.update(numerator=row.get("numerator", ""), denominator=row.get("denominator", ""))
            else:
                spec["reverse"] = True
            contrasts.append(spec)
    screening = {"enabled": False}
    if a.get("screening") == "yes":
        host = "genome" if quantifier == "star" else "transcriptome"
        refs = [dict(name="Host", role="expected", fasta=host),
                dict(name="PhiX", role="possible_contaminant", fasta=PHIX_URL)]
        for part in filter(None, (x.strip() for x in str(a.get("screening_extra") or "").split(";"))):
            name, _, source = part.partition("=")
            if not name.strip() or not source.strip():
                fail("screening_extra", "use name=FASTA pairs separated by ';'")
            refs.append(dict(name=name.strip(), role="possible_contaminant", fasta=absolute(source.strip())))
        screening = dict(enabled=True, fragments=int(a.get("screening_fragments") or 100000), seed=1, references=refs)
        from screen_reads import validate
        try:
            validate(screening)
        except ValueError as exc:
            fail("screening_extra", str(exc))
    record["context_only"] = ["contact_email", "tissue_type", "library_type", "library_kit"]
    record["limitations"] = [
        "A named protocol sets adapters and strand only; it is not a kit-level scientific qualification.",
        "Contact email is recorded but does not configure scheduler notifications.",
        "Multiple libraries per sample and non-bulk workbook execution remain pending assay qualification.",
    ]
    if record["schema_version"] == 1:
        record["limitations"].append("Legacy workbook: assumes raw non-UMI bulk input. Confirm protocol or use the current template with explicit eligibility fields.")
    tables = record.get("metadata_tables")
    if tables:
        units = [row.get("biological_unit_id") or row["sample_id"] for row in tables["samples"]]
        if len(set(units)) != len(units) and not subject:
            fail("repeated_measures", "Repeated biological units require repeated measures = yes and a subject column")
        for read in tables["reads"]:
            read["uri"] = absolute(read["uri"])
        if kit and any(lib["strandedness"] not in ("unknown", kit["strandedness"]) for lib in tables["libraries"]):
            fail("library_protocol", f"Libraries strandedness conflicts with {protocol} ({kit['strandedness']})")
    return dict(project_name=a["project_name"], tax_id=tax_id, seq_type=layout,
                analysis=dict(quantifier=quantifier, backend=backend), preprocessing=policy,
                reference_overrides=custom_ref, reference_decoys=decoys,
                metadata_tables=tables,
                fastq_source=absolute(str(a["fastq_path"])) if not tables else None,
                metadata_file=absolute(str(a["metadata_file"])) if not tables else None,
                sample_id_column=str(a["sample_id_column"]), primary_factor=factor,
                model_fixed_effects=fixed, model_random_effects=random, contrasts=contrasts or None,
                objectives=objectives, screening=screening,
                expected_libtype=supplied or expected, description=a["biological_question"],
                run_enrichment=enrichment, run_wgcna=a["run_wgcna"] == "yes", subject_column=subject,
                orgdb_strategy=a["orgdb_strategy"] if enrichment else "skip",
                orgdb_cache_dir=absolute(str(a["orgdb_cache_dir"])) if a.get("orgdb_cache_dir") else None,
                expected_min_samples=a["min_samples_per_group"], storage_budget_gb=None)


# First stage each input affects; Snakemake's own rerun identity decides what
# actually reruns. project_name changes the result root, so it needs a new project.
IMPACT = [
    (("tax_id", "reference_overrides", "reference_decoys"), "reference, quantification, counts, statistics, report"),
    (("seq_type", "expected_libtype", "preprocessing", "metadata_tables", "fastq_source",
      "metadata_file", "sample_id_column"), "preprocessing, screening, quantification, counts, statistics, report"),
    (("analysis.quantifier",), "quantification, counts, statistics, report"),
    (("screening",), "screening, report"),
    (("objectives",), "requested targets"),
    (("analysis.backend", "model_fixed_effects", "model_random_effects", "primary_factor", "contrasts",
      "subject_column", "expected_min_samples"), "differential expression, enrichment, report"),
    (("run_enrichment", "orgdb_strategy", "orgdb_cache_dir"), "enrichment, report"),
    (("run_wgcna",), "coexpression, report"),
]


def revision_impact(old, new):
    """Changed setup inputs and the stages they reach."""
    def flat(d):
        return {f"{k}.{s}": v for k, sub in d.items() if isinstance(sub, dict) and k in ("analysis",)
                for s, v in sub.items()} | {k: v for k, v in d.items() if k != "analysis"}
    a, b = flat(old), flat(new)
    changes = []
    for key in sorted(set(a) | set(b)):
        if a.get(key) != b.get(key):
            stages = next((s for keys, s in IMPACT if key in keys), "report")
            changes.append(dict(field=key, stages=stages))
    blocked = [c["field"] for c in changes if c["field"] == "project_name"]
    return dict(changes=changes, blocked=blocked)
