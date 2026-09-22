#!/usr/bin/env python3
"""qc_report.py — auto-generate comparative QC + alignment charts.

Takes the sample list from --samples (the run's sample sheet), so sample
directories left behind by an earlier run on a different cohort are ignored.
The report title comes from --project/--species. The reference line comes from
--reference-lock (reference.lock.json) and from nothing else: R17c requires the
actual construction to be reported from provenance, and RF02 forbids describing
a transcriptome-only index as decoy-aware.

For each sample it reads, under <quant>/<sample>/:
  fastp.json                  (cleaning: before vs after filtering)
  aux_info/meta_info.json     (salmon alignment stats)
  lib_format_counts.json      (detected library/strand type)

and writes, into <out>/:
  alignment_summary.tsv       per-sample table; missing values are "NA"
  comparative_charts.png      before/after cleaning + alignment-fate charts
  qc_charts.html              self-contained report embedding the PNG + table

A metric that is absent is carried as None all the way through: it is written
as "NA", and no bar is drawn for it. It is never coerced to 0, which would be
indistinguishable from a sample that genuinely had zero reads.

matplotlib is imported inside make_charts so the collection and table logic can
be imported and checked without it (see tests/check_qc_report.py).
"""
import argparse, base64, json, os


def load_json(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except Exception:
        return None


def read_sample_ids(sheet_path):
    """Sample IDs from the run's sample sheet, in sheet order.

    This is the authoritative list. Globbing the quant directory instead would
    pick up directories from a previous run on a different cohort.
    """
    ids = []
    with open(sheet_path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                if "sample_id" not in header:
                    raise SystemExit(
                        f"qc_report: {sheet_path} has no 'sample_id' column")
                continue
            ids.append(fields[header.index("sample_id")])
    if not ids:
        raise SystemExit(f"qc_report: no samples listed in {sheet_path}")
    return ids


def collect(quant_dir, sample_ids):
    """One record per sample in sample_ids. Absent metrics stay None."""
    recs = []
    for s in sample_ids:
        sd = os.path.join(quant_dir, s)
        fp = load_json(os.path.join(sd, "fastp.json"))
        mi = load_json(os.path.join(sd, "aux_info", "meta_info.json"))
        lf = load_json(os.path.join(sd, "lib_format_counts.json"))
        r = {"sample": s, "has_quant": os.path.isfile(os.path.join(sd, "quant.sf"))}
        for c in COLS[1:]:
            r[c] = None
        if fp:
            b = fp["summary"]["before_filtering"]
            a = fp["summary"]["after_filtering"]
            ac = fp.get("adapter_cutting", {})
            r["reads_before"] = b["total_reads"]
            r["reads_after"] = a["total_reads"]
            r["q30_before"] = round(b["q30_rate"] * 100, 2)
            r["q30_after"] = round(a["q30_rate"] * 100, 2)
            # Absent only when no adapters were looked for; 0 here is a real 0.
            r["adapter_trimmed_reads"] = ac.get("adapter_trimmed_reads", 0)
        if mi:
            r["num_processed"] = mi.get("num_processed")
            r["num_mapped"] = mi.get("num_mapped")
            r["num_decoy"] = mi.get("num_decoy_fragments")
            pm = mi.get("percent_mapped")
            r["percent_mapped"] = round(pm, 2) if pm is not None else None
            fl = mi.get("frag_length_mean")
            r["frag_len_mean"] = round(fl, 1) if fl is not None else None
        if lf:
            r["library_type"] = lf.get("expected_format")
        recs.append(r)
    return recs


COLS = ["sample", "percent_mapped", "num_mapped", "num_processed", "num_decoy",
        "frag_len_mean", "library_type", "reads_before", "reads_after",
        "q30_before", "q30_after", "adapter_trimmed_reads"]


def cell(v):
    """Table cell text. None becomes NA so it cannot be read as a zero."""
    return "NA" if v is None else str(v)


def write_tsv(recs, out):
    with open(os.path.join(out, "alignment_summary.tsv"), "w") as o:
        o.write("\t".join(COLS) + "\n")
        for r in recs:
            o.write("\t".join(cell(r.get(c)) for c in COLS) + "\n")


def present(recs, key, scale=1.0):
    """(x positions, heights) for records that actually have `key`.

    Records missing it are omitted rather than plotted at 0, so a gap in the
    chart means "not measured" and a zero-height bar means "measured as zero".
    """
    xs, ys = [], []
    for i, r in enumerate(recs):
        v = r.get(key)
        if v is not None:
            xs.append(i)
            ys.append(v / scale)
    return xs, ys


def unmapped(recs):
    """(x, height) for processed - mapped - decoy, only where all three exist."""
    xs, ys = [], []
    for i, r in enumerate(recs):
        p, m, d = r.get("num_processed"), r.get("num_mapped"), r.get("num_decoy")
        if None not in (p, m, d):
            xs.append(i)
            ys.append(max(0, p - m - d) / 1e6)
    return xs, ys


def tick_labels(recs):
    """Sample names, marking the ones with no metrics so gaps are readable."""
    out = []
    for r in recs:
        missing = []
        if r.get("reads_before") is None:
            missing.append("no fastp")
        if r.get("num_processed") is None:
            missing.append("no salmon")
        out.append(r["sample"] + (f" ({', '.join(missing)})" if missing else ""))
    return out


def make_charts(recs, out, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(recs)
    x = list(range(n))
    labels = tick_labels(recs)
    w = 0.4
    fig, axes = plt.subplots(1, 3, figsize=(max(10, n * 1.1 + 6), 5))

    def style(ax, ttl, ylab):
        ax.set_title(ttl)
        ax.set_ylabel(ylab)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.legend(fontsize=8)

    # 1) reads before vs after cleaning
    ax = axes[0]
    xb, rb = present(recs, "reads_before", 1e6)
    xa, ra = present(recs, "reads_after", 1e6)
    ax.bar([i - w / 2 for i in xb], rb, w, label="before clean", color="#9ecae1")
    ax.bar([i + w / 2 for i in xa], ra, w, label="after clean", color="#3182bd")
    style(ax, "Reads: before vs after cleaning", "reads, R1+R2 (millions)")

    # 2) Q30 before vs after cleaning
    ax = axes[1]
    xqb, qb = present(recs, "q30_before")
    xqa, qa = present(recs, "q30_after")
    ax.bar([i - w / 2 for i in xqb], qb, w, label="before", color="#a1d99b")
    ax.bar([i + w / 2 for i in xqa], qa, w, label="after", color="#31a354")
    ax.set_ylim(0, 100)
    style(ax, "Q30% before vs after cleaning", "% bases >= Q30")

    # 3) alignment fate: mapped / decoy / unmapped (stacked)
    # Each segment is drawn only where its own count exists; a stacked segment
    # sits on whichever lower segments exist for that sample.
    ax = axes[2]
    xm, mapped = present(recs, "num_mapped", 1e6)
    xd, decoy = present(recs, "num_decoy", 1e6)
    xu, unmap = unmapped(recs)
    mapped_at = dict(zip(xm, mapped))
    decoy_at = dict(zip(xd, decoy))
    ax.bar(xm, mapped, label="mapped to genes", color="#31a354")
    ax.bar(xd, decoy, bottom=[mapped_at.get(i, 0) for i in xd],
           label="genome / decoy", color="#fec44f")
    ax.bar(xu, unmap,
           bottom=[mapped_at.get(i, 0) + decoy_at.get(i, 0) for i in xu],
           label="unmapped", color="#de2d26")
    style(ax, "Alignment fate (fragments)", "read pairs (millions)")

    fig.suptitle(title, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    png = os.path.join(out, "comparative_charts.png")
    fig.savefig(png, dpi=120)
    plt.close(fig)
    return png


def reference_description(lock_path):
    """How the reference was actually built, read from reference.lock.json.

    R17c requires the report to state the actual reference construction from
    provenance. RF02 forbids calling a transcriptome-only index decoy-aware,
    and the previous version of this file printed "decoy-aware salmon"
    unconditionally while retrieve.smk built no decoys at all.

    When no lock is available the construction is reported as unrecorded. It is
    never guessed, and no decoy claim is made in that case.
    """
    if not lock_path or not os.path.isfile(lock_path):
        return None, "reference construction unrecorded (no reference.lock.json)"

    lock = load_json(lock_path)
    if not lock:
        return None, "reference construction unrecorded (lock unreadable)"

    asm = lock.get("assembly") or {}
    idx = lock.get("index") or {}
    decoys = idx.get("decoys") or {}

    acc = asm.get("accession") or "unspecified accession"
    name = asm.get("assembly_name")
    ref = f"{acc} ({name})" if name else acc

    tool = idx.get("tool") or "unknown quantifier"
    ver = idx.get("tool_version")
    tool_s = f"{tool} {ver}" if ver else tool
    mode = idx.get("quantification_mode") or "unrecorded mode"
    label = decoys.get("label") or "decoy status unrecorded"

    return lock, f"{ref}, {tool_s}, {mode}, {label}"


def identity(project, species):
    """Normalise the run's identity strings, filling blanks with honest text."""
    name = project.strip() or "RNA-Seq run"
    species = species.strip()
    title = f"{name} — {species}" if species else name
    return name, species, f"{title}: cleaning & alignment QC"


def write_html(recs, png, out, name, species, reference):
    b64 = base64.b64encode(open(png, "rb").read()).decode()
    rows = "\n".join(
        "<tr>" + "".join(f"<td>{cell(r.get(c))}</td>" for c in COLS) + "</tr>"
        for r in recs)
    head = "".join(f"<th>{c}</th>" for c in COLS)
    heading = f"{name} — <i>{species}</i>" if species else name
    n_no_fastp = sum(1 for r in recs if r.get("reads_before") is None)
    n_no_salmon = sum(1 for r in recs if r.get("num_processed") is None)
    gaps = ""
    if n_no_fastp or n_no_salmon:
        gaps = (f" <b>{n_no_fastp} sample(s) without fastp metrics and "
                f"{n_no_salmon} without salmon metrics are shown as NA and "
                f"left blank in the charts, not as zero.</b>")
    html = f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<title>{name} QC &amp; Alignment Report</title>
<style>
 body{{font-family:-apple-system,Arial,sans-serif;margin:30px;color:#222}}
 h1{{font-size:20px}} h2{{font-size:16px;margin-top:28px}}
 img{{max-width:100%;border:1px solid #ddd}}
 table{{border-collapse:collapse;margin-top:10px;font-size:13px}}
 th,td{{border:1px solid #ccc;padding:5px 9px;text-align:right}}
 th:first-child,td:first-child{{text-align:left}} th{{background:#f4f4f4}}
 .note{{color:#666;font-size:13px}}
</style></head><body>
<h1>{heading}: cleaning &amp; alignment QC</h1>
<p class="note">{len(recs)} sample(s). Reference: {reference}.
Auto-generated by the pipeline (rule qc_report).{gaps}</p>
<h2>Comparative charts</h2>
<img src="data:image/png;base64,{b64}" alt="comparative charts">
<h2>Per-sample summary</h2>
<table><tr>{head}</tr>
{rows}
</table>
<p class="note">Full interactive QC: see <code>multiqc/multiqc_report.html</code> in this folder.</p>
</body></html>"""
    with open(os.path.join(out, "qc_charts.html"), "w") as o:
        o.write(html)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quant", required=True, help="quant/ directory (per-sample subdirs)")
    ap.add_argument("--out", required=True, help="output directory for the report")
    ap.add_argument("--samples", required=True,
                    help="sample sheet; its sample_id column is the authoritative list")
    ap.add_argument("--project", default="", help="project name for the report title")
    ap.add_argument("--species", default="", help="organism name for the report title")
    ap.add_argument("--reference-lock", default="",
                    help="reference.lock.json; the only permitted source for "
                         "statements about how the reference was built")
    a = ap.parse_args()

    sample_ids = read_sample_ids(a.samples)
    name, species, title = identity(a.project, a.species)
    _lock, reference = reference_description(a.reference_lock)

    os.makedirs(a.out, exist_ok=True)
    recs = collect(a.quant, sample_ids)

    missing = [r["sample"] for r in recs if not r["has_quant"]]
    if len(missing) == len(recs):
        raise SystemExit(
            f"qc_report: none of the {len(recs)} samples in the sheet have "
            f"quant.sf under {a.quant}")
    if missing:
        print(f"qc_report: WARNING {len(missing)} sample(s) in the sheet have no "
              f"quant.sf and are reported as NA: {', '.join(missing)}")

    stale = sorted(set(os.listdir(a.quant)) - set(sample_ids)) if os.path.isdir(a.quant) else []
    if stale:
        print(f"qc_report: ignoring {len(stale)} directory/ies not in the sample "
              f"sheet: {', '.join(stale)}")

    write_tsv(recs, a.out)
    png = make_charts(recs, a.out, title)
    write_html(recs, png, a.out, name, species, reference)
    print(f"qc_report: wrote summary + charts for {len(recs)} sample(s) to {a.out}")


if __name__ == "__main__":
    main()
