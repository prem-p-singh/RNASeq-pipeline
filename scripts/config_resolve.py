#!/usr/bin/env python3
"""config_resolve.py — one configuration resolution path (R17e, R05).

Master plan section 6.3 states the contract: "Resolve schema defaults,
compatible assay/organism presets, project choices, analysis choices, and
execution settings into one immutable snapshot. Record the origin of every
nontrivial decision. Deep-merge supported mappings; replace lists explicitly;
reject unknown keys and incompatible types."

Before this module there were three separate paths and they disagreed:

  Snakefile     merged config/thresholds.yaml one level deep
  preflight.py  checked keys and types against config.template.yaml
  setup.py      rendered a project config from the template

The one-level merge had a live consequence. config.template.yaml carried its
own `thresholds.wgcna` block that omitted merge_cut_height_default and
merge_cut_height_bumped. Because `thresholds.wgcna` was already present, the
richer block in config/thresholds.yaml was never merged, so 05_wgcna.R read
NULL for both heights. WGCNA does not fail on that: it catches the error and
returns unmerged modules, so the run completed while the requested merging
silently never happened.

Defaults therefore have a single source (master plan 5.2): config.template.yaml
for the main configuration, config/thresholds.yaml for the `thresholds`
subtree. This module composes them and is the only place that does.

Pure functions, no I/O beyond reading the two schema files, so setup, preflight
and the Snakefile can all call it and cannot drift apart.
"""
from __future__ import annotations

from pathlib import Path

import yaml

# Subtrees whose keys are chosen by the user, not fixed by the schema. Their
# contents are carried through untouched and never reported as unknown.
OPAQUE: tuple[tuple[str, ...], ...] = (
    ("contrasts",),
    ("reference", "annotation_tsv", "columns"),
    ("samples", "columns"),
)

# Injected by submit.sh at launch rather than authored by a user.
INJECTED: tuple[tuple[str, ...], ...] = (("repo_dir",),)


def _is_opaque(path: tuple[str, ...]) -> bool:
    return any(path[: len(o)] == o for o in OPAQUE + INJECTED)


def load_schema(repo_root) -> dict:
    """The composed default/schema document.

    config.template.yaml supplies everything except `thresholds`, which comes
    from config/thresholds.yaml. A `thresholds` block in the template is
    ignored deliberately: two sources of the same defaults is what produced the
    silent WGCNA shadowing, and master plan 5.2 requires a single source.
    """
    root = Path(repo_root)
    template = yaml.safe_load((root / "config" / "config.template.yaml").read_text()) or {}
    thresholds = yaml.safe_load((root / "config" / "thresholds.yaml").read_text()) or {}
    schema = {k: v for k, v in template.items() if k != "thresholds"}
    schema["thresholds"] = thresholds
    return schema


def deep_merge(base, override):
    """Mappings merge recursively; everything else replaces.

    Lists replace rather than concatenate, per section 6.3 ("replace lists
    explicitly"): appending to a default list would silently keep entries the
    user meant to drop.
    """
    if not isinstance(base, dict) or not isinstance(override, dict):
        return override
    out = dict(base)
    for k, v in override.items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def walk(node, prefix=()):
    """Yield (path, value) for every key in a nested mapping."""
    if isinstance(node, dict):
        for k, v in node.items():
            path = prefix + (k,)
            yield path, v
            if isinstance(v, dict) and not _is_opaque(path):
                yield from walk(v, path)


def unknown_keys(cfg: dict, schema: dict) -> list[str]:
    """Dotted paths present in cfg but absent from the schema."""
    known = {p for p, _ in walk(schema)}
    out = []
    for path, _ in walk(cfg):
        if _is_opaque(path) or path in known:
            continue
        # Report the shallowest unknown path only, so one stray section does
        # not produce a finding per leaf beneath it.
        if any(path[:i] not in known for i in range(1, len(path))):
            continue
        out.append(".".join(path))
    return sorted(out)


def _kind(v):
    if isinstance(v, bool):
        return "boolean"
    if isinstance(v, (int, float)):
        return "number"
    if isinstance(v, str):
        return "string"
    if isinstance(v, list):
        return "list"
    if isinstance(v, dict):
        return "mapping"
    return None


def type_conflicts(cfg: dict, schema: dict) -> list[str]:
    """Paths whose value kind disagrees with the schema's.

    Only checked where the schema commits to a concrete value: a null default
    declares the key without constraining its type.
    """
    by_path = {p: v for p, v in walk(schema)}
    out = []
    for path, value in walk(cfg):
        if _is_opaque(path) or path not in by_path:
            continue
        want, got = _kind(by_path[path]), _kind(value)
        if want is None or got is None or by_path[path] is None or value is None:
            continue
        if want != got:
            out.append(f"{'.'.join(path)}: expected {want}, got {got}")
    return sorted(out)


# Keys kept in the schema only so a legacy project is not rejected outright.
# RS12: a legacy storage cap must raise a migration warning, not an error.
DEPRECATED = {("hpc", "storage_budget_gb"): "CFG005"}


def resolve(repo_root, layers):
    """Merge ordered layers over the schema defaults.

    `layers` is [(origin_name, mapping), ...], applied left to right. Returns
    (config, origins, issues) where origins maps each dotted path to the name
    of the layer that last set it, satisfying "record the origin of every
    nontrivial decision", and issues is a list of (code, detail) using the
    catalogue in config/spec.yaml.
    """
    schema = load_schema(repo_root)
    config = schema
    origins = {".".join(p): "default" for p, _ in walk(schema)}

    issues = []
    for name, layer in layers:
        if not layer:
            continue
        for code, details in (("CFG001", unknown_keys(layer, schema)),
                              ("CFG003", type_conflicts(layer, schema))):
            for d in details:
                issues.append((code, f"{d} (from {name})"))
        for path, value in walk(layer):
            code = DEPRECATED.get(path)
            if code and value not in (None, ""):
                issues.append((code, f"{'.'.join(path)}={value!r} (from {name})"))
        config = deep_merge(config, layer)
        for p, _ in walk(layer):
            if not _is_opaque(p):
                origins[".".join(p)] = name

    return config, origins, issues


def blocking(issues, repo_root) -> list:
    """Issues whose severity in config/spec.yaml is `error`.

    A deprecation is a warning and must not stop a run: RS12 requires a
    migration warning for a legacy storage cap, not a rejection.
    """
    try:
        catalog = yaml.safe_load(
            (Path(repo_root) / "config" / "spec.yaml").read_text())["issues"]
    except (OSError, KeyError, ValueError):
        return list(issues)          # cannot classify: treat all as blocking
    return [(c, d) for c, d in issues
            if catalog.get(c, {}).get("severity", "error") == "error"]


def demo():
    """Self-check: the shadowing this module exists to prevent."""
    schema = {"thresholds": {"wgcna": {"min_samples": 15, "target_r2": 0.85,
                                       "merge_cut_height_default": 0.25}}}
    user = {"thresholds": {"wgcna": {"min_samples": 20}}}
    merged = deep_merge(schema, user)
    w = merged["thresholds"]["wgcna"]
    assert w["min_samples"] == 20, w
    assert w["merge_cut_height_default"] == 0.25, "sibling default was dropped"
    assert w["target_r2"] == 0.85, w

    # lists replace, they do not append
    assert deep_merge({"a": [1, 2]}, {"a": [3]})["a"] == [3]

    assert unknown_keys({"bogus": 1}, schema) == ["bogus"]
    assert unknown_keys({"thresholds": {"wgcna": {"min_samples": 20}}}, schema) == []
    assert type_conflicts({"thresholds": {"wgcna": {"min_samples": "many"}}}, schema)
    assert type_conflicts({"thresholds": {"wgcna": {"min_samples": 20}}}, schema) == []
    print("config_resolve.demo: ok")


if __name__ == "__main__":
    demo()
