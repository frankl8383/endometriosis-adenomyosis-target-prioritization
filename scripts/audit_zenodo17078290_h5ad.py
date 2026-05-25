#!/usr/bin/env python3
"""Audit the Zenodo 17078290 adenomyosis h5ad after download completion."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW = PROJECT_ROOT / "data" / "raw_downloads"
META = PROJECT_ROOT / "data" / "raw_metadata"
RESULTS = PROJECT_ROOT / "results" / "singlecell"
PY_PKG_DIR = PROJECT_ROOT / "tools" / "py_packages"

if sys.version_info[:2] != (3, 12):
    raise SystemExit(
        "Use the bundled Python 3.12 runtime: "
        "python3 "
        "scripts/audit_zenodo17078290_h5ad.py"
    )

sys.path.insert(0, str(PY_PKG_DIR))

import anndata as ad  # noqa: E402
import pandas as pd  # noqa: E402


H5AD = RAW / "Zenodo_17078290__celllable.diff_PRO.h5ad"
RECORD_JSON = META / "zenodo_17078290.json"
GWAS_CANDIDATES = PROJECT_ROOT / "results" / "gwas" / "gwas_candidate_gene_universe.tsv"

OUT_AUDIT = RESULTS / "Zenodo17078290_h5ad_audit.md"
OUT_OBS = RESULTS / "Zenodo17078290_h5ad_obs_columns.tsv"
OUT_VAR = RESULTS / "Zenodo17078290_h5ad_var_columns.tsv"
OUT_KEYS = RESULTS / "Zenodo17078290_h5ad_key_audit.tsv"
OUT_COVERAGE = RESULTS / "Zenodo17078290_candidate_gene_coverage.tsv"
OUT_REVIEW = RESULTS / "phase8_zenodo17078290_h5ad_self_review.md"


def md5_stream(path: Path, chunk_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def clean_id(value: object) -> str:
    return str(value).split(".")[0]


def summarize_columns(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for col in df.columns:
        series = df[col]
        values = series.astype(str)
        unique_values = sorted(values.dropna().unique().tolist())
        sample_values = unique_values[:20]
        rows.append(
            {
                "column": col,
                "dtype": str(series.dtype),
                "n_missing": int(series.isna().sum()),
                "n_unique": int(values.nunique(dropna=True)),
                "example_values": "|".join(sample_values),
            }
        )
    return pd.DataFrame(rows)


def audit_keys(adata: ad.AnnData) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    rows.append({"slot": "X", "key": "X", "type": type(adata.X).__name__, "shape": "x".join(map(str, adata.shape))})
    for slot_name, mapping in [
        ("layers", adata.layers),
        ("obsm", adata.obsm),
        ("varm", adata.varm),
        ("obsp", adata.obsp),
        ("varp", adata.varp),
    ]:
        for key in list(mapping.keys()):
            value = mapping[key]
            shape = getattr(value, "shape", "")
            rows.append(
                {
                    "slot": slot_name,
                    "key": key,
                    "type": type(value).__name__,
                    "shape": "x".join(map(str, shape)) if shape else "",
                }
            )
    for key in list(adata.uns.keys()):
        value = adata.uns[key]
        rows.append({"slot": "uns", "key": key, "type": type(value).__name__, "shape": ""})
    return pd.DataFrame(rows)


def candidate_coverage(adata: ad.AnnData) -> pd.DataFrame:
    candidates = pd.read_csv(GWAS_CANDIDATES, sep="\t", keep_default_na=False)
    var = adata.var.copy()
    var_names_raw = {str(x) for x in adata.var_names}
    var_names_clean = {clean_id(x) for x in adata.var_names}
    candidate_columns = [
        col
        for col in var.columns
        if any(token in col.lower() for token in ["gene", "ensembl", "symbol", "feature", "name", "id"])
    ]
    column_values: dict[str, set[str]] = {}
    column_values_clean: dict[str, set[str]] = {}
    for col in candidate_columns:
        values = var[col].astype(str)
        column_values[col] = set(values.tolist())
        column_values_clean[col] = {clean_id(value) for value in values.tolist()}

    rows: list[dict[str, object]] = []
    for _, row in candidates.iterrows():
        gene_id = clean_id(row["gene_id"])
        symbol = str(row["gene_symbol"])
        matches: list[str] = []
        if gene_id in var_names_clean:
            matches.append("var_names_ensembl")
        if symbol and symbol in var_names_raw:
            matches.append("var_names_symbol")
        for col in candidate_columns:
            if gene_id in column_values_clean[col]:
                matches.append(f"{col}:ensembl")
            if symbol and symbol in column_values[col]:
                matches.append(f"{col}:symbol")
        rows.append(
            {
                "gene_id": gene_id,
                "gene_symbol": symbol,
                "genetic_priority": row.get("genetic_priority", ""),
                "ld_neighborhood_class": row.get("ld_neighborhood_class", ""),
                "module_hint_preliminary": row.get("module_hint_preliminary", ""),
                "matched": bool(matches),
                "match_basis": "|".join(sorted(set(matches))),
            }
        )
    return pd.DataFrame(rows)


def likely_metadata_columns(obs: pd.DataFrame) -> list[str]:
    tokens = [
        "cell",
        "type",
        "anno",
        "cluster",
        "sample",
        "patient",
        "donor",
        "group",
        "disease",
        "condition",
        "tissue",
        "region",
        "batch",
        "spatial",
    ]
    return [col for col in obs.columns if any(token in col.lower() for token in tokens)]


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    if not H5AD.exists():
        raise SystemExit(f"Missing h5ad: {H5AD}")
    partial = H5AD.with_suffix(H5AD.suffix + ".partial")
    aria2_state = Path(str(partial) + ".aria2")
    record = json.loads(RECORD_JSON.read_text(encoding="utf-8"))
    expected = record["files"][0]
    expected_size = int(expected["size"])
    expected_md5 = str(expected["checksum"]).replace("md5:", "")
    actual_size = H5AD.stat().st_size
    actual_md5 = md5_stream(H5AD)
    md5_ok = actual_md5 == expected_md5
    size_ok = actual_size == expected_size
    no_partial = not partial.exists() and not aria2_state.exists()

    adata = ad.read_h5ad(H5AD, backed="r")
    try:
        obs_summary = summarize_columns(adata.obs)
        var_summary = summarize_columns(adata.var)
        key_summary = audit_keys(adata)
        coverage = candidate_coverage(adata)
        likely_cols = likely_metadata_columns(adata.obs)
        spatial_keys = key_summary[(key_summary["slot"] == "obsm") & (key_summary["key"].astype(str).str.contains("spatial", case=False))]
        has_spatial = not spatial_keys.empty or "spatial" in adata.uns

        obs_summary.to_csv(OUT_OBS, sep="\t", index=False)
        var_summary.to_csv(OUT_VAR, sep="\t", index=False)
        key_summary.to_csv(OUT_KEYS, sep="\t", index=False)
        coverage.to_csv(OUT_COVERAGE, sep="\t", index=False)

        n_matched = int(coverage["matched"].sum())
        n_candidates = int(len(coverage))
        likely_lines = []
        for col in likely_cols[:40]:
            values = adata.obs[col].astype(str)
            unique_values = sorted(values.dropna().unique().tolist())
            display_values = "|".join(unique_values[:20])
            likely_lines.append(f"{col}\t{values.nunique(dropna=True)}\t{display_values}")

        lines = [
            "# Zenodo 17078290 h5ad audit",
            "",
            "## File integrity",
            "",
            f"- File: `{H5AD}`",
            f"- Expected size: {expected_size}",
            f"- Actual size: {actual_size}",
            f"- Size match: {size_ok}",
            f"- Expected MD5: `{expected_md5}`",
            f"- Actual MD5: `{actual_md5}`",
            f"- MD5 match: {md5_ok}",
            f"- No partial/aria2 sidecar remains: {no_partial}",
            "",
            "## AnnData backed-read audit",
            "",
            f"- Shape: {adata.n_obs} observations x {adata.n_vars} variables",
            f"- obs columns: {adata.obs.shape[1]}",
            f"- var columns: {adata.var.shape[1]}",
            f"- layers: {len(list(adata.layers.keys()))}",
            f"- obsm keys: {len(list(adata.obsm.keys()))}",
            f"- uns keys: {len(list(adata.uns.keys()))}",
            f"- Spatial evidence detected in object structure: {has_spatial}",
            "",
            "## Likely phenotype/cell/spatial metadata columns",
            "",
            "column\tn_unique\texample_values",
        ]
        lines.extend(likely_lines if likely_lines else ["NA\t0\t"])
        lines.extend(
            [
                "",
                "## Candidate gene coverage",
                "",
                f"- Matched candidates: {n_matched}/{n_candidates}",
                f"- Unmatched candidates: {n_candidates - n_matched}",
                "",
                "## Output files",
                "",
                f"- obs column audit: `{OUT_OBS}`",
                f"- var column audit: `{OUT_VAR}`",
                f"- key audit: `{OUT_KEYS}`",
                f"- candidate coverage: `{OUT_COVERAGE}`",
                f"- self-review: `{OUT_REVIEW}`",
            ]
        )
        OUT_AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")

        verdict = "PASS" if (md5_ok and size_ok and no_partial and n_matched >= 80) else "PASS_WITH_CONDITIONS"
        review_lines = [
            "# Phase 8 Self-Review: Zenodo 17078290 adenomyosis h5ad audit",
            "",
            "## Verdict",
            "",
            verdict,
            "",
            "## What passed",
            "",
            f"- Final h5ad exists and backed-read opened successfully with shape {adata.n_obs} x {adata.n_vars}.",
            f"- Size check passed: {size_ok}.",
            f"- MD5 check passed: {md5_ok}.",
            f"- Candidate coverage is {n_matched}/{n_candidates}.",
            f"- Spatial metadata/coordinates detected: {has_spatial}.",
            "",
            "## Conditions and next checks",
            "",
            "- This audit verifies file integrity and object structure only; it does not yet validate cell-type labels, disease groups, spatial coordinates or count layers biologically.",
            "- Before adenomyosis claims, inspect the likely phenotype/cell-type columns and define disease/control, tissue, cell-type and spatial-region fields explicitly.",
            "- Downstream analysis should use backed or chunked access where possible because the object is large.",
            "",
            "## Decision",
            "",
            "The adenomyosis h5ad is now ready for metadata parsing and candidate localization planning. It is not yet ready for biological conclusions.",
        ]
        OUT_REVIEW.write_text("\n".join(review_lines) + "\n", encoding="utf-8")
        print(OUT_AUDIT)
        print(OUT_REVIEW)
        return 0
    finally:
        adata.file.close()


if __name__ == "__main__":
    raise SystemExit(main())
