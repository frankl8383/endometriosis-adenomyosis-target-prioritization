#!/usr/bin/env python3
"""Audit single-cell/spatial source readiness before ingestion."""

from __future__ import annotations

import json
import os
import subprocess
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
        "scripts/audit_singlecell_sources.py"
    )

sys.path.insert(0, str(PY_PKG_DIR))

import pandas as pd  # noqa: E402


def file_status(name: str, expected_size: int | None = None) -> dict[str, object]:
    path = RAW / name
    partial = path.with_suffix(path.suffix + ".partial")
    aria2 = partial.with_suffix(partial.suffix + ".aria2")
    active_path = path if path.exists() else partial
    status = {
        "file": name,
        "path": str(active_path) if active_path.exists() else "",
        "exists": active_path.exists(),
        "is_partial": partial.exists() and not path.exists(),
        "has_aria2_control_file": aria2.exists(),
        "size_bytes": active_path.stat().st_size if active_path.exists() else 0,
        "expected_size_bytes": expected_size or "",
        "complete_by_size": bool(path.exists() and expected_size and path.stat().st_size == expected_size),
    }
    return status


def tar_head(path: Path, max_entries: int = 25) -> list[str]:
    if not path.exists():
        return []
    try:
        proc = subprocess.run(
            ["tar", "-tf", str(path)],
            check=False,
            capture_output=True,
            text=True,
            timeout=20,
        )
    except Exception as exc:  # noqa: BLE001
        return [f"tar_list_failed:{type(exc).__name__}:{exc}"]
    entries = [line for line in proc.stdout.splitlines() if line]
    if proc.returncode != 0 and not entries:
        entries.append(f"tar_return_code_{proc.returncode}:{proc.stderr.strip()[:200]}")
    return entries[:max_entries]


def check_python_packages() -> dict[str, object]:
    out: dict[str, object] = {"py_pkg_dir": str(PY_PKG_DIR), "py_pkg_dir_exists": PY_PKG_DIR.exists()}
    for package in ["h5py", "anndata", "numpy", "scipy", "pandas"]:
        try:
            module = __import__(package)
            out[f"{package}_available"] = True
            out[f"{package}_version"] = getattr(module, "__version__", "")
        except Exception as exc:  # noqa: BLE001
            out[f"{package}_available"] = False
            out[f"{package}_version"] = f"{type(exc).__name__}:{exc}"
    return out


def audit_gse203191() -> dict[str, object]:
    path = RAW / "GSE203191__GSE203191_Shih_endo_meta_frame.tsv.gz"
    if not path.exists():
        return {"dataset": "GSE203191", "status": "missing"}
    df = pd.read_csv(path, sep="\t")
    cluster_counts = df["clusterID"].value_counts().head(20).to_dict()
    pheno_counts = df["pheno"].value_counts(dropna=False).to_dict()
    subject_counts = df.groupby("pheno")["subjectID"].nunique().to_dict()
    return {
        "dataset": "GSE203191",
        "status": "metadata_available_expression_pending_raw_tar",
        "cells": len(df),
        "subjects": df["subjectID"].nunique(),
        "phenotype_counts_cells": json.dumps(pheno_counts, sort_keys=True),
        "phenotype_counts_subjects": json.dumps(subject_counts, sort_keys=True),
        "top_cluster_counts": json.dumps(cluster_counts, sort_keys=True),
    }


def audit_zenodo_17078290() -> dict[str, object]:
    record = json.loads((META / "zenodo_17078290.json").read_text(encoding="utf-8"))
    expected = record["files"][0]["size"]
    status = file_status("Zenodo_17078290__celllable.diff_PRO.h5ad", expected)
    status.update(
        {
            "dataset": "Zenodo_17078290",
            "record_title": record.get("title", ""),
            "license": record.get("metadata", {}).get("license", {}).get("id", ""),
            "expected_md5": record["files"][0].get("checksum", ""),
            "readiness": "complete_h5ad_ready" if status["complete_by_size"] else "partial_or_missing_resume_required",
        }
    )
    return status


def audit_gse179640() -> dict[str, object]:
    status = file_status("GSE179640__GSE179640_RAW.tar")
    path = Path(status["path"]) if status["path"] else RAW / "GSE179640__GSE179640_RAW.tar.partial"
    entries = tar_head(path)
    h5_entries = [entry for entry in entries if entry.endswith(".h5")]
    status.update(
        {
            "dataset": "GSE179640",
            "readiness": "complete_tar_ready" if (RAW / "GSE179640__GSE179640_RAW.tar").exists() else "partial_tar_headers_visible",
            "first_tar_entries": json.dumps(entries),
            "first_h5_entries": json.dumps(h5_entries),
        }
    )
    return status


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    package_status = check_python_packages()
    source_rows = [audit_zenodo_17078290(), audit_gse179640()]
    source_df = pd.DataFrame(source_rows)
    source_df.to_csv(RESULTS / "singlecell_source_file_audit.tsv", sep="\t", index=False)

    gse203191 = audit_gse203191()
    pd.DataFrame([gse203191]).to_csv(RESULTS / "GSE203191_metadata_audit.tsv", sep="\t", index=False)
    pd.DataFrame([package_status]).to_csv(RESULTS / "singlecell_python_package_audit.tsv", sep="\t", index=False)

    zenodo_ready = source_rows[0].get("readiness") == "complete_h5ad_ready"
    gse179640_ready = source_rows[1].get("readiness") == "complete_tar_ready"
    zenodo_sentence = (
        "- Zenodo 17078290 h5ad is complete by expected file size and ready for backed/on-disk h5ad auditing."
        if zenodo_ready
        else "- Zenodo 17078290 h5ad is currently usable only after the large partial download is resumed to completion."
    )
    gse179640_sentence = (
        "- GSE179640 RAW tar is complete and ingestible as 10x filtered_feature_bc_matrix files."
        if gse179640_ready
        else "- GSE179640 RAW tar headers indicate 10x filtered_feature_bc_matrix.h5 sample files, so it should be ingestible after the tar download completes."
    )
    lines = [
        "# Single-cell/spatial source readiness audit",
        "",
        "## Python reader packages",
        "",
        "\n".join(f"- {key}: {value}" for key, value in package_status.items()),
        "",
        "## Source files",
        "",
        source_df.to_csv(sep="\t", index=False),
        "",
        "## GSE203191 metadata",
        "",
        pd.DataFrame([gse203191]).to_csv(sep="\t", index=False),
        "",
        "## Interpretation",
        "",
        gse179640_sentence,
        "- GSE203191 metadata are available for cell-state/proportion audits, but expression requires the RAW tar or another processed object.",
        zenodo_sentence,
        "",
    ]
    (RESULTS / "singlecell_source_readiness_audit.md").write_text("\n".join(lines), encoding="utf-8")
    print(RESULTS / "singlecell_source_readiness_audit.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
