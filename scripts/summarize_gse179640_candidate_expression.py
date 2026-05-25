#!/usr/bin/env python3
"""Summarize GWAS candidate expression/prevalence in GSE179640 10x h5 files."""

from __future__ import annotations

import csv
import json
import shutil
import sys
import tarfile
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW = PROJECT_ROOT / "data" / "raw_downloads"
RESULTS = PROJECT_ROOT / "results" / "singlecell"
PY_PKG_DIR = PROJECT_ROOT / "tools" / "py_packages"

if sys.version_info[:2] != (3, 12):
    raise SystemExit(
        "Use the bundled Python 3.12 runtime: "
        "python3 "
        "scripts/summarize_gse179640_candidate_expression.py"
    )

sys.path.insert(0, str(PY_PKG_DIR))

import h5py  # noqa: E402
import numpy as np  # noqa: E402
from scipy import sparse  # noqa: E402


def decode_array(values) -> list[str]:
    out: list[str] = []
    for value in values:
        if isinstance(value, bytes):
            out.append(value.decode("utf-8", errors="replace"))
        else:
            out.append(str(value))
    return out


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_candidates() -> list[dict[str, str]]:
    gwas_rows = load_tsv(PROJECT_ROOT / "results" / "gwas" / "gwas_candidate_gene_universe.tsv")
    bulk_path = PROJECT_ROOT / "results" / "bulk" / "bulk_candidate_expression_support_scores.tsv"
    bulk_by_id: dict[str, dict[str, str]] = {}
    if bulk_path.exists():
        bulk_by_id = {row["gene_id"]: row for row in load_tsv(bulk_path)}
    candidates: list[dict[str, str]] = []
    for row in gwas_rows:
        bulk = bulk_by_id.get(row["gene_id"], {})
        candidates.append(
            {
                "gene_id": row["gene_id"].split(".")[0],
                "gene_symbol": row.get("gene_symbol", ""),
                "genetic_priority": row.get("genetic_priority", ""),
                "ld_neighborhood_class": row.get("ld_neighborhood_class", ""),
                "module_hint_preliminary": row.get("module_hint_preliminary", ""),
                "bulk_expression_support_score_20": bulk.get("bulk_expression_support_score_20", ""),
                "bulk_support_class": bulk.get("bulk_support_class", ""),
            }
        )
    return candidates


def find_feature_index(
    candidate: dict[str, str],
    feature_id_to_index: dict[str, int],
    feature_name_to_indices: dict[str, list[int]],
) -> tuple[int | None, str]:
    gene_id = candidate["gene_id"].split(".")[0]
    if gene_id in feature_id_to_index:
        return feature_id_to_index[gene_id], "ensembl_id"
    symbol = candidate["gene_symbol"]
    indices = feature_name_to_indices.get(symbol, [])
    if len(indices) == 1:
        return indices[0], "gene_symbol"
    if len(indices) > 1:
        return indices[0], "gene_symbol_first_of_multiple"
    return None, "missing"


def h5_candidate_rows(path: str, candidates: list[dict[str, str]]) -> tuple[list[dict[str, object]], dict[str, object]]:
    with h5py.File(path, "r") as h5:
        matrix = h5["matrix"]
        shape = tuple(int(x) for x in matrix["shape"][()])
        feature_ids = [item.split(".")[0] for item in decode_array(matrix["features"]["id"][:])]
        feature_names = decode_array(matrix["features"]["name"][:])
        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]
        mat = sparse.csc_matrix((data, indices, indptr), shape=shape)

    feature_id_to_index = {gene_id: idx for idx, gene_id in enumerate(feature_ids)}
    feature_name_to_indices: dict[str, list[int]] = defaultdict(list)
    for idx, name in enumerate(feature_names):
        feature_name_to_indices[name].append(idx)

    n_features, n_cells = shape
    total_umi = int(np.asarray(mat.sum()).ravel()[0])
    rows: list[dict[str, object]] = []
    match_basis_counts: Counter[str] = Counter()
    for candidate in candidates:
        feature_index, match_basis = find_feature_index(candidate, feature_id_to_index, feature_name_to_indices)
        match_basis_counts[match_basis] += 1
        total_counts = 0
        expressing_cells = 0
        if feature_index is not None:
            gene_vec = mat.getrow(feature_index)
            total_counts = int(np.asarray(gene_vec.sum()).ravel()[0])
            expressing_cells = int(gene_vec.getnnz())
        rows.append(
            {
                **candidate,
                "feature_index": "" if feature_index is None else feature_index,
                "match_basis": match_basis,
                "n_cells": n_cells,
                "candidate_total_counts": total_counts,
                "candidate_expressing_cells": expressing_cells,
                "candidate_prevalence": expressing_cells / n_cells if n_cells else 0,
                "candidate_mean_counts_all_cells": total_counts / n_cells if n_cells else 0,
            }
        )
    sample_metrics = {
        "n_features": n_features,
        "n_cells": n_cells,
        "total_umi": total_umi,
        "match_basis_counts": json.dumps(match_basis_counts, sort_keys=True),
    }
    return rows, sample_metrics


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    tar_path = RAW / "GSE179640__GSE179640_RAW.tar"
    manifest_rows = load_tsv(RESULTS / "GSE179640_tar_manifest.tsv")
    manifest_by_entry = {row["tar_entry"]: row for row in manifest_rows}
    candidates = load_candidates()

    expression_rows: list[dict[str, object]] = []
    sample_rows: list[dict[str, object]] = []
    with tarfile.open(tar_path, "r") as tar:
        h5_members = [
            member
            for member in tar.getmembers()
            if member.isfile()
            and member.name.endswith("_filtered_feature_bc_matrix.h5")
            and manifest_by_entry.get(member.name, {}).get("tissue") != "Patient-Derived Organoid"
        ]
        for member in h5_members:
            extracted = tar.extractfile(member)
            if extracted is None:
                raise RuntimeError(f"Could not extract {member.name}")
            with extracted, tempfile.NamedTemporaryFile(suffix=".h5") as tmp:
                shutil.copyfileobj(extracted, tmp)
                tmp.flush()
                rows, metrics = h5_candidate_rows(tmp.name, candidates)

            meta = manifest_by_entry[member.name]
            sample_summary = {
                "tar_entry": member.name,
                "geo_accession": meta.get("geo_accession", ""),
                "subject_code": meta.get("source_name", ""),
                "condition": meta.get("condition", ""),
                "sample_location": meta.get("sample_location", ""),
                "tissue": meta.get("tissue", ""),
                "title": meta.get("title", ""),
                **metrics,
            }
            sample_rows.append(sample_summary)
            for row in rows:
                row.update(sample_summary)
                expression_rows.append(row)

    write_tsv(RESULTS / "GSE179640_candidate_expression_by_sample.tsv", expression_rows)
    write_tsv(RESULTS / "GSE179640_sample_qc_expression_summary.tsv", sample_rows)

    grouped: dict[tuple[str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in expression_rows:
        key = (
            str(row["gene_id"]),
            str(row["gene_symbol"]),
            str(row["sample_location"]),
            str(row["tissue"]),
        )
        grouped[key].append(row)

    summary_rows: list[dict[str, object]] = []
    for (gene_id, gene_symbol, sample_location, tissue), rows in grouped.items():
        prevalences = np.array([float(row["candidate_prevalence"]) for row in rows], dtype=float)
        counts = np.array([float(row["candidate_total_counts"]) for row in rows], dtype=float)
        expressing = np.array([float(row["candidate_expressing_cells"]) for row in rows], dtype=float)
        detected_samples = sum(1 for row in rows if int(row["candidate_expressing_cells"]) > 0)
        template = rows[0]
        summary_rows.append(
            {
                "gene_id": gene_id,
                "gene_symbol": gene_symbol,
                "sample_location": sample_location,
                "tissue": tissue,
                "n_samples": len(rows),
                "n_detected_samples": detected_samples,
                "sample_detection_fraction": detected_samples / len(rows) if rows else 0,
                "median_prevalence": float(np.median(prevalences)) if len(prevalences) else 0,
                "mean_prevalence": float(np.mean(prevalences)) if len(prevalences) else 0,
                "median_total_counts": float(np.median(counts)) if len(counts) else 0,
                "median_expressing_cells": float(np.median(expressing)) if len(expressing) else 0,
                "genetic_priority": template.get("genetic_priority", ""),
                "ld_neighborhood_class": template.get("ld_neighborhood_class", ""),
                "module_hint_preliminary": template.get("module_hint_preliminary", ""),
                "bulk_expression_support_score_20": template.get("bulk_expression_support_score_20", ""),
                "bulk_support_class": template.get("bulk_support_class", ""),
            }
        )
    summary_rows.sort(
        key=lambda row: (
            str(row["sample_location"]),
            -float(row["median_prevalence"]),
            str(row["gene_symbol"]),
        )
    )
    write_tsv(RESULTS / "GSE179640_candidate_expression_summary_by_location.tsv", summary_rows)

    location_counts = Counter(row["sample_location"] for row in sample_rows)
    subject_by_location: dict[str, set[str]] = defaultdict(set)
    for row in sample_rows:
        subject_by_location[str(row["sample_location"])].add(str(row["subject_code"]))
    high_bulk = [
        row
        for row in summary_rows
        if row["bulk_support_class"] == "high_bulk_support"
        and row["sample_location"] in {"Ectopic", "Ectopic Adjacent", "Ectopic Ovary"}
    ]
    high_bulk_top = sorted(high_bulk, key=lambda row: -float(row["median_prevalence"]))[:15]

    lines = [
        "# GSE179640 candidate expression pre-QC summary",
        "",
        f"- Primary tissue 10x h5 samples summarized: {len(sample_rows)}",
        f"- Candidate genes summarized per sample: {len(candidates)}",
        f"- Expression rows: {len(expression_rows)}",
        f"- Sample-location counts: `{json.dumps(location_counts, sort_keys=True)}`",
        f"- Subject counts by sample location: `{json.dumps({k: len(v) for k, v in subject_by_location.items()}, sort_keys=True)}`",
        "",
        "## High-bulk-support candidates with highest lesion/adjacent/ovary single-cell detectability",
        "",
        "gene_symbol\tgene_id\tsample_location\ttissue\tn_samples\tmedian_prevalence\tbulk_score",
    ]
    for row in high_bulk_top:
        lines.append(
            "\t".join(
                [
                    str(row["gene_symbol"]),
                    str(row["gene_id"]),
                    str(row["sample_location"]),
                    str(row["tissue"]),
                    str(row["n_samples"]),
                    f"{float(row['median_prevalence']):.4f}",
                    str(row["bulk_expression_support_score_20"]),
                ]
            )
        )
    lines.extend(
        [
            "",
            "## Self-review",
            "",
            "Verdict: **PASS_WITH_CONDITIONS**",
            "",
            "- This is a raw 10x h5 detectability/prevalence screen, not a cell-type-specific localization result.",
            "- Organoid samples are excluded from the primary summary to avoid mixing in vitro and tissue-derived states.",
            "- The output is suitable for deciding which candidates are measurable enough for downstream Seurat/annotation/pseudobulk analysis.",
            "- Biological conclusions require downstream cell QC, doublet handling, cell-type annotation, and donor-aware modeling.",
            "",
        ]
    )
    (RESULTS / "GSE179640_candidate_expression_pre_qc_summary.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )
    print(RESULTS / "GSE179640_candidate_expression_pre_qc_summary.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
