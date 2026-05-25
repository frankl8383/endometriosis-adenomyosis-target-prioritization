#!/usr/bin/env python3
"""Summarize GSE179640 candidate expression after adaptive tissue QC."""

from __future__ import annotations

import csv
import gzip
import shutil
import sys
import tarfile
import tempfile
from collections import defaultdict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW = PROJECT_ROOT / "data" / "raw_downloads"
RESULTS = PROJECT_ROOT / "results" / "singlecell"
PY_PKG_DIR = PROJECT_ROOT / "tools" / "py_packages"

if sys.version_info[:2] != (3, 12):
    raise SystemExit(
        "Use the bundled Python 3.12 runtime: "
        "python3 "
        "scripts/summarize_gse179640_candidate_expression_qc_filtered.py"
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
    bulk_rows = load_tsv(PROJECT_ROOT / "results" / "bulk" / "bulk_candidate_expression_support_scores.tsv")
    bulk_by_id = {row["gene_id"].split(".")[0]: row for row in bulk_rows}
    out: list[dict[str, str]] = []
    for row in gwas_rows:
        gene_id = row["gene_id"].split(".")[0]
        bulk = bulk_by_id.get(gene_id, {})
        out.append(
            {
                "gene_id": gene_id,
                "gene_symbol": row.get("gene_symbol", ""),
                "genetic_priority": row.get("genetic_priority", ""),
                "ld_neighborhood_class": row.get("ld_neighborhood_class", ""),
                "module_hint_preliminary": row.get("module_hint_preliminary", ""),
                "bulk_expression_support_score_20": bulk.get("bulk_expression_support_score_20", ""),
                "bulk_support_class": bulk.get("bulk_support_class", ""),
            }
        )
    return out


def load_pass_barcodes() -> dict[str, set[str]]:
    path = RESULTS / "GSE179640_cell_qc_metrics.tsv.gz"
    out: dict[str, set[str]] = defaultdict(set)
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["pass_adaptive_qc"] == "1":
                out[row["geo_accession"]].add(row["barcode"])
    return out


def find_feature_index(
    candidate: dict[str, str],
    feature_id_to_index: dict[str, int],
    feature_name_to_indices: dict[str, list[int]],
) -> tuple[int | None, str]:
    gene_id = candidate["gene_id"]
    if gene_id in feature_id_to_index:
        return feature_id_to_index[gene_id], "ensembl_id"
    symbol = candidate["gene_symbol"]
    indices = feature_name_to_indices.get(symbol, [])
    if len(indices) == 1:
        return indices[0], "gene_symbol"
    if len(indices) > 1:
        return indices[0], "gene_symbol_first_of_multiple"
    return None, "missing"


def summarize_h5(path: str, candidates: list[dict[str, str]], pass_barcodes: set[str]) -> list[dict[str, object]]:
    with h5py.File(path, "r") as h5:
        matrix = h5["matrix"]
        shape = tuple(int(x) for x in matrix["shape"][()])
        feature_ids = [item.split(".")[0] for item in decode_array(matrix["features"]["id"][:])]
        feature_names = decode_array(matrix["features"]["name"][:])
        barcodes = decode_array(matrix["barcodes"][:])
        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]
        mat = sparse.csc_matrix((data, indices, indptr), shape=shape)

    keep_idx = np.array([idx for idx, barcode in enumerate(barcodes) if barcode in pass_barcodes], dtype=int)
    mat_keep = mat[:, keep_idx] if keep_idx.size else mat[:, []]
    feature_id_to_index = {gene_id: idx for idx, gene_id in enumerate(feature_ids)}
    feature_name_to_indices: dict[str, list[int]] = defaultdict(list)
    for idx, name in enumerate(feature_names):
        feature_name_to_indices[name].append(idx)

    rows: list[dict[str, object]] = []
    n_cells = int(keep_idx.size)
    for candidate in candidates:
        feature_index, match_basis = find_feature_index(candidate, feature_id_to_index, feature_name_to_indices)
        total_counts = 0
        expressing_cells = 0
        if feature_index is not None and n_cells:
            gene_vec = mat_keep.getrow(feature_index)
            total_counts = int(np.asarray(gene_vec.sum()).ravel()[0])
            expressing_cells = int(gene_vec.getnnz())
        rows.append(
            {
                **candidate,
                "match_basis": match_basis,
                "n_qc_cells": n_cells,
                "candidate_total_counts_qc": total_counts,
                "candidate_expressing_cells_qc": expressing_cells,
                "candidate_prevalence_qc": expressing_cells / n_cells if n_cells else 0,
                "candidate_mean_counts_all_qc_cells": total_counts / n_cells if n_cells else 0,
            }
        )
    return rows


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    candidates = load_candidates()
    pass_barcodes = load_pass_barcodes()
    manifest_rows = load_tsv(RESULTS / "GSE179640_tar_manifest.tsv")
    manifest_by_entry = {row["tar_entry"]: row for row in manifest_rows}
    tar_path = RAW / "GSE179640__GSE179640_RAW.tar"

    expression_rows: list[dict[str, object]] = []
    with tarfile.open(tar_path, "r") as tar:
        h5_members = [
            member
            for member in tar.getmembers()
            if member.isfile()
            and member.name.endswith("_filtered_feature_bc_matrix.h5")
            and manifest_by_entry.get(member.name, {}).get("tissue") != "Patient-Derived Organoid"
        ]
        for member in h5_members:
            meta = manifest_by_entry[member.name]
            sample_id = meta["geo_accession"]
            extracted = tar.extractfile(member)
            if extracted is None:
                raise RuntimeError(f"Could not extract {member.name}")
            with extracted, tempfile.NamedTemporaryFile(suffix=".h5") as tmp:
                shutil.copyfileobj(extracted, tmp)
                tmp.flush()
                rows = summarize_h5(tmp.name, candidates, pass_barcodes.get(sample_id, set()))
            for row in rows:
                row.update(
                    {
                        "tar_entry": member.name,
                        "geo_accession": sample_id,
                        "subject_code": meta.get("source_name", ""),
                        "condition": meta.get("condition", ""),
                        "sample_location": meta.get("sample_location", ""),
                        "tissue": meta.get("tissue", ""),
                    }
                )
                expression_rows.append(row)

    write_tsv(RESULTS / "GSE179640_candidate_expression_by_sample_adaptive_qc.tsv", expression_rows)

    grouped: dict[tuple[str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in expression_rows:
        grouped[
            (
                str(row["gene_id"]),
                str(row["gene_symbol"]),
                str(row["sample_location"]),
                str(row["tissue"]),
            )
        ].append(row)

    summary_rows: list[dict[str, object]] = []
    for (gene_id, gene_symbol, sample_location, tissue), rows in grouped.items():
        prevalences = np.array([float(row["candidate_prevalence_qc"]) for row in rows], dtype=float)
        counts = np.array([float(row["candidate_total_counts_qc"]) for row in rows], dtype=float)
        expressing = np.array([float(row["candidate_expressing_cells_qc"]) for row in rows], dtype=float)
        qc_cells = np.array([float(row["n_qc_cells"]) for row in rows], dtype=float)
        detected_samples = sum(1 for row in rows if int(row["candidate_expressing_cells_qc"]) > 0)
        template = rows[0]
        summary_rows.append(
            {
                "gene_id": gene_id,
                "gene_symbol": gene_symbol,
                "sample_location": sample_location,
                "tissue": tissue,
                "n_samples": len(rows),
                "median_qc_cells": float(np.median(qc_cells)) if len(qc_cells) else 0,
                "n_detected_samples": detected_samples,
                "sample_detection_fraction_qc": detected_samples / len(rows) if rows else 0,
                "median_prevalence_qc": float(np.median(prevalences)) if len(prevalences) else 0,
                "mean_prevalence_qc": float(np.mean(prevalences)) if len(prevalences) else 0,
                "median_total_counts_qc": float(np.median(counts)) if len(counts) else 0,
                "median_expressing_cells_qc": float(np.median(expressing)) if len(expressing) else 0,
                "genetic_priority": template.get("genetic_priority", ""),
                "ld_neighborhood_class": template.get("ld_neighborhood_class", ""),
                "module_hint_preliminary": template.get("module_hint_preliminary", ""),
                "bulk_expression_support_score_20": template.get("bulk_expression_support_score_20", ""),
                "bulk_support_class": template.get("bulk_support_class", ""),
            }
        )
    summary_rows.sort(key=lambda row: (str(row["sample_location"]), -float(row["median_prevalence_qc"]), str(row["gene_symbol"])))
    write_tsv(RESULTS / "GSE179640_candidate_expression_summary_by_location_adaptive_qc.tsv", summary_rows)

    high_bulk = [
        row
        for row in summary_rows
        if row["bulk_support_class"] == "high_bulk_support"
        and row["sample_location"] in {"Ectopic", "Ectopic Adjacent", "Ectopic Ovary"}
    ]
    high_bulk_top = sorted(high_bulk, key=lambda row: -float(row["median_prevalence_qc"]))[:15]
    lines = [
        "# GSE179640 adaptive-QC candidate expression summary",
        "",
        f"- Candidate genes summarized: {len(candidates)}",
        f"- Expression rows after adaptive QC: {len(expression_rows)}",
        "- This summary uses only cells passing `pass_adaptive_qc == 1` in `GSE179640_cell_qc_metrics.tsv.gz`.",
        "",
        "## High-bulk-support candidates with highest QC-filtered lesion/adjacent/ovary detectability",
        "",
        "gene_symbol\tgene_id\tsample_location\ttissue\tn_samples\tmedian_qc_cells\tmedian_prevalence_qc\tbulk_score",
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
                    f"{float(row['median_qc_cells']):.0f}",
                    f"{float(row['median_prevalence_qc']):.4f}",
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
            "- QC-filtered detectability is more defensible than raw barcode detectability, but it is still not cell-type localization.",
            "- The result supports prioritizing measurable candidate genes for downstream annotated-cell-state analysis.",
            "- Genes with low prevalence after QC should be treated carefully in single-cell localization and may rely more on bulk/spatial or adenomyosis evidence.",
            "",
        ]
    )
    (RESULTS / "GSE179640_candidate_expression_adaptive_qc_summary.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )
    print(RESULTS / "GSE179640_candidate_expression_adaptive_qc_summary.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
