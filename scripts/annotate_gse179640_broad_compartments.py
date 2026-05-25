#!/usr/bin/env python3
"""Conservative marker-panel broad-compartment annotation for GSE179640."""

from __future__ import annotations

import csv
import gzip
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
        "scripts/annotate_gse179640_broad_compartments.py"
    )

sys.path.insert(0, str(PY_PKG_DIR))

import h5py  # noqa: E402
import numpy as np  # noqa: E402
from scipy import sparse  # noqa: E402


MARKER_PANELS: dict[str, list[str]] = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "MUC1", "MSLN"],
    "stromal_fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "PDGFRA"],
    "endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "KDR", "EMCN", "RAMP2"],
    "mural_smooth_muscle": ["RGS5", "ACTA2", "TAGLN", "MYH11", "MCAM", "PDGFRB", "CSPG4"],
    "t_nk": ["CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "PRF1", "KLRD1"],
    "myeloid_macrophage": ["LST1", "C1QA", "C1QB", "C1QC", "CD68", "MS4A7", "CSF1R"],
    "b_plasma": ["MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "IGKC"],
    "mast": ["TPSAB1", "TPSB2", "CPA3", "KIT"],
    "cycling": ["MKI67", "TOP2A", "STMN1", "PCNA", "UBE2C"],
}


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


def load_pass_barcodes() -> dict[str, set[str]]:
    path = RESULTS / "GSE179640_cell_qc_metrics.tsv.gz"
    out: dict[str, set[str]] = defaultdict(set)
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["pass_adaptive_qc"] == "1":
                out[row["geo_accession"]].add(row["barcode"])
    return out


def load_candidates() -> list[dict[str, str]]:
    gwas_rows = load_tsv(PROJECT_ROOT / "results" / "gwas" / "gwas_candidate_gene_universe.tsv")
    bulk_rows = load_tsv(PROJECT_ROOT / "results" / "bulk" / "bulk_candidate_expression_support_scores.tsv")
    bulk_by_id = {row["gene_id"].split(".")[0]: row for row in bulk_rows}
    candidates: list[dict[str, str]] = []
    for row in gwas_rows:
        gene_id = row["gene_id"].split(".")[0]
        bulk = bulk_by_id.get(gene_id, {})
        candidates.append(
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
    return candidates


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


def load_10x_h5(path: str, pass_barcodes: set[str]) -> dict[str, object]:
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
    kept_barcodes = [barcodes[idx] for idx in keep_idx]
    return {
        "mat": mat[:, keep_idx] if keep_idx.size else mat[:, []],
        "feature_ids": feature_ids,
        "feature_names": feature_names,
        "barcodes": kept_barcodes,
    }


def score_panels(mat: sparse.csc_matrix, feature_names: list[str]) -> tuple[np.ndarray, list[str], list[dict[str, object]]]:
    feature_name_to_indices: dict[str, list[int]] = defaultdict(list)
    for idx, name in enumerate(feature_names):
        feature_name_to_indices[name].append(idx)
    n_cells = mat.shape[1]
    library_size = np.asarray(mat.sum(axis=0)).ravel().astype(float)
    scale = np.divide(10000.0, library_size, out=np.zeros_like(library_size), where=library_size > 0)
    panel_names = list(MARKER_PANELS)
    scores = np.zeros((n_cells, len(panel_names)), dtype=float)
    coverage_rows: list[dict[str, object]] = []
    for panel_idx, panel in enumerate(panel_names):
        markers = MARKER_PANELS[panel]
        indices = [feature_name_to_indices[marker][0] for marker in markers if marker in feature_name_to_indices]
        coverage_rows.append(
            {
                "panel": panel,
                "markers_requested": "|".join(markers),
                "n_markers_requested": len(markers),
                "markers_present": "|".join(marker for marker in markers if marker in feature_name_to_indices),
                "n_markers_present": len(indices),
            }
        )
        if not indices or n_cells == 0:
            continue
        sub = mat[indices, :].toarray().astype(float)
        normalized = np.log1p(sub * scale.reshape(1, -1))
        scores[:, panel_idx] = normalized.mean(axis=0)
    return scores, panel_names, coverage_rows


def assign_labels(scores: np.ndarray, panel_names: list[str]) -> tuple[list[str], np.ndarray, np.ndarray]:
    if scores.shape[0] == 0:
        return [], np.array([]), np.array([])
    order = np.argsort(-scores, axis=1)
    top_idx = order[:, 0]
    second_idx = order[:, 1] if scores.shape[1] > 1 else top_idx
    top_score = scores[np.arange(scores.shape[0]), top_idx]
    second_score = scores[np.arange(scores.shape[0]), second_idx]
    labels = []
    for idx, score, margin in zip(top_idx, top_score, top_score - second_score):
        if score < 0.05:
            labels.append("low_marker_unassigned")
        elif margin < 0.01:
            labels.append(f"ambiguous_{panel_names[int(idx)]}")
        else:
            labels.append(panel_names[int(idx)])
    return labels, top_score, second_score


def candidate_expression_by_label(
    mat: sparse.csc_matrix,
    labels: list[str],
    feature_ids: list[str],
    feature_names: list[str],
    candidates: list[dict[str, str]],
) -> list[dict[str, object]]:
    feature_id_to_index = {gene_id: idx for idx, gene_id in enumerate(feature_ids)}
    feature_name_to_indices: dict[str, list[int]] = defaultdict(list)
    for idx, name in enumerate(feature_names):
        feature_name_to_indices[name].append(idx)
    label_array = np.array(labels, dtype=object)
    rows: list[dict[str, object]] = []
    for label in sorted(set(labels)):
        mask = label_array == label
        n_cells = int(np.sum(mask))
        cell_indices = np.where(mask)[0]
        mat_label = mat[:, cell_indices] if n_cells else mat[:, []]
        for candidate in candidates:
            feature_index, match_basis = find_feature_index(candidate, feature_id_to_index, feature_name_to_indices)
            total_counts = 0
            expressing_cells = 0
            if feature_index is not None and n_cells:
                gene_vec = mat_label.getrow(feature_index)
                total_counts = int(np.asarray(gene_vec.sum()).ravel()[0])
                expressing_cells = int(gene_vec.getnnz())
            rows.append(
                {
                    **candidate,
                    "broad_compartment": label,
                    "n_cells": n_cells,
                    "candidate_total_counts": total_counts,
                    "candidate_expressing_cells": expressing_cells,
                    "candidate_prevalence": expressing_cells / n_cells if n_cells else 0,
                    "match_basis": match_basis,
                }
            )
    return rows


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    tar_path = RAW / "GSE179640__GSE179640_RAW.tar"
    manifest_rows = load_tsv(RESULTS / "GSE179640_tar_manifest.tsv")
    manifest_by_entry = {row["tar_entry"]: row for row in manifest_rows}
    pass_barcodes = load_pass_barcodes()
    candidates = load_candidates()

    cell_annotation_path = RESULTS / "GSE179640_broad_compartment_cell_annotations.tsv.gz"
    cell_summary_rows: list[dict[str, object]] = []
    candidate_sample_rows: list[dict[str, object]] = []
    coverage_rows_all: list[dict[str, object]] = []

    with gzip.open(cell_annotation_path, "wt", encoding="utf-8", newline="") as cell_handle:
        cell_writer = csv.DictWriter(
            cell_handle,
            fieldnames=[
                "cell_id",
                "barcode",
                "geo_accession",
                "subject_code",
                "sample_location",
                "tissue",
                "broad_compartment",
                "top_marker_score",
                "second_marker_score",
            ],
            delimiter="\t",
        )
        cell_writer.writeheader()

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
                    loaded = load_10x_h5(tmp.name, pass_barcodes.get(sample_id, set()))

                mat = loaded["mat"]
                barcodes = loaded["barcodes"]
                feature_names = loaded["feature_names"]
                feature_ids = loaded["feature_ids"]
                scores, panel_names, coverage_rows = score_panels(mat, feature_names)
                labels, top_score, second_score = assign_labels(scores, panel_names)
                counts = Counter(labels)
                for row in coverage_rows:
                    row.update({"geo_accession": sample_id, "tar_entry": member.name})
                    coverage_rows_all.append(row)

                for idx, barcode in enumerate(barcodes):
                    cell_writer.writerow(
                        {
                            "cell_id": f"{sample_id}_{barcode}",
                            "barcode": barcode,
                            "geo_accession": sample_id,
                            "subject_code": meta.get("source_name", ""),
                            "sample_location": meta.get("sample_location", ""),
                            "tissue": meta.get("tissue", ""),
                            "broad_compartment": labels[idx],
                            "top_marker_score": f"{float(top_score[idx]):.5f}",
                            "second_marker_score": f"{float(second_score[idx]):.5f}",
                        }
                    )

                for label, n_cells in sorted(counts.items()):
                    cell_summary_rows.append(
                        {
                            "geo_accession": sample_id,
                            "subject_code": meta.get("source_name", ""),
                            "sample_location": meta.get("sample_location", ""),
                            "tissue": meta.get("tissue", ""),
                            "broad_compartment": label,
                            "n_cells": n_cells,
                            "fraction": n_cells / len(labels) if labels else 0,
                        }
                    )

                candidate_rows = candidate_expression_by_label(mat, labels, feature_ids, feature_names, candidates)
                for row in candidate_rows:
                    row.update(
                        {
                            "geo_accession": sample_id,
                            "subject_code": meta.get("source_name", ""),
                            "sample_location": meta.get("sample_location", ""),
                            "tissue": meta.get("tissue", ""),
                        }
                    )
                    candidate_sample_rows.append(row)

    write_tsv(RESULTS / "GSE179640_broad_compartment_marker_coverage.tsv", coverage_rows_all)
    write_tsv(RESULTS / "GSE179640_broad_compartment_summary_by_sample.tsv", cell_summary_rows)
    write_tsv(RESULTS / "GSE179640_candidate_expression_by_broad_compartment_sample.tsv", candidate_sample_rows)

    location_compartment_counts: dict[tuple[str, str, str], Counter[str]] = defaultdict(Counter)
    for row in cell_summary_rows:
        key = (str(row["sample_location"]), str(row["tissue"]), str(row["broad_compartment"]))
        location_compartment_counts[key]["n_samples_with_compartment"] += 1
        location_compartment_counts[key]["n_cells"] += int(row["n_cells"])
    total_by_location: Counter[tuple[str, str]] = Counter()
    for (sample_location, tissue, _compartment), counts in location_compartment_counts.items():
        total_by_location[(sample_location, tissue)] += int(counts["n_cells"])
    location_compartment_rows = []
    for (sample_location, tissue, compartment), counts in sorted(location_compartment_counts.items()):
        total = total_by_location[(sample_location, tissue)]
        location_compartment_rows.append(
            {
                "sample_location": sample_location,
                "tissue": tissue,
                "broad_compartment": compartment,
                "n_samples_with_compartment": counts["n_samples_with_compartment"],
                "n_cells": counts["n_cells"],
                "fraction_within_location": counts["n_cells"] / total if total else 0,
            }
        )
    write_tsv(RESULTS / "GSE179640_broad_compartment_summary_by_location.tsv", location_compartment_rows)

    grouped: dict[tuple[str, str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in candidate_sample_rows:
        grouped[
            (
                str(row["gene_id"]),
                str(row["gene_symbol"]),
                str(row["sample_location"]),
                str(row["tissue"]),
                str(row["broad_compartment"]),
            )
        ].append(row)

    summary_rows: list[dict[str, object]] = []
    for (gene_id, gene_symbol, sample_location, tissue, broad_compartment), rows in grouped.items():
        prevalences = np.array([float(row["candidate_prevalence"]) for row in rows], dtype=float)
        n_cells = np.array([float(row["n_cells"]) for row in rows], dtype=float)
        template = rows[0]
        summary_rows.append(
            {
                "gene_id": gene_id,
                "gene_symbol": gene_symbol,
                "sample_location": sample_location,
                "tissue": tissue,
                "broad_compartment": broad_compartment,
                "n_samples_with_compartment": len(rows),
                "median_compartment_cells": float(np.median(n_cells)) if len(n_cells) else 0,
                "median_prevalence": float(np.median(prevalences)) if len(prevalences) else 0,
                "mean_prevalence": float(np.mean(prevalences)) if len(prevalences) else 0,
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
            str(row["broad_compartment"]),
            -float(row["median_prevalence"]),
            str(row["gene_symbol"]),
        )
    )
    write_tsv(RESULTS / "GSE179640_candidate_expression_summary_by_broad_compartment.tsv", summary_rows)

    high_bulk = [
        row
        for row in summary_rows
        if row["bulk_support_class"] == "high_bulk_support"
        and row["sample_location"] in {"Ectopic", "Ectopic Adjacent", "Ectopic Ovary"}
        and float(row["median_compartment_cells"]) >= 20
    ]
    top_high_bulk = sorted(high_bulk, key=lambda row: -float(row["median_prevalence"]))[:25]
    compartment_totals = Counter()
    for row in cell_summary_rows:
        compartment_totals[str(row["broad_compartment"])] += int(row["n_cells"])

    lines = [
        "# GSE179640 broad-compartment marker annotation",
        "",
        f"- QC-passing cells annotated: {sum(compartment_totals.values())}",
        f"- Broad-compartment totals: `{json.dumps(compartment_totals, sort_keys=True)}`",
        f"- Cell annotation table: `{cell_annotation_path}`",
        "",
        "## Top high-bulk-support candidates by QC-filtered broad compartment",
        "",
        "gene_symbol\tgene_id\tsample_location\tbroad_compartment\tn_samples\tmedian_cells\tmedian_prevalence\tbulk_score",
    ]
    for row in top_high_bulk:
        lines.append(
            "\t".join(
                [
                    str(row["gene_symbol"]),
                    str(row["gene_id"]),
                    str(row["sample_location"]),
                    str(row["broad_compartment"]),
                    str(row["n_samples_with_compartment"]),
                    f"{float(row['median_compartment_cells']):.0f}",
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
            "- This is a conservative marker-panel broad-compartment annotation, not a final cluster-level atlas.",
            "- Ambiguous and low-marker cells are kept visible instead of being forced into confident labels.",
            "- Results can support candidate triage and figure planning, but final cell-state localization still requires clustering/author-label comparison or independent validation.",
            "- Candidate summaries are sample-aware by computing gene prevalence per sample-compartment before across-sample medians.",
            "",
        ]
    )
    (RESULTS / "GSE179640_broad_compartment_annotation_summary.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )
    print(RESULTS / "GSE179640_broad_compartment_annotation_summary.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
