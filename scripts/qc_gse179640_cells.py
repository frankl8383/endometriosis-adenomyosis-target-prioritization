#!/usr/bin/env python3
"""Compute cell-level QC metrics for primary-tissue GSE179640 10x h5 files."""

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
        "scripts/qc_gse179640_cells.py"
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


def quantile(values: np.ndarray, q: float) -> float:
    if values.size == 0:
        return 0.0
    return float(np.quantile(values, q))


def summarize(values: np.ndarray, prefix: str) -> dict[str, object]:
    return {
        f"{prefix}_min": float(np.min(values)) if values.size else 0,
        f"{prefix}_q01": quantile(values, 0.01),
        f"{prefix}_q05": quantile(values, 0.05),
        f"{prefix}_median": quantile(values, 0.50),
        f"{prefix}_q95": quantile(values, 0.95),
        f"{prefix}_q99": quantile(values, 0.99),
        f"{prefix}_q995": quantile(values, 0.995),
        f"{prefix}_max": float(np.max(values)) if values.size else 0,
    }


def qc_from_10x_h5(path: str) -> dict[str, object]:
    with h5py.File(path, "r") as h5:
        matrix = h5["matrix"]
        shape = tuple(int(x) for x in matrix["shape"][()])
        feature_names = decode_array(matrix["features"]["name"][:])
        barcodes = decode_array(matrix["barcodes"][:])
        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]
        mat = sparse.csc_matrix((data, indices, indptr), shape=shape)

    n_counts = np.asarray(mat.sum(axis=0)).ravel().astype(float)
    n_features = np.diff(mat.indptr).astype(float)
    mt_gene_indices = np.array([idx for idx, name in enumerate(feature_names) if name.startswith("MT-")], dtype=int)
    ribo_gene_indices = np.array(
        [idx for idx, name in enumerate(feature_names) if name.startswith("RPS") or name.startswith("RPL")],
        dtype=int,
    )
    if mt_gene_indices.size:
        mt_counts = np.asarray(mat[mt_gene_indices, :].sum(axis=0)).ravel().astype(float)
    else:
        mt_counts = np.zeros_like(n_counts)
    if ribo_gene_indices.size:
        ribo_counts = np.asarray(mat[ribo_gene_indices, :].sum(axis=0)).ravel().astype(float)
    else:
        ribo_counts = np.zeros_like(n_counts)
    pct_mt = np.divide(mt_counts, n_counts, out=np.zeros_like(mt_counts), where=n_counts > 0) * 100
    pct_ribo = np.divide(ribo_counts, n_counts, out=np.zeros_like(ribo_counts), where=n_counts > 0) * 100
    return {
        "barcodes": barcodes,
        "n_counts": n_counts,
        "n_features": n_features,
        "pct_mt": pct_mt,
        "pct_ribo": pct_ribo,
        "n_mt_genes": int(mt_gene_indices.size),
        "n_ribo_genes": int(ribo_gene_indices.size),
    }


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    tar_path = RAW / "GSE179640__GSE179640_RAW.tar"
    manifest_rows = load_tsv(RESULTS / "GSE179640_tar_manifest.tsv")
    manifest_by_entry = {row["tar_entry"]: row for row in manifest_rows}
    metrics_path = RESULTS / "GSE179640_cell_qc_metrics.tsv.gz"

    sample_rows: list[dict[str, object]] = []
    tissue_rows: list[dict[str, object]] = []
    tissue_accumulators: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)

    with gzip.open(metrics_path, "wt", encoding="utf-8", newline="") as metrics_handle:
        fieldnames = [
            "cell_id",
            "barcode",
            "tar_entry",
            "geo_accession",
            "subject_code",
            "condition",
            "sample_location",
            "tissue",
            "n_counts",
            "n_features",
            "pct_mt",
            "pct_ribo",
            "pass_lenient_qc",
            "pass_standard_qc",
            "pass_adaptive_qc",
        ]
        writer = csv.DictWriter(metrics_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

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
                    qc = qc_from_10x_h5(tmp.name)

                n_counts = qc["n_counts"]
                n_features = qc["n_features"]
                pct_mt = qc["pct_mt"]
                pct_ribo = qc["pct_ribo"]
                upper_features_adaptive = max(5000.0, quantile(n_features, 0.995))
                upper_counts_adaptive = max(25000.0, quantile(n_counts, 0.995))
                pass_lenient = (n_features >= 200) & (n_counts >= 300) & (pct_mt <= 25)
                pass_standard = (n_features >= 300) & (n_features <= 7500) & (n_counts >= 500) & (pct_mt <= 20)
                pass_adaptive = (
                    (n_features >= 300)
                    & (n_features <= upper_features_adaptive)
                    & (n_counts >= 500)
                    & (n_counts <= upper_counts_adaptive)
                    & (pct_mt <= 25)
                )
                meta = manifest_by_entry[member.name]
                sample_id = str(meta.get("geo_accession", ""))
                for idx, barcode in enumerate(qc["barcodes"]):
                    row = {
                        "cell_id": f"{sample_id}_{barcode}",
                        "barcode": barcode,
                        "tar_entry": member.name,
                        "geo_accession": sample_id,
                        "subject_code": meta.get("source_name", ""),
                        "condition": meta.get("condition", ""),
                        "sample_location": meta.get("sample_location", ""),
                        "tissue": meta.get("tissue", ""),
                        "n_counts": int(n_counts[idx]),
                        "n_features": int(n_features[idx]),
                        "pct_mt": f"{pct_mt[idx]:.4f}",
                        "pct_ribo": f"{pct_ribo[idx]:.4f}",
                        "pass_lenient_qc": int(bool(pass_lenient[idx])),
                        "pass_standard_qc": int(bool(pass_standard[idx])),
                        "pass_adaptive_qc": int(bool(pass_adaptive[idx])),
                    }
                    writer.writerow(row)

                sample_row: dict[str, object] = {
                    "tar_entry": member.name,
                    "geo_accession": sample_id,
                    "subject_code": meta.get("source_name", ""),
                    "condition": meta.get("condition", ""),
                    "sample_location": meta.get("sample_location", ""),
                    "tissue": meta.get("tissue", ""),
                    "n_cells_raw": int(n_counts.size),
                    "n_cells_lenient_qc": int(np.sum(pass_lenient)),
                    "n_cells_standard_qc": int(np.sum(pass_standard)),
                    "n_cells_adaptive_qc": int(np.sum(pass_adaptive)),
                    "lenient_qc_fraction": float(np.mean(pass_lenient)) if n_counts.size else 0,
                    "standard_qc_fraction": float(np.mean(pass_standard)) if n_counts.size else 0,
                    "adaptive_qc_fraction": float(np.mean(pass_adaptive)) if n_counts.size else 0,
                    "adaptive_upper_features": upper_features_adaptive,
                    "adaptive_upper_counts": upper_counts_adaptive,
                    "n_mt_genes": qc["n_mt_genes"],
                    "n_ribo_genes": qc["n_ribo_genes"],
                }
                sample_row.update(summarize(n_counts, "n_counts"))
                sample_row.update(summarize(n_features, "n_features"))
                sample_row.update(summarize(pct_mt, "pct_mt"))
                sample_row.update(summarize(pct_ribo, "pct_ribo"))
                sample_rows.append(sample_row)

                key = (str(meta.get("sample_location", "")), str(meta.get("tissue", "")))
                tissue_accumulators[key]["n_samples"] += 1
                tissue_accumulators[key]["n_cells_raw"] += int(n_counts.size)
                tissue_accumulators[key]["n_cells_lenient_qc"] += int(np.sum(pass_lenient))
                tissue_accumulators[key]["n_cells_standard_qc"] += int(np.sum(pass_standard))
                tissue_accumulators[key]["n_cells_adaptive_qc"] += int(np.sum(pass_adaptive))

    for (sample_location, tissue), counts in sorted(tissue_accumulators.items()):
        raw = counts["n_cells_raw"]
        tissue_rows.append(
            {
                "sample_location": sample_location,
                "tissue": tissue,
                "n_samples": counts["n_samples"],
                "n_cells_raw": raw,
                "n_cells_lenient_qc": counts["n_cells_lenient_qc"],
                "n_cells_standard_qc": counts["n_cells_standard_qc"],
                "n_cells_adaptive_qc": counts["n_cells_adaptive_qc"],
                "lenient_qc_fraction": counts["n_cells_lenient_qc"] / raw if raw else 0,
                "standard_qc_fraction": counts["n_cells_standard_qc"] / raw if raw else 0,
                "adaptive_qc_fraction": counts["n_cells_adaptive_qc"] / raw if raw else 0,
            }
        )

    write_tsv(RESULTS / "GSE179640_cell_qc_summary_by_sample.tsv", sample_rows)
    write_tsv(RESULTS / "GSE179640_cell_qc_summary_by_location.tsv", tissue_rows)

    total_raw = sum(int(row["n_cells_raw"]) for row in sample_rows)
    total_adaptive = sum(int(row["n_cells_adaptive_qc"]) for row in sample_rows)
    location_counts = Counter(str(row["sample_location"]) for row in sample_rows)
    subject_by_location: dict[str, set[str]] = defaultdict(set)
    for row in sample_rows:
        subject_by_location[str(row["sample_location"])].add(str(row["subject_code"]))

    lines = [
        "# GSE179640 cell-level QC summary",
        "",
        f"- Primary tissue samples processed: {len(sample_rows)}",
        f"- Raw Cell Ranger filtered barcodes: {total_raw}",
        f"- Barcodes passing adaptive tissue QC: {total_adaptive}",
        f"- Adaptive QC retained fraction: {total_adaptive / total_raw:.4f}" if total_raw else "- Adaptive QC retained fraction: NA",
        f"- Sample-location counts: `{json.dumps(location_counts, sort_keys=True)}`",
        f"- Subject counts by sample location: `{json.dumps({k: len(v) for k, v in subject_by_location.items()}, sort_keys=True)}`",
        f"- Cell-level QC metrics: `{metrics_path}`",
        "",
        "## QC definitions",
        "",
        "- Lenient QC: n_features >= 200, n_counts >= 300, pct_mt <= 25.",
        "- Standard QC: 300 <= n_features <= 7500, n_counts >= 500, pct_mt <= 20.",
        "- Adaptive tissue QC: n_features >= 300, n_counts >= 500, pct_mt <= 25, with per-sample upper filters at max(5000, q99.5 n_features) and max(25000, q99.5 n_counts).",
        "",
        "## By-location QC summary",
        "",
        "sample_location\ttissue\tn_samples\tn_cells_raw\tn_cells_adaptive_qc\tadaptive_qc_fraction",
    ]
    for row in tissue_rows:
        lines.append(
            "\t".join(
                [
                    str(row["sample_location"]),
                    str(row["tissue"]),
                    str(row["n_samples"]),
                    str(row["n_cells_raw"]),
                    str(row["n_cells_adaptive_qc"]),
                    f"{float(row['adaptive_qc_fraction']):.4f}",
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
            "- This QC layer uses Cell Ranger filtered barcodes as input and calculates independent per-cell QC metrics before any clustering.",
            "- The adaptive QC gate is intentionally conservative for heterogeneous surgical tissues; it is a proposed analysis gate, not a claim that all other cells are artifacts.",
            "- Doublet detection and cell-type annotation are still missing; downstream cell-state localization must not proceed without those steps or an explicitly conservative alternative.",
            "- Donor/sample metadata are preserved for later pseudobulk and donor-aware summaries.",
            "",
        ]
    )
    (RESULTS / "GSE179640_cell_qc_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print(RESULTS / "GSE179640_cell_qc_summary.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
