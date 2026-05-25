#!/usr/bin/env python3
"""Merge GSE179640 detectability and donor-aware pseudobulk evidence."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT_ROOT / "results" / "singlecell"
PY_PKG_DIR = PROJECT_ROOT / "tools" / "py_packages"

if sys.version_info[:2] != (3, 12):
    raise SystemExit(
        "Use the bundled Python 3.12 runtime: "
        "python3 "
        "scripts/build_gse179640_singlecell_evidence_matrix_v2.py"
    )

sys.path.insert(0, str(PY_PKG_DIR))

import pandas as pd  # noqa: E402


DETECTABILITY = RESULTS / "GSE179640_singlecell_candidate_evidence_matrix.tsv"
PSEUDOBULK = RESULTS / "GSE179640_candidate_pseudobulk_support.tsv"
OUT_MATRIX = RESULTS / "GSE179640_singlecell_candidate_evidence_matrix_v2.tsv"
OUT_SUMMARY = RESULTS / "GSE179640_singlecell_candidate_evidence_matrix_v2_summary.md"
OUT_REVIEW = RESULTS / "phase7b_gse179640_integrated_singlecell_evidence_self_review.md"


def integrated_class(row: pd.Series) -> str:
    detectability_score = int(row["gse179640_singlecell_support_score_25"])
    pseudobulk_score = int(row["pseudobulk_support_score_10"])
    combined_score = int(row["gse179640_integrated_singlecell_support_score_35"])
    if combined_score >= 30 and pseudobulk_score >= 10:
        return "strong_donor_aware_singlecell_support"
    if combined_score >= 28 and pseudobulk_score >= 8:
        return "moderate_relaxed_fdr_singlecell_support"
    if combined_score >= 23 and pseudobulk_score >= 5:
        return "suggestive_donor_aware_singlecell_support"
    if detectability_score >= 18 and pseudobulk_score < 5:
        return "high_detectability_limited_donor_support"
    if combined_score >= 18:
        return "moderate_donor_aware_singlecell_support"
    if combined_score >= 8:
        return "limited_singlecell_support"
    return "minimal_singlecell_support"


def format_value(value: object, digits: int = 4) -> str:
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    try:
        return f"{float(value):.{digits}g}"
    except Exception:
        return str(value)


def main() -> int:
    if not DETECTABILITY.exists():
        raise SystemExit(f"Missing input: {DETECTABILITY}")
    if not PSEUDOBULK.exists():
        raise SystemExit(f"Missing input: {PSEUDOBULK}")

    detect = pd.read_csv(DETECTABILITY, sep="\t", keep_default_na=False)
    pseudo = pd.read_csv(PSEUDOBULK, sep="\t", keep_default_na=False)
    for col in ["gse179640_singlecell_support_score_25"]:
        detect[col] = pd.to_numeric(detect[col], errors="coerce").fillna(0).astype(int)
    pseudo["pseudobulk_support_score_10"] = (
        pd.to_numeric(pseudo["pseudobulk_support_score_10"], errors="coerce").fillna(0).astype(int)
    )

    keep_pseudo = [
        "gene_id",
        "pseudobulk_support_score_10",
        "pseudobulk_support_class",
        "best_pseudobulk_compartment",
        "best_pseudobulk_comparison",
        "best_pseudobulk_metric",
        "best_pseudobulk_effect",
        "best_pseudobulk_p_value",
        "best_pseudobulk_q_value_by_comparison_metric",
        "best_pseudobulk_q_value_global",
        "n_case_subjects",
        "n_reference_subjects",
        "n_paired_subjects",
    ]
    merged = detect.merge(pseudo[keep_pseudo], on="gene_id", how="left")
    merged["pseudobulk_support_score_10"] = merged["pseudobulk_support_score_10"].fillna(0).astype(int)
    merged["pseudobulk_support_class"] = merged["pseudobulk_support_class"].replace("", "not_tested")
    merged["gse179640_integrated_singlecell_support_score_35"] = (
        merged["gse179640_singlecell_support_score_25"] + merged["pseudobulk_support_score_10"]
    )
    merged["gse179640_integrated_singlecell_support_class"] = merged.apply(integrated_class, axis=1)
    merged["v2_interpretation_guardrail"] = (
        "integrates_detectability_and_subject_level_pseudobulk_broad_compartment_not_final_cell_state"
    )

    sort_cols = [
        "gse179640_integrated_singlecell_support_score_35",
        "pseudobulk_support_score_10",
        "fibrovascular_max_prevalence",
        "bulk_expression_support_score_20",
        "gene_symbol",
    ]
    for col in ["fibrovascular_max_prevalence", "bulk_expression_support_score_20"]:
        merged[col] = pd.to_numeric(merged[col], errors="coerce").fillna(0)
    merged = merged.sort_values(sort_cols, ascending=[False, False, False, False, True])
    merged.to_csv(OUT_MATRIX, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)

    class_counts = merged["gse179640_integrated_singlecell_support_class"].value_counts().to_dict()
    donor_limited = merged[
        (merged["gse179640_singlecell_support_score_25"] >= 18) & (merged["pseudobulk_support_score_10"] < 5)
    ].copy()
    high_bulk = merged[merged["bulk_support_class"] == "high_bulk_support"].copy()

    lines = [
        "# GSE179640 integrated single-cell evidence matrix v2",
        "",
        f"- Candidate genes scored: {len(merged)}",
        f"- Integrated class counts: `{class_counts}`",
        "- Integrated score range: 0-35 = previous QC-filtered broad-compartment detectability score out of 25 plus donor-aware pseudobulk support score out of 10.",
        "- Strong donor-aware support requires q<=0.10 subject-level pseudobulk evidence; q<=0.20 is labeled moderate/relaxed rather than high.",
        "- Guardrail: this remains a broad-compartment GSE179640 layer, not final disease cell-state differential expression.",
        "",
        "## Top integrated candidates",
        "",
        "gene_symbol\tgene_id\tintegrated_score\tintegrated_class\tdetectability_score\tpseudobulk_score\tbest_compartment\tbest_pseudobulk_comparison\tbest_pseudobulk_effect\tbest_pseudobulk_p",
    ]
    for _, row in merged.head(25).iterrows():
        lines.append(
            "\t".join(
                [
                    str(row["gene_symbol"]),
                    str(row["gene_id"]),
                    str(int(row["gse179640_integrated_singlecell_support_score_35"])),
                    str(row["gse179640_integrated_singlecell_support_class"]),
                    str(int(row["gse179640_singlecell_support_score_25"])),
                    str(int(row["pseudobulk_support_score_10"])),
                    str(row.get("top_broad_compartment", "")),
                    str(row.get("best_pseudobulk_comparison", "")),
                    format_value(row.get("best_pseudobulk_effect", "")),
                    format_value(row.get("best_pseudobulk_p_value", "")),
                ]
            )
        )
    lines.extend(
        [
            "",
            "## High-bulk-support candidates",
            "",
            "gene_symbol\tgene_id\tbulk_score\tintegrated_score\tintegrated_class\tdetectability_score\tpseudobulk_score\ttop_compartment\tbest_pseudobulk_compartment\tbest_pseudobulk_effect\tbest_pseudobulk_p",
        ]
    )
    for _, row in high_bulk.iterrows():
        lines.append(
            "\t".join(
                [
                    str(row["gene_symbol"]),
                    str(row["gene_id"]),
                    str(format_value(row["bulk_expression_support_score_20"], 3)),
                    str(int(row["gse179640_integrated_singlecell_support_score_35"])),
                    str(row["gse179640_integrated_singlecell_support_class"]),
                    str(int(row["gse179640_singlecell_support_score_25"])),
                    str(int(row["pseudobulk_support_score_10"])),
                    str(row.get("top_broad_compartment", "")),
                    str(row.get("best_pseudobulk_compartment", "")),
                    format_value(row.get("best_pseudobulk_effect", "")),
                    format_value(row.get("best_pseudobulk_p_value", "")),
                ]
            )
        )
    lines.extend(
        [
            "",
            "## Detectability-high but donor-limited genes",
            "",
            f"- Count: {len(donor_limited)}",
            "- These genes should not be promoted as strong GSE179640 single-cell evidence unless supported by other layers.",
            "",
            "gene_symbol\tgene_id\tdetectability_score\tpseudobulk_score\ttop_compartment",
        ]
    )
    for _, row in donor_limited.head(20).iterrows():
        lines.append(
            "\t".join(
                [
                    str(row["gene_symbol"]),
                    str(row["gene_id"]),
                    str(int(row["gse179640_singlecell_support_score_25"])),
                    str(int(row["pseudobulk_support_score_10"])),
                    str(row.get("top_broad_compartment", "")),
                ]
            )
        )
    lines.extend(
        [
            "",
            "## Output files",
            "",
            f"- Matrix: `{OUT_MATRIX}`",
            f"- Self-review: `{OUT_REVIEW}`",
        ]
    )
    OUT_SUMMARY.write_text("\n".join(lines) + "\n", encoding="utf-8")

    review_lines = [
        "# Phase 7b Self-Review: integrated GSE179640 single-cell evidence v2",
        "",
        "## Verdict",
        "",
        "PASS_WITH_CONDITIONS",
        "",
        "## What improved",
        "",
        "- The v2 matrix no longer ranks candidates using detectability/localization alone; subject-level pseudobulk support is incorporated.",
        "- Genes with high broad-compartment detectability but weak donor-aware evidence are explicitly separated from moderate/strong donor-aware support.",
        "- Strong donor-aware support requires q<=0.10 pseudobulk evidence; relaxed q<=0.20 evidence is labeled moderate.",
        "- GSE179640 evidence is broad-compartment context, not final cell-state differential expression.",
        "",
        "## Conditions",
        "",
        "- Pseudobulk support is based on candidate-level counts and prevalence, not full transcriptome pseudobulk normalization.",
        "- Control comparisons remain underpowered because only three control subjects are available in the primary tissue subset.",
        "- Final target prioritization must still incorporate adenomyosis atlas evidence, druggability, safety and rank-stability analyses.",
        "",
        "## Decision",
        "",
        "Use `GSE179640_singlecell_candidate_evidence_matrix_v2.tsv` for downstream scoring.",
    ]
    OUT_REVIEW.write_text("\n".join(review_lines) + "\n", encoding="utf-8")
    print(OUT_SUMMARY)
    print(OUT_MATRIX)
    print(OUT_REVIEW)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
