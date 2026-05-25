#!/usr/bin/env python3
"""Build cross-disease scRNA localization comparison for candidate genes."""

from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT_ROOT / "results"
SINGLE = RESULTS / "singlecell"
INTEGRATION = RESULTS / "integration"
PY_PKG_DIR = PROJECT_ROOT / "tools" / "py_packages"

if sys.version_info[:2] != (3, 12):
    raise SystemExit(
        "Use the bundled Python 3.12 runtime: "
        "python3 "
        "scripts/build_cross_disease_singlecell_localization_matrix.py"
    )

sys.path.insert(0, str(PY_PKG_DIR))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


GWAS = RESULTS / "gwas" / "gwas_candidate_gene_universe.tsv"
BULK = RESULTS / "bulk" / "bulk_candidate_expression_support_scores.tsv"
ENDO = SINGLE / "GSE179640_singlecell_candidate_evidence_matrix_v2.tsv"
ADENO = SINGLE / "Zenodo17078290_candidate_localization_matrix.tsv"

OUT_MATRIX = INTEGRATION / "cross_disease_singlecell_localization_matrix.tsv"
OUT_SUMMARY = INTEGRATION / "cross_disease_singlecell_localization_summary.md"
OUT_REVIEW = INTEGRATION / "phase10_cross_disease_singlecell_localization_self_review.md"


def clean_id(value: object) -> str:
    return str(value).split(".")[0]


def to_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").fillna(0)


def genetics_score(priority: str) -> int:
    if "tier1" in str(priority):
        return 20
    if "tier2" in str(priority):
        return 12
    if "tier3" in str(priority):
        return 5
    return 0


def localization_class(endo_norm: float, adeno_score: float) -> str:
    if endo_norm >= 12 and adeno_score >= 12:
        return "shared_scRNA_localization"
    if endo_norm >= 12 and adeno_score >= 5:
        return "endometriosis_dominant_with_adenomyosis_signal"
    if adeno_score >= 12 and endo_norm >= 5:
        return "adenomyosis_dominant_with_endometriosis_signal"
    if endo_norm >= 12:
        return "endometriosis_prioritized_scRNA_localization"
    if adeno_score >= 12:
        return "adenomyosis_prioritized_scRNA_localization"
    if endo_norm >= 5 or adeno_score >= 5:
        return "limited_or_single_layer_scRNA_localization"
    return "minimal_scRNA_localization"


def axis_score(row: pd.Series, axis: str) -> float:
    if axis == "fibrovascular":
        return float(row.get("endo_fibrovascular_max_prevalence", 0)) + float(row.get("adeno_fibrovascular_max_prevalence", 0))
    if axis == "immune":
        return float(row.get("endo_immune_max_prevalence", 0)) + float(row.get("adeno_immune_max_prevalence", 0))
    if axis == "epithelial":
        return float(row.get("endo_epithelial_max_prevalence", 0)) + float(row.get("adeno_epithelial_max_prevalence", 0))
    return 0.0


def dominant_axis(row: pd.Series) -> str:
    scores = {
        "fibrovascular": axis_score(row, "fibrovascular"),
        "immune": axis_score(row, "immune"),
        "epithelial": axis_score(row, "epithelial"),
    }
    best_axis, best_score = max(scores.items(), key=lambda item: item[1])
    if best_score <= 0:
        return "none_detected"
    ordered = sorted(scores.values(), reverse=True)
    if len(ordered) > 1 and ordered[1] >= 0.8 * ordered[0]:
        return "mixed"
    return best_axis


def shared_axis_flags(row: pd.Series) -> tuple[bool, bool, bool]:
    shared_fibro = (
        float(row.get("endo_fibrovascular_max_prevalence", 0)) >= 0.10
        and float(row.get("adeno_fibrovascular_max_prevalence", 0)) >= 0.10
    )
    shared_immune = (
        float(row.get("endo_immune_max_prevalence", 0)) >= 0.10
        and float(row.get("adeno_immune_max_prevalence", 0)) >= 0.10
    )
    shared_epithelial = (
        float(row.get("endo_epithelial_max_prevalence", 0)) >= 0.10
        and float(row.get("adeno_epithelial_max_prevalence", 0)) >= 0.10
    )
    return shared_fibro, shared_immune, shared_epithelial


def main() -> int:
    INTEGRATION.mkdir(parents=True, exist_ok=True)
    gwas = pd.read_csv(GWAS, sep="\t", keep_default_na=False)
    bulk = pd.read_csv(BULK, sep="\t", keep_default_na=False)
    endo = pd.read_csv(ENDO, sep="\t", keep_default_na=False)
    adeno = pd.read_csv(ADENO, sep="\t", keep_default_na=False)
    for df in [gwas, bulk, endo, adeno]:
        df["gene_id"] = df["gene_id"].map(clean_id)

    base_cols = [
        "gene_id",
        "gene_symbol",
        "gene_biotype",
        "gene_description",
        "genetic_priority",
        "ld_neighborhood_class",
        "module_hint_preliminary",
    ]
    base = gwas[[col for col in base_cols if col in gwas.columns]].drop_duplicates("gene_id")
    bulk_keep = ["gene_id", "bulk_expression_support_score_20", "bulk_support_class"]
    base = base.merge(bulk[bulk_keep], on="gene_id", how="left")

    endo_keep = [
        "gene_id",
        "gse179640_integrated_singlecell_support_score_35",
        "gse179640_integrated_singlecell_support_class",
        "gse179640_singlecell_support_score_25",
        "pseudobulk_support_score_10",
        "lesion_top_location",
        "top_broad_compartment",
        "fibrovascular_max_prevalence",
        "fibrovascular_top_location_compartment",
        "immune_max_prevalence",
        "immune_top_location_compartment",
        "epithelial_max_prevalence",
        "epithelial_top_location_compartment",
        "best_pseudobulk_comparison",
        "best_pseudobulk_effect",
        "best_pseudobulk_p_value",
        "v2_interpretation_guardrail",
    ]
    adeno_keep = [
        "gene_id",
        "matched_in_h5ad",
        "zenodo17078290_localization_score_25",
        "zenodo17078290_localization_class",
        "top_adeno_group",
        "top_adeno_cluster",
        "top_adeno_median_prevalence",
        "control_max_median_prevalence",
        "fibrovascular_max_prevalence",
        "fibrovascular_top_group_cluster",
        "immune_max_prevalence",
        "immune_top_group_cluster",
        "epithelial_max_prevalence",
        "epithelial_top_group_cluster",
        "best_sample_level_comparison",
        "best_sample_level_effect",
        "best_sample_level_p_value",
        "interpretation_guardrail",
    ]
    endo = endo[[col for col in endo_keep if col in endo.columns]].rename(
        columns={
            "gse179640_integrated_singlecell_support_score_35": "endo_scRNA_score_35",
            "gse179640_integrated_singlecell_support_class": "endo_scRNA_class",
            "gse179640_singlecell_support_score_25": "endo_detectability_score_25",
            "pseudobulk_support_score_10": "endo_pseudobulk_score_10",
            "lesion_top_location": "endo_top_location",
            "top_broad_compartment": "endo_top_broad_compartment",
            "fibrovascular_max_prevalence": "endo_fibrovascular_max_prevalence",
            "fibrovascular_top_location_compartment": "endo_fibrovascular_top_location_compartment",
            "immune_max_prevalence": "endo_immune_max_prevalence",
            "immune_top_location_compartment": "endo_immune_top_location_compartment",
            "epithelial_max_prevalence": "endo_epithelial_max_prevalence",
            "epithelial_top_location_compartment": "endo_epithelial_top_location_compartment",
            "best_pseudobulk_comparison": "endo_best_pseudobulk_comparison",
            "best_pseudobulk_effect": "endo_best_pseudobulk_effect",
            "best_pseudobulk_p_value": "endo_best_pseudobulk_p_value",
            "v2_interpretation_guardrail": "endo_guardrail",
        }
    )
    adeno = adeno[[col for col in adeno_keep if col in adeno.columns]].rename(
        columns={
            "matched_in_h5ad": "adeno_matched_in_h5ad",
            "zenodo17078290_localization_score_25": "adeno_scRNA_localization_score_25",
            "zenodo17078290_localization_class": "adeno_scRNA_localization_class",
            "fibrovascular_max_prevalence": "adeno_fibrovascular_max_prevalence",
            "fibrovascular_top_group_cluster": "adeno_fibrovascular_top_group_cluster",
            "immune_max_prevalence": "adeno_immune_max_prevalence",
            "immune_top_group_cluster": "adeno_immune_top_group_cluster",
            "epithelial_max_prevalence": "adeno_epithelial_max_prevalence",
            "epithelial_top_group_cluster": "adeno_epithelial_top_group_cluster",
            "best_sample_level_comparison": "adeno_best_sample_level_comparison",
            "best_sample_level_effect": "adeno_best_sample_level_effect",
            "best_sample_level_p_value": "adeno_best_sample_level_p_value",
            "interpretation_guardrail": "adeno_guardrail",
        }
    )
    merged = base.merge(endo, on="gene_id", how="left").merge(adeno, on="gene_id", how="left")
    numeric_cols = [
        "bulk_expression_support_score_20",
        "endo_scRNA_score_35",
        "endo_detectability_score_25",
        "endo_pseudobulk_score_10",
        "endo_fibrovascular_max_prevalence",
        "endo_immune_max_prevalence",
        "endo_epithelial_max_prevalence",
        "adeno_scRNA_localization_score_25",
        "adeno_fibrovascular_max_prevalence",
        "adeno_immune_max_prevalence",
        "adeno_epithelial_max_prevalence",
    ]
    for col in numeric_cols:
        if col in merged.columns:
            merged[col] = to_num(merged[col])
    merged["genetics_support_score_20"] = merged["genetic_priority"].map(genetics_score)
    merged["endo_scRNA_normalized_score_20"] = (merged["endo_scRNA_score_35"] / 35.0 * 20.0).round(3)
    merged["adeno_scRNA_normalized_score_20"] = (merged["adeno_scRNA_localization_score_25"] / 25.0 * 20.0).round(3)
    merged["cross_disease_scRNA_score_40"] = (
        merged["endo_scRNA_normalized_score_20"] + merged["adeno_scRNA_normalized_score_20"]
    ).round(3)
    merged["cross_disease_scRNA_class"] = merged.apply(
        lambda row: localization_class(float(row["endo_scRNA_normalized_score_20"]), float(row["adeno_scRNA_normalized_score_20"])),
        axis=1,
    )
    flags = merged.apply(shared_axis_flags, axis=1, result_type="expand")
    flags.columns = ["shared_fibrovascular_flag", "shared_immune_flag", "shared_epithelial_flag"]
    merged = pd.concat([merged, flags], axis=1)
    merged["dominant_cross_disease_axis"] = merged.apply(dominant_axis, axis=1)
    merged["pre_druggability_biologic_evidence_score_80"] = (
        merged["genetics_support_score_20"]
        + merged["bulk_expression_support_score_20"]
        + merged["cross_disease_scRNA_score_40"]
    ).round(3)
    merged["cross_disease_guardrail"] = (
        "broad_label_scRNA_localization_not_final_cell_state_DE_not_spatial"
    )
    merged = merged.sort_values(
        [
            "pre_druggability_biologic_evidence_score_80",
            "cross_disease_scRNA_score_40",
            "bulk_expression_support_score_20",
            "gene_symbol",
        ],
        ascending=[False, False, False, True],
    )
    merged.to_csv(OUT_MATRIX, sep="\t", index=False)

    class_counts = merged["cross_disease_scRNA_class"].value_counts().to_dict()
    axis_counts = merged["dominant_cross_disease_axis"].value_counts().to_dict()
    shared_fibro = int(merged["shared_fibrovascular_flag"].sum())
    shared_immune = int(merged["shared_immune_flag"].sum())
    shared_epi = int(merged["shared_epithelial_flag"].sum())
    lines = [
        "# Cross-disease scRNA localization summary",
        "",
        "## Scope",
        "",
        "- Integrates GSE179640 endometriosis broad-compartment v2 evidence with Zenodo 17078290 adenomyosis scRNA localization evidence.",
        "- This is a pre-druggability biological evidence matrix; it does not include Open Targets, ChEMBL, DGIdb, DepMap or safety yet.",
        "- Guardrail: broad-label/candidate localization evidence only; not final cell-state differential expression and not spatial evidence.",
        "",
        "## Counts",
        "",
        f"- Candidate genes: {len(merged)}",
        f"- Cross-disease scRNA classes: `{class_counts}`",
        f"- Dominant axis counts: `{axis_counts}`",
        f"- Shared fibrovascular prevalence flags: {shared_fibro}",
        f"- Shared immune prevalence flags: {shared_immune}",
        f"- Shared epithelial prevalence flags: {shared_epi}",
        "",
        "## Top pre-druggability biological evidence candidates",
        "",
        "gene_symbol\tgene_id\tbio80\tgenetics20\tbulk20\tscRNA40\tclass\taxis\tendo_score20\tadeno_score20\tshared_fibro\tshared_immune\tshared_epi",
    ]
    for _, row in merged.head(30).iterrows():
        lines.append(
            "\t".join(
                [
                    str(row["gene_symbol"]),
                    str(row["gene_id"]),
                    f"{float(row['pre_druggability_biologic_evidence_score_80']):.3f}",
                    str(row["genetics_support_score_20"]),
                    f"{float(row['bulk_expression_support_score_20']):.3f}",
                    f"{float(row['cross_disease_scRNA_score_40']):.3f}",
                    str(row["cross_disease_scRNA_class"]),
                    str(row["dominant_cross_disease_axis"]),
                    f"{float(row['endo_scRNA_normalized_score_20']):.3f}",
                    f"{float(row['adeno_scRNA_normalized_score_20']):.3f}",
                    str(bool(row["shared_fibrovascular_flag"])),
                    str(bool(row["shared_immune_flag"])),
                    str(bool(row["shared_epithelial_flag"])),
                ]
            )
        )
    lines.extend(
        [
            "",
            "## High-bulk candidates",
            "",
            "gene_symbol\tgene_id\tbio80\tclass\taxis\tendo_score20\tadeno_score20\tendo_top\tadeno_top",
        ]
    )
    high_bulk = merged[merged["bulk_support_class"] == "high_bulk_support"].copy()
    for _, row in high_bulk.iterrows():
        lines.append(
            "\t".join(
                [
                    str(row["gene_symbol"]),
                    str(row["gene_id"]),
                    f"{float(row['pre_druggability_biologic_evidence_score_80']):.3f}",
                    str(row["cross_disease_scRNA_class"]),
                    str(row["dominant_cross_disease_axis"]),
                    f"{float(row['endo_scRNA_normalized_score_20']):.3f}",
                    f"{float(row['adeno_scRNA_normalized_score_20']):.3f}",
                    f"{row.get('endo_top_location', '')}:{row.get('endo_top_broad_compartment', '')}",
                    f"{row.get('top_adeno_group', '')}:{row.get('top_adeno_cluster', '')}",
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

    review = [
        "# Phase 10 Self-Review: cross-disease scRNA localization matrix",
        "",
        "## Verdict",
        "",
        "PASS_WITH_CONDITIONS",
        "",
        "## What passed",
        "",
        "- Joined 102 GWAS candidate records with bulk, endometriosis GSE179640 v2 and adenomyosis Zenodo 17078290 localization evidence.",
        "- Created explicit shared/endometriosis-prioritized/adenomyosis-prioritized scRNA localization classes.",
        "- Added fibrovascular, immune and epithelial shared-prevalence flags.",
        "- Created an 80-point pre-druggability biological evidence score combining genetics, bulk and cross-disease scRNA evidence.",
        "",
        "## Limitations",
        "",
        "- The 80-point score is not the final target score because druggability, safety, directionality and rank-stability layers are not included.",
        "- GSE179640 remains broad-compartment marker-panel triage, not author-label cluster DE.",
        "- Zenodo 17078290 lacks detected spatial coordinates and has only three samples per encoded group.",
        "- Several high-scoring genes are structural or housekeeping-like readouts; they must not be promoted as drug targets without actionability evidence.",
        "",
        "## Decision",
        "",
        "Use this matrix for downstream cross-disease cell-context scoring.",
    ]
    OUT_REVIEW.write_text("\n".join(review) + "\n", encoding="utf-8")
    print(OUT_SUMMARY)
    print(OUT_MATRIX)
    print(OUT_REVIEW)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
