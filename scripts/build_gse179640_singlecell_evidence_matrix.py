#!/usr/bin/env python3
"""Build a conservative GSE179640 single-cell evidence matrix for candidate genes."""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT_ROOT / "results" / "singlecell"
BULK = PROJECT_ROOT / "results" / "bulk" / "bulk_candidate_expression_support_scores.tsv"
GWAS = PROJECT_ROOT / "results" / "gwas" / "gwas_candidate_gene_universe.tsv"

LOCATION_SUMMARY = RESULTS / "GSE179640_candidate_expression_summary_by_location_adaptive_qc.tsv"
COMPARTMENT_SUMMARY = RESULTS / "GSE179640_candidate_expression_summary_by_broad_compartment.tsv"
OUT_MATRIX = RESULTS / "GSE179640_singlecell_candidate_evidence_matrix.tsv"
OUT_SUMMARY = RESULTS / "GSE179640_singlecell_candidate_evidence_matrix_summary.md"

LESION_LOCATIONS = ["Ectopic", "Ectopic Adjacent", "Ectopic Ovary"]
PRIMARY_LOCATIONS = ["Eutopic", *LESION_LOCATIONS]
INTERPRETABLE_COMPARTMENTS = [
    "epithelial",
    "stromal_fibroblast",
    "endothelial",
    "mural_smooth_muscle",
    "t_nk",
    "myeloid_macrophage",
    "b_plasma",
    "mast",
]
FIBROVASCULAR_COMPARTMENTS = {"stromal_fibroblast", "endothelial", "mural_smooth_muscle"}
IMMUNE_COMPARTMENTS = {"t_nk", "myeloid_macrophage", "b_plasma", "mast"}
EPITHELIAL_COMPARTMENTS = {"epithelial"}


def load_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def to_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def to_int(value: object) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def prevalence_band(value: float) -> str:
    if value >= 0.20:
        return "high_detectability"
    if value >= 0.10:
        return "moderate_detectability"
    if value >= 0.05:
        return "low_to_moderate_detectability"
    if value > 0:
        return "low_detectability"
    return "not_detected"


def location_detectability_score(max_prev: float) -> int:
    if max_prev >= 0.20:
        return 5
    if max_prev >= 0.10:
        return 4
    if max_prev >= 0.05:
        return 3
    if max_prev > 0:
        return 1
    return 0


def robustness_score(row: pd.Series | None) -> int:
    if row is None:
        return 0
    detection_fraction = to_float(row.get("sample_detection_fraction_qc", 0))
    n_samples = to_int(row.get("n_samples", 0))
    if detection_fraction >= 0.80 and n_samples >= 4:
        return 5
    if detection_fraction >= 0.50 and n_samples >= 3:
        return 3
    if detection_fraction > 0:
        return 1
    return 0


def compartment_score(row: pd.Series | None) -> int:
    if row is None:
        return 0
    median_prev = to_float(row.get("median_prevalence", 0))
    median_cells = to_float(row.get("median_compartment_cells", 0))
    n_samples = to_int(row.get("n_samples_with_compartment", 0))
    if median_prev >= 0.10 and median_cells >= 50 and n_samples >= 3:
        return 5
    if median_prev >= 0.05 and median_cells >= 20 and n_samples >= 2:
        return 3
    if median_prev > 0 and median_cells >= 10:
        return 1
    return 0


def compartment_group_score(max_prev: float) -> int:
    if max_prev >= 0.20:
        return 5
    if max_prev >= 0.10:
        return 4
    if max_prev >= 0.05:
        return 3
    if max_prev > 0:
        return 1
    return 0


def consistency_score(location_prev: dict[str, float]) -> int:
    n_supported = sum(1 for loc in LESION_LOCATIONS if location_prev.get(loc, 0) >= 0.05)
    if n_supported == 3:
        return 5
    if n_supported == 2:
        return 3
    if n_supported == 1:
        return 1
    return 0


def evidence_class(score: int) -> str:
    if score >= 18:
        return "high_singlecell_support"
    if score >= 11:
        return "moderate_singlecell_support"
    if score >= 5:
        return "limited_singlecell_support"
    return "minimal_singlecell_support"


def first_row(df: pd.DataFrame) -> pd.Series | None:
    if df.empty:
        return None
    return df.iloc[0]


def top_compartment_row(comp: pd.DataFrame) -> pd.Series | None:
    eligible = comp[
        comp["sample_location"].isin(LESION_LOCATIONS)
        & comp["broad_compartment"].isin(INTERPRETABLE_COMPARTMENTS)
    ].copy()
    if eligible.empty:
        return None
    eligible["_prev"] = eligible["median_prevalence"].map(to_float)
    eligible["_cells"] = eligible["median_compartment_cells"].map(to_float)
    eligible["_n_samples"] = eligible["n_samples_with_compartment"].map(to_int)
    eligible = eligible[(eligible["_cells"] >= 20) & (eligible["_n_samples"] >= 2)]
    if eligible.empty:
        return None
    eligible = eligible.sort_values(["_prev", "_cells", "_n_samples"], ascending=[False, False, False])
    return eligible.iloc[0]


def max_compartment_prev(comp: pd.DataFrame, compartments: set[str]) -> tuple[float, str]:
    eligible = comp[
        comp["sample_location"].isin(LESION_LOCATIONS)
        & comp["broad_compartment"].isin(compartments)
    ].copy()
    if eligible.empty:
        return 0.0, ""
    eligible["_prev"] = eligible["median_prevalence"].map(to_float)
    eligible["_cells"] = eligible["median_compartment_cells"].map(to_float)
    eligible["_n_samples"] = eligible["n_samples_with_compartment"].map(to_int)
    eligible = eligible[(eligible["_cells"] >= 20) & (eligible["_n_samples"] >= 2)]
    if eligible.empty:
        return 0.0, ""
    eligible = eligible.sort_values(["_prev", "_cells"], ascending=[False, False])
    row = eligible.iloc[0]
    return to_float(row["_prev"]), f"{row['sample_location']}:{row['broad_compartment']}"


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    gwas = load_tsv(GWAS)
    bulk = load_tsv(BULK)
    loc = load_tsv(LOCATION_SUMMARY)
    comp = load_tsv(COMPARTMENT_SUMMARY)

    bulk_by_gene = {row["gene_id"].split(".")[0]: row for _, row in bulk.iterrows()}
    loc["gene_id_clean"] = loc["gene_id"].map(lambda x: str(x).split(".")[0])
    comp["gene_id_clean"] = comp["gene_id"].map(lambda x: str(x).split(".")[0])

    rows: list[dict[str, object]] = []
    for _, gene in gwas.iterrows():
        gene_id = str(gene["gene_id"]).split(".")[0]
        gene_loc = loc[loc["gene_id_clean"] == gene_id]
        gene_comp = comp[comp["gene_id_clean"] == gene_id]
        bulk_row = bulk_by_gene.get(gene_id, {})

        location_prev = {
            location: to_float(
                first_row(gene_loc[gene_loc["sample_location"] == location]).get("median_prevalence_qc", 0)
                if first_row(gene_loc[gene_loc["sample_location"] == location]) is not None
                else 0
            )
            for location in PRIMARY_LOCATIONS
        }
        lesion_location_rows = gene_loc[gene_loc["sample_location"].isin(LESION_LOCATIONS)].copy()
        lesion_location_rows["_prev"] = lesion_location_rows["median_prevalence_qc"].map(to_float)
        lesion_location_rows = lesion_location_rows.sort_values("_prev", ascending=False)
        top_loc = first_row(lesion_location_rows)
        max_lesion_prev = to_float(top_loc.get("_prev", 0)) if top_loc is not None else 0.0

        top_comp = top_compartment_row(gene_comp)
        max_fibrovascular, top_fibrovascular = max_compartment_prev(gene_comp, FIBROVASCULAR_COMPARTMENTS)
        max_immune, top_immune = max_compartment_prev(gene_comp, IMMUNE_COMPARTMENTS)
        max_epithelial, top_epithelial = max_compartment_prev(gene_comp, EPITHELIAL_COMPARTMENTS)

        loc_score = location_detectability_score(max_lesion_prev)
        robust_score = robustness_score(top_loc)
        comp_score = compartment_score(top_comp)
        fibro_score = compartment_group_score(max_fibrovascular)
        consist_score = consistency_score(location_prev)
        singlecell_score = loc_score + robust_score + comp_score + fibro_score + consist_score

        rows.append(
            {
                "gene_id": gene_id,
                "gene_symbol": gene.get("gene_symbol", ""),
                "genetic_priority": gene.get("genetic_priority", ""),
                "ld_neighborhood_class": gene.get("ld_neighborhood_class", ""),
                "module_hint_preliminary": gene.get("module_hint_preliminary", ""),
                "gene_description": gene.get("gene_description", ""),
                "bulk_expression_support_score_20": bulk_row.get("bulk_expression_support_score_20", ""),
                "bulk_support_class": bulk_row.get("bulk_support_class", ""),
                "lesion_top_location": top_loc.get("sample_location", "") if top_loc is not None else "",
                "lesion_top_location_median_prevalence_qc": max_lesion_prev,
                "lesion_detectability_band": prevalence_band(max_lesion_prev),
                "eutopic_median_prevalence_qc": location_prev.get("Eutopic", 0),
                "ectopic_peritoneum_median_prevalence_qc": location_prev.get("Ectopic", 0),
                "ectopic_adjacent_median_prevalence_qc": location_prev.get("Ectopic Adjacent", 0),
                "ectopic_ovary_median_prevalence_qc": location_prev.get("Ectopic Ovary", 0),
                "top_broad_compartment": top_comp.get("broad_compartment", "") if top_comp is not None else "",
                "top_broad_compartment_location": top_comp.get("sample_location", "") if top_comp is not None else "",
                "top_broad_compartment_median_prevalence": to_float(top_comp.get("median_prevalence", 0)) if top_comp is not None else 0,
                "top_broad_compartment_median_cells": to_float(top_comp.get("median_compartment_cells", 0)) if top_comp is not None else 0,
                "fibrovascular_max_prevalence": max_fibrovascular,
                "fibrovascular_top_location_compartment": top_fibrovascular,
                "immune_max_prevalence": max_immune,
                "immune_top_location_compartment": top_immune,
                "epithelial_max_prevalence": max_epithelial,
                "epithelial_top_location_compartment": top_epithelial,
                "singlecell_location_score_5": loc_score,
                "singlecell_sample_robustness_score_5": robust_score,
                "singlecell_compartment_score_5": comp_score,
                "singlecell_fibrovascular_score_5": fibro_score,
                "singlecell_cross_lesion_location_score_5": consist_score,
                "gse179640_singlecell_support_score_25": singlecell_score,
                "gse179640_singlecell_support_class": evidence_class(singlecell_score),
                "interpretation_guardrail": "broad_compartment_triage_not_final_cell_state",
            }
        )

    out = pd.DataFrame(rows)
    out = out.sort_values(
        [
            "gse179640_singlecell_support_score_25",
            "fibrovascular_max_prevalence",
            "lesion_top_location_median_prevalence_qc",
            "gene_symbol",
        ],
        ascending=[False, False, False, True],
    )
    out.to_csv(OUT_MATRIX, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)

    class_counts = out["gse179640_singlecell_support_class"].value_counts().to_dict()
    top_overall = out.head(20)
    high_bulk = out[out["bulk_support_class"] == "high_bulk_support"].copy()
    high_bulk = high_bulk.sort_values(
        ["gse179640_singlecell_support_score_25", "fibrovascular_max_prevalence"],
        ascending=[False, False],
    )

    lines = [
        "# GSE179640 single-cell evidence matrix",
        "",
        f"- Candidate genes scored: {len(out)}",
        f"- Single-cell support class counts: `{class_counts}`",
        "- Score range: 0-25, using QC-filtered lesion detectability, sample robustness, broad-compartment localization, fibrovascular compartment prevalence and cross-lesion-location consistency.",
        "- This score is only the endometriosis GSE179640 single-cell layer; it is not the final target score.",
        "",
        "## Top candidates by GSE179640 single-cell support",
        "",
        "gene_symbol\tgene_id\tsc_score\tbulk_score\ttop_location\ttop_compartment\tfibrovascular_prev\tguardrail",
    ]
    for _, row in top_overall.iterrows():
        lines.append(
            "\t".join(
                [
                    str(row["gene_symbol"]),
                    str(row["gene_id"]),
                    str(row["gse179640_singlecell_support_score_25"]),
                    str(row["bulk_expression_support_score_20"]),
                    str(row["lesion_top_location"]),
                    str(row["top_broad_compartment"]),
                    f"{to_float(row['fibrovascular_max_prevalence']):.4f}",
                    str(row["interpretation_guardrail"]),
                ]
            )
        )
    lines.extend(
        [
            "",
            "## High-bulk-support candidates in the GSE179640 single-cell layer",
            "",
            "gene_symbol\tgene_id\tsc_score\tsc_class\ttop_compartment\tfibrovascular_prev\timmune_prev",
        ]
    )
    for _, row in high_bulk.iterrows():
        lines.append(
            "\t".join(
                [
                    str(row["gene_symbol"]),
                    str(row["gene_id"]),
                    str(row["gse179640_singlecell_support_score_25"]),
                    str(row["gse179640_singlecell_support_class"]),
                    str(row["top_broad_compartment"]),
                    f"{to_float(row['fibrovascular_max_prevalence']):.4f}",
                    f"{to_float(row['immune_max_prevalence']):.4f}",
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
            "- Strength: the matrix joins GWAS, bulk and single-cell evidence with explicit sample-aware QC-filtered summaries.",
            "- Strength: broad-compartment localization is constrained to interpretable marker compartments with minimum cell/sample support.",
            "- Limitation: no doublet-removal model, graph clustering, author-label reconciliation or donor-aware cell-type DE has been completed.",
            "- Limitation: high expression in structural genes such as `HSPG2` can indicate niche localization but does not by itself make a drug target.",
            "- Required next gate: compare these broad-compartment signals with adenomyosis h5ad once downloaded and with cluster/marker validation before final target prioritization.",
            "",
        ]
    )
    OUT_SUMMARY.write_text("\n".join(lines), encoding="utf-8")
    print(OUT_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
