#!/usr/bin/env python3
"""Create a transparent candidate-level bulk expression support score."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
BULK = PROJECT_ROOT / "results" / "bulk"
MODELS = BULK / "models"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def fdr_score(adj_p: float | None, p: float | None, high: int, mid: int, low: int, abs_logfc: float | None = None) -> int:
    if abs_logfc is not None and (pd.isna(abs_logfc) or abs(abs_logfc) < 0.25):
        return 0
    if pd.notna(adj_p) and adj_p < 0.10:
        return high
    if pd.notna(p) and p < 0.05:
        return mid
    if pd.notna(p) and p < 0.10:
        return low
    return 0


def best_by_abs_effect(df: pd.DataFrame, contrasts: list[str]) -> pd.Series | None:
    sub = df[(df["contrast"].isin(contrasts)) & (df["analysis_status"] == "tested")].copy()
    if sub.empty:
        return None
    sub["rank_adj"] = sub["adj.P.Val"].fillna(1.0)
    sub["rank_abs_logfc"] = sub["logFC"].abs().fillna(0.0)
    sub = sub.sort_values(["rank_adj", "rank_abs_logfc"], ascending=[True, False])
    return sub.iloc[0]


def get_one(df: pd.DataFrame, **query) -> pd.Series | None:
    sub = df.copy()
    for key, value in query.items():
        sub = sub[sub[key] == value]
    if sub.empty:
        return None
    return sub.iloc[0]


def dataframe_to_markdown(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for _, row in df.iterrows():
        values = []
        for col in cols:
            value = row[col]
            if isinstance(value, float):
                value = f"{value:.4g}"
            values.append(str(value).replace("|", "\\|"))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def main() -> int:
    candidates = read_tsv(BULK / "candidate_genes_unique.tsv")
    cycle = read_tsv(MODELS / "GSE234354_cycle_candidate_results.tsv")
    immune = read_tsv(MODELS / "GSE313775_th_subset_candidate_results.tsv")
    tissue = read_tsv(MODELS / "GSE141549_tissue_candidate_results.tsv")
    gse519_all = read_tsv(MODELS / "GSE51981_endometrium_candidate_results.tsv")
    gse519_single = read_tsv(MODELS / "GSE51981_endometrium_candidate_results_single_gene_probe_sensitivity.tsv")

    rows: list[dict[str, object]] = []
    for _, cand in candidates.iterrows():
        gene_id = cand["gene_id"]
        symbol = cand["gene_symbol"]
        row: dict[str, object] = cand.to_dict()

        cycle_row = get_one(cycle, gene_id=gene_id)
        cycle_adj = cycle_row["adj.P.Val"] if cycle_row is not None and "adj.P.Val" in cycle_row else np.nan
        row["cycle_adj_p"] = cycle_adj
        row["cycle_driven_flag"] = bool(pd.notna(cycle_adj) and cycle_adj < 0.10)
        if pd.isna(cycle_adj):
            cycle_score = 2
        elif cycle_adj >= 0.10:
            cycle_score = 4
        elif cycle_adj >= 0.05:
            cycle_score = 2
        else:
            cycle_score = 0
        row["cycle_nonconfounded_score_4"] = cycle_score

        tissue_sub = tissue[tissue["gene_symbol"] == symbol]
        lesion_best = best_by_abs_effect(
            tissue_sub,
            [
                "lesion-control_endometrium",
                "lesion-patient_eutopic_endometrium",
                "lesion-patient_peritoneum",
            ],
        )
        if lesion_best is None:
            lesion_score = 0
            row["lesion_top_contrast"] = ""
            row["lesion_logFC"] = np.nan
            row["lesion_adj_p"] = np.nan
        else:
            lesion_score = fdr_score(lesion_best["adj.P.Val"], lesion_best["P.Value"], 6, 4, 2, lesion_best["logFC"])
            row["lesion_top_contrast"] = lesion_best["contrast"]
            row["lesion_logFC"] = lesion_best["logFC"]
            row["lesion_adj_p"] = lesion_best["adj.P.Val"]
        row["lesion_support_score_6"] = lesion_score

        all_row = get_one(gse519_all, gene_symbol=symbol)
        single_row = get_one(gse519_single, gene_symbol=symbol)
        all_score = 0
        if all_row is not None:
            all_score = fdr_score(all_row["adj.P.Val"], all_row["P.Value"], 3, 2, 1, all_row["logFC"])
            row["gse51981_all_logFC"] = all_row["logFC"]
            row["gse51981_all_adj_p"] = all_row["adj.P.Val"]
        else:
            row["gse51981_all_logFC"] = np.nan
            row["gse51981_all_adj_p"] = np.nan
        single_score = 0
        if single_row is not None and single_row.get("analysis_status") == "tested":
            single_score = fdr_score(single_row["adj.P.Val"], single_row["P.Value"], 1, 1, 0, single_row["logFC"])
            row["gse51981_single_logFC"] = single_row["logFC"]
            row["gse51981_single_adj_p"] = single_row["adj.P.Val"]
            row["gse51981_single_probe_direction_consistent"] = (
                pd.notna(row["gse51981_all_logFC"])
                and pd.notna(row["gse51981_single_logFC"])
                and np.sign(row["gse51981_all_logFC"]) == np.sign(row["gse51981_single_logFC"])
            )
        else:
            row["gse51981_single_logFC"] = np.nan
            row["gse51981_single_adj_p"] = np.nan
            row["gse51981_single_probe_direction_consistent"] = False
        validation_score = all_score + (single_score if row["gse51981_single_probe_direction_consistent"] else 0)
        row["independent_endometrium_validation_score_4"] = min(validation_score, 4)

        immune_sub = immune[(immune["gene_id"] == gene_id) & (immune["analysis_status"] == "tested")].copy()
        disease_contrasts = [
            "endometriosis_vs_control_Th1",
            "endometriosis_vs_control_Th1_17",
            "endometriosis_vs_control_Th17",
        ]
        immune_best = best_by_abs_effect(immune_sub, disease_contrasts)
        if immune_best is None:
            immune_score = 0
            row["immune_top_contrast"] = ""
            row["immune_logFC"] = np.nan
            row["immune_adj_p"] = np.nan
        else:
            immune_score = fdr_score(immune_best["adj.P.Val"], immune_best["P.Value"], 4, 2, 1, immune_best["logFC"])
            row["immune_top_contrast"] = immune_best["contrast"]
            row["immune_logFC"] = immune_best["logFC"]
            row["immune_adj_p"] = immune_best["adj.P.Val"]
        row["immune_support_score_4"] = immune_score

        eutopic = get_one(tissue_sub, contrast="patient_eutopic_endometrium-control_endometrium")
        consistency_score = 0
        if eutopic is not None and all_row is not None:
            row["eutopic_vs_control_logFC"] = eutopic["logFC"]
            row["eutopic_vs_control_adj_p"] = eutopic["adj.P.Val"]
            if (
                pd.notna(eutopic["P.Value"])
                and pd.notna(all_row["P.Value"])
                and eutopic["P.Value"] < 0.05
                and all_row["P.Value"] < 0.05
                and np.sign(eutopic["logFC"]) == np.sign(all_row["logFC"])
            ):
                consistency_score = 2
        else:
            row["eutopic_vs_control_logFC"] = np.nan
            row["eutopic_vs_control_adj_p"] = np.nan
        row["cross_bulk_direction_consistency_score_2"] = consistency_score

        score = (
            row["lesion_support_score_6"]
            + row["cycle_nonconfounded_score_4"]
            + row["independent_endometrium_validation_score_4"]
            + row["immune_support_score_4"]
            + row["cross_bulk_direction_consistency_score_2"]
        )
        row["bulk_expression_support_score_20"] = score
        if score >= 14:
            priority = "high_bulk_support"
        elif score >= 8:
            priority = "moderate_bulk_support"
        else:
            priority = "limited_bulk_support"
        row["bulk_support_class"] = priority
        rows.append(row)

    out = pd.DataFrame(rows)
    sort_cols = ["bulk_expression_support_score_20", "lesion_support_score_6", "independent_endometrium_validation_score_4", "immune_support_score_4"]
    out = out.sort_values(sort_cols, ascending=False)
    out.to_csv(BULK / "bulk_candidate_expression_support_scores.tsv", sep="\t", index=False)

    class_counts = out["bulk_support_class"].value_counts().to_dict()
    top_cols = [
        "gene_symbol",
        "gene_id",
        "genetic_priority",
        "bulk_expression_support_score_20",
        "bulk_support_class",
        "lesion_top_contrast",
        "lesion_logFC",
        "lesion_adj_p",
        "immune_top_contrast",
        "immune_logFC",
        "immune_adj_p",
        "cycle_driven_flag",
    ]
    lines = [
        "# Bulk expression support score summary",
        "",
        "This is an intermediate evidence layer, not the final therapeutic target score.",
        "",
        f"Candidate rows scored: {len(out)}",
        f"High bulk support: {class_counts.get('high_bulk_support', 0)}",
        f"Moderate bulk support: {class_counts.get('moderate_bulk_support', 0)}",
        f"Limited bulk support: {class_counts.get('limited_bulk_support', 0)}",
        "",
        "## Top candidates by bulk evidence",
        "",
        dataframe_to_markdown(out[top_cols].head(20)),
        "",
        "## Scoring notes",
        "",
        "- 6 points: GSE141549 lesion/tissue support.",
        "- 4 points: low cycle-confounding risk from GSE234354; strongly cycle-dependent genes receive 0 here.",
        "- 4 points: independent GSE51981 endometrium validation, with single-gene-probe sensitivity contributing only when direction is consistent.",
        "- 4 points: GSE313775 Th-subset immune support.",
        "- 2 points: direction consistency between GSE141549 eutopic endometrium and GSE51981 endometrium.",
        "",
    ]
    (BULK / "bulk_candidate_expression_support_score_summary.md").write_text("\n".join(lines), encoding="utf-8")
    self_review = [
        "# Phase 5 bulk expression support score self-review",
        "",
        "Verdict: PASS_WITH_CONDITIONS",
        "",
        "Checks passed:",
        "",
        "- Combined prespecified bulk model outputs into a candidate-level 20-point expression support layer.",
        "- Kept cycle dependence as a penalty/control component rather than interpreting it as disease evidence.",
        "- Preserved separate lesion, independent endometrium, immune and cross-bulk consistency fields for auditability.",
        "- Labeled this output as an intermediate evidence layer rather than a therapeutic target shortlist.",
        "",
        "Conditions before final scoring:",
        "",
        "- Bulk evidence must be integrated with single-cell/spatial localization before prioritizing targets.",
        "- Some candidates are strongly cycle-dependent; their bulk support must be interpreted cautiously unless single-cell/spatial lesion localization is strong.",
        "- Tissue-code assumptions for GSE141549 and donor-order assumptions for GSE313775 remain documented caveats.",
        "",
    ]
    (BULK / "phase5_bulk_expression_support_score_self_review.md").write_text("\n".join(self_review), encoding="utf-8")
    print(f"Wrote {BULK / 'bulk_candidate_expression_support_scores.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
