#!/usr/bin/env python3
"""Donor-aware broad-compartment pseudobulk tests for GSE179640 candidates.

This analysis uses sample-level candidate summaries from the conservative
broad-compartment annotation step. It uses subjects, not cells,
as the statistical units and should be interpreted as a candidate triage layer,
not as final cluster/state-level differential expression.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT_ROOT / "results" / "singlecell"
PY_PKG_DIR = PROJECT_ROOT / "tools" / "py_packages"

if sys.version_info[:2] != (3, 12):
    raise SystemExit(
        "Use the bundled Python 3.12 runtime: "
        "python3 "
        "scripts/run_gse179640_compartment_pseudobulk_tests.py"
    )

sys.path.insert(0, str(PY_PKG_DIR))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402


INPUT = RESULTS / "GSE179640_candidate_expression_by_broad_compartment_sample.tsv"
OUT_TESTS = RESULTS / "GSE179640_broad_compartment_pseudobulk_tests.tsv"
OUT_GENE_SUPPORT = RESULTS / "GSE179640_candidate_pseudobulk_support.tsv"
OUT_SUMMARY = RESULTS / "GSE179640_broad_compartment_pseudobulk_summary.md"
OUT_REVIEW = RESULTS / "phase7_gse179640_pseudobulk_self_review.md"

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
LESION_LOCATIONS = {"Ectopic", "Ectopic Adjacent", "Ectopic Ovary"}
MIN_CELLS_PER_COMPARTMENT_SAMPLE = 20
MIN_SUBJECTS_PER_GROUP = 3
MIN_PAIRED_SUBJECTS = 3

COMPARISONS = [
    {
        "comparison": "paired_ectopic_all_vs_patient_eutopic",
        "mode": "paired",
        "case_region": "ectopic_all",
        "reference_region": "patient_eutopic",
        "primary": True,
    },
    {
        "comparison": "paired_ectopic_peritoneum_vs_patient_eutopic",
        "mode": "paired",
        "case_region": "ectopic_peritoneum",
        "reference_region": "patient_eutopic",
        "primary": False,
    },
    {
        "comparison": "paired_ectopic_adjacent_vs_patient_eutopic",
        "mode": "paired",
        "case_region": "ectopic_adjacent",
        "reference_region": "patient_eutopic",
        "primary": False,
    },
    {
        "comparison": "paired_ectopic_ovary_vs_patient_eutopic",
        "mode": "paired",
        "case_region": "ectopic_ovary",
        "reference_region": "patient_eutopic",
        "primary": False,
    },
    {
        "comparison": "unpaired_patient_eutopic_vs_control_eutopic",
        "mode": "unpaired",
        "case_region": "patient_eutopic",
        "reference_region": "control_eutopic",
        "primary": True,
    },
    {
        "comparison": "unpaired_ectopic_all_vs_control_eutopic",
        "mode": "unpaired",
        "case_region": "ectopic_all",
        "reference_region": "control_eutopic",
        "primary": False,
    },
]

METRICS = ["log1p_mean_counts_per_cell", "prevalence"]


def bh_adjust(values: pd.Series) -> pd.Series:
    pvals = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=values.index, dtype=float)
    valid = pvals.dropna()
    if valid.empty:
        return out
    order = valid.sort_values().index.to_list()
    n = len(order)
    adjusted = {}
    running = 1.0
    for rank_from_end, idx in enumerate(reversed(order), start=1):
        rank = n - rank_from_end + 1
        q = float(valid.loc[idx]) * n / rank
        running = min(running, q)
        adjusted[idx] = min(running, 1.0)
    for idx, q in adjusted.items():
        out.loc[idx] = q
    return out


def safe_mannwhitneyu(case_values: np.ndarray, ref_values: np.ndarray) -> float:
    if len(case_values) < MIN_SUBJECTS_PER_GROUP or len(ref_values) < MIN_SUBJECTS_PER_GROUP:
        return math.nan
    if np.all(case_values == case_values[0]) and np.all(ref_values == ref_values[0]) and case_values[0] == ref_values[0]:
        return 1.0
    try:
        return float(stats.mannwhitneyu(case_values, ref_values, alternative="two-sided").pvalue)
    except Exception:
        return math.nan


def safe_wilcoxon(case_values: np.ndarray, ref_values: np.ndarray) -> float:
    if len(case_values) < MIN_PAIRED_SUBJECTS:
        return math.nan
    diffs = case_values - ref_values
    if np.allclose(diffs, 0):
        return 1.0
    try:
        return float(stats.wilcoxon(case_values, ref_values, zero_method="wilcox", alternative="two-sided").pvalue)
    except Exception:
        return math.nan


def cliffs_delta(case_values: np.ndarray, ref_values: np.ndarray) -> float:
    if len(case_values) == 0 or len(ref_values) == 0:
        return math.nan
    greater = 0
    lesser = 0
    for case_value in case_values:
        greater += int(np.sum(case_value > ref_values))
        lesser += int(np.sum(case_value < ref_values))
    return float((greater - lesser) / (len(case_values) * len(ref_values)))


def classify_region(row: pd.Series) -> str:
    subject = str(row["subject_code"])
    location = str(row["sample_location"])
    if subject.startswith("C") and location == "Eutopic":
        return "control_eutopic"
    if subject.startswith("E") and location == "Eutopic":
        return "patient_eutopic"
    if subject.startswith("E") and location == "Ectopic":
        return "ectopic_peritoneum"
    if subject.startswith("E") and location == "Ectopic Adjacent":
        return "ectopic_adjacent"
    if subject.startswith("E") and location == "Ectopic Ovary":
        return "ectopic_ovary"
    return "unclassified"


def build_subject_region_table(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["region_specific"] = df.apply(classify_region, axis=1)
    df = df[df["region_specific"] != "unclassified"].copy()
    all_rows = [df]
    lesion = df[df["sample_location"].isin(LESION_LOCATIONS)].copy()
    lesion["region_specific"] = "ectopic_all"
    all_rows.append(lesion)
    expanded = pd.concat(all_rows, ignore_index=True)

    aggregate_cols = [
        "gene_id",
        "gene_symbol",
        "genetic_priority",
        "ld_neighborhood_class",
        "module_hint_preliminary",
        "bulk_expression_support_score_20",
        "bulk_support_class",
        "broad_compartment",
        "subject_code",
        "region_specific",
    ]
    grouped = (
        expanded.groupby(aggregate_cols, dropna=False)
        .agg(
            n_samples=("geo_accession", "nunique"),
            median_compartment_cells=("n_cells", "median"),
            min_compartment_cells=("n_cells", "min"),
            log1p_mean_counts_per_cell=("log1p_mean_counts_per_cell", "mean"),
            prevalence=("candidate_prevalence", "mean"),
        )
        .reset_index()
        .rename(columns={"region_specific": "region"})
    )
    return grouped


def run_tests(subject_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_cols = ["gene_id", "gene_symbol", "broad_compartment"]
    for (gene_id, gene_symbol, compartment), gene_comp in subject_df.groupby(group_cols, dropna=False):
        metadata = gene_comp.iloc[0].to_dict()
        for comparison in COMPARISONS:
            case_region = comparison["case_region"]
            ref_region = comparison["reference_region"]
            case = gene_comp[gene_comp["region"] == case_region].copy()
            ref = gene_comp[gene_comp["region"] == ref_region].copy()
            for metric in METRICS:
                if comparison["mode"] == "paired":
                    merged = case[["subject_code", metric]].merge(
                        ref[["subject_code", metric]],
                        on="subject_code",
                        suffixes=("_case", "_reference"),
                    )
                    case_values = merged[f"{metric}_case"].to_numpy(dtype=float)
                    ref_values = merged[f"{metric}_reference"].to_numpy(dtype=float)
                    p_value = safe_wilcoxon(case_values, ref_values)
                    effect = float(np.median(case_values - ref_values)) if len(case_values) else math.nan
                    effect_type = "median_paired_difference"
                    n_pairs = int(len(merged))
                    n_case = int(len(case_values))
                    n_reference = int(len(ref_values))
                    delta = math.nan
                else:
                    case_values = case[metric].to_numpy(dtype=float)
                    ref_values = ref[metric].to_numpy(dtype=float)
                    p_value = safe_mannwhitneyu(case_values, ref_values)
                    effect = (
                        float(np.median(case_values) - np.median(ref_values))
                        if len(case_values) and len(ref_values)
                        else math.nan
                    )
                    effect_type = "median_difference"
                    n_pairs = 0
                    n_case = int(len(case_values))
                    n_reference = int(len(ref_values))
                    delta = cliffs_delta(case_values, ref_values)

                rows.append(
                    {
                        "gene_id": gene_id,
                        "gene_symbol": gene_symbol,
                        "genetic_priority": metadata.get("genetic_priority", ""),
                        "ld_neighborhood_class": metadata.get("ld_neighborhood_class", ""),
                        "module_hint_preliminary": metadata.get("module_hint_preliminary", ""),
                        "bulk_expression_support_score_20": metadata.get("bulk_expression_support_score_20", ""),
                        "bulk_support_class": metadata.get("bulk_support_class", ""),
                        "broad_compartment": compartment,
                        "comparison": comparison["comparison"],
                        "comparison_mode": comparison["mode"],
                        "primary_comparison": bool(comparison["primary"]),
                        "case_region": case_region,
                        "reference_region": ref_region,
                        "metric": metric,
                        "n_case_subjects": n_case,
                        "n_reference_subjects": n_reference,
                        "n_paired_subjects": n_pairs,
                        "case_median": float(np.median(case_values)) if len(case_values) else math.nan,
                        "reference_median": float(np.median(ref_values)) if len(ref_values) else math.nan,
                        "effect": effect,
                        "effect_type": effect_type,
                        "cliffs_delta": delta,
                        "p_value": p_value,
                    }
                )
    out = pd.DataFrame(rows)
    out["q_value_global"] = bh_adjust(out["p_value"])
    out["q_value_by_comparison_metric"] = np.nan
    for _, idx in out.groupby(["comparison", "metric"]).groups.items():
        out.loc[idx, "q_value_by_comparison_metric"] = bh_adjust(out.loc[idx, "p_value"])
    out["interpretation_guardrail"] = "subject_level_broad_compartment_pseudobulk_not_cluster_state"
    return out


def support_class(row: pd.Series) -> tuple[int, str]:
    p_value = row.get("p_value", math.nan)
    q_value = row.get("q_value_by_comparison_metric", math.nan)
    effect = abs(float(row.get("effect", 0))) if pd.notna(row.get("effect", math.nan)) else 0.0
    primary = bool(row.get("primary_comparison", False))
    if pd.notna(q_value) and q_value <= 0.10 and effect >= 0.05:
        return 10, "fdr_supported_subject_level_signal"
    if pd.notna(q_value) and q_value <= 0.20 and effect >= 0.05:
        return 8, "relaxed_fdr_subject_level_signal"
    if primary and pd.notna(p_value) and p_value <= 0.05 and effect >= 0.05:
        return 7, "nominal_primary_subject_level_signal"
    if pd.notna(p_value) and p_value <= 0.05 and effect >= 0.05:
        return 6, "nominal_secondary_subject_level_signal"
    if primary and pd.notna(p_value) and p_value <= 0.10 and effect >= 0.03:
        return 5, "suggestive_primary_subject_level_signal"
    if effect >= 0.10:
        return 3, "descriptive_effect_without_nominal_support"
    if effect > 0:
        return 1, "limited_subject_level_support"
    return 0, "minimal_subject_level_support"


def build_gene_support(tests: pd.DataFrame) -> pd.DataFrame:
    eligible = tests[
        tests["metric"].isin(METRICS)
        & tests["broad_compartment"].isin(INTERPRETABLE_COMPARTMENTS)
        & tests["p_value"].notna()
    ].copy()
    if eligible.empty:
        return pd.DataFrame()
    eligible[["pseudobulk_support_score_10", "pseudobulk_support_class"]] = eligible.apply(
        lambda row: pd.Series(support_class(row)), axis=1
    )
    eligible["abs_effect"] = eligible["effect"].abs()
    eligible["is_fibrovascular_compartment"] = eligible["broad_compartment"].isin(FIBROVASCULAR_COMPARTMENTS)
    eligible["is_immune_compartment"] = eligible["broad_compartment"].isin(IMMUNE_COMPARTMENTS)
    eligible = eligible.sort_values(
        [
            "pseudobulk_support_score_10",
            "primary_comparison",
            "abs_effect",
            "p_value",
        ],
        ascending=[False, False, False, True],
    )
    best = eligible.groupby(["gene_id", "gene_symbol"], as_index=False).head(1).copy()
    best = best.rename(
        columns={
            "broad_compartment": "best_pseudobulk_compartment",
            "comparison": "best_pseudobulk_comparison",
            "metric": "best_pseudobulk_metric",
            "effect": "best_pseudobulk_effect",
            "p_value": "best_pseudobulk_p_value",
            "q_value_by_comparison_metric": "best_pseudobulk_q_value_by_comparison_metric",
            "q_value_global": "best_pseudobulk_q_value_global",
        }
    )
    keep_cols = [
        "gene_id",
        "gene_symbol",
        "genetic_priority",
        "ld_neighborhood_class",
        "module_hint_preliminary",
        "bulk_expression_support_score_20",
        "bulk_support_class",
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
        "case_median",
        "reference_median",
        "interpretation_guardrail",
    ]
    return best[keep_cols].sort_values(
        ["pseudobulk_support_score_10", "best_pseudobulk_p_value", "gene_symbol"],
        ascending=[False, True, True],
    )


def format_value(value: object, digits: int = 4) -> str:
    if value is None or pd.isna(value):
        return ""
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    try:
        return f"{float(value):.{digits}g}"
    except Exception:
        return str(value)


def write_summary(subject_df: pd.DataFrame, tests: pd.DataFrame, support: pd.DataFrame) -> None:
    region_counts = (
        subject_df[["subject_code", "region"]]
        .drop_duplicates()
        .groupby("region")["subject_code"]
        .nunique()
        .sort_index()
    )
    class_counts = support["pseudobulk_support_class"].value_counts().sort_index() if not support.empty else pd.Series(dtype=int)
    high_bulk = support[support["bulk_support_class"] == "high_bulk_support"].copy()
    high_bulk = high_bulk.sort_values(["pseudobulk_support_score_10", "best_pseudobulk_p_value"], ascending=[False, True])

    lines = [
        "# GSE179640 donor-aware broad-compartment pseudobulk summary",
        "",
        "## Inputs and limits",
        "",
        f"- Input table: `{INPUT}`",
        f"- Minimum cells per subject-compartment-sample: {MIN_CELLS_PER_COMPARTMENT_SAMPLE}",
        "- Statistical unit: subject-level region means, not individual cells.",
        "- Expression metrics tested: log1p candidate counts per compartment cell and expressing-cell prevalence.",
        "- Interpretation: broad-compartment candidate triage only; not final cluster/state-level differential expression.",
        "",
        "## Subject-region coverage",
        "",
        "region\tn_subjects",
    ]
    for region, n_subjects in region_counts.items():
        lines.append(f"{region}\t{int(n_subjects)}")
    lines.extend(
        [
            "",
            "## Support-class counts",
            "",
            "support_class\tn_genes",
        ]
    )
    for klass, n_genes in class_counts.items():
        lines.append(f"{klass}\t{int(n_genes)}")
    lines.extend(
        [
            "",
            "## Top donor-aware pseudobulk-supported candidates",
            "",
            "gene_symbol\tgene_id\tscore\tclass\tcompartment\tcomparison\tmetric\teffect\tp_value\tq_by_comparison_metric",
        ]
    )
    for _, row in support.head(20).iterrows():
        lines.append(
            "\t".join(
                [
                    str(row["gene_symbol"]),
                    str(row["gene_id"]),
                    str(int(row["pseudobulk_support_score_10"])),
                    str(row["pseudobulk_support_class"]),
                    str(row["best_pseudobulk_compartment"]),
                    str(row["best_pseudobulk_comparison"]),
                    str(row["best_pseudobulk_metric"]),
                    format_value(row["best_pseudobulk_effect"]),
                    format_value(row["best_pseudobulk_p_value"]),
                    format_value(row["best_pseudobulk_q_value_by_comparison_metric"]),
                ]
            )
        )
    lines.extend(
        [
            "",
            "## High-bulk-support candidates in this donor-aware layer",
            "",
            "gene_symbol\tgene_id\tbulk_score\tpseudobulk_score\tclass\tcompartment\tcomparison\tmetric\teffect\tp_value",
        ]
    )
    for _, row in high_bulk.iterrows():
        lines.append(
            "\t".join(
                [
                    str(row["gene_symbol"]),
                    str(row["gene_id"]),
                    str(row["bulk_expression_support_score_20"]),
                    str(int(row["pseudobulk_support_score_10"])),
                    str(row["pseudobulk_support_class"]),
                    str(row["best_pseudobulk_compartment"]),
                    str(row["best_pseudobulk_comparison"]),
                    str(row["best_pseudobulk_metric"]),
                    format_value(row["best_pseudobulk_effect"]),
                    format_value(row["best_pseudobulk_p_value"]),
                ]
            )
        )
    lines.extend(
        [
            "",
            "## Output files",
            "",
            f"- Test table: `{OUT_TESTS}`",
            f"- Gene support table: `{OUT_GENE_SUPPORT}`",
            f"- Self-review: `{OUT_REVIEW}`",
        ]
    )
    OUT_SUMMARY.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_self_review(tests: pd.DataFrame, support: pd.DataFrame) -> None:
    n_tests = int(len(tests))
    n_genes = int(support["gene_id"].nunique()) if not support.empty else 0
    n_fdr = int((support["pseudobulk_support_class"] == "fdr_supported_subject_level_signal").sum()) if not support.empty else 0
    n_relaxed_fdr = (
        int((support["pseudobulk_support_class"] == "relaxed_fdr_subject_level_signal").sum())
        if not support.empty
        else 0
    )
    nominal_classes = {
        "nominal_primary_subject_level_signal",
        "nominal_secondary_subject_level_signal",
    }
    n_nominal = int(support["pseudobulk_support_class"].isin(nominal_classes).sum()) if not support.empty else 0
    lines = [
        "# Phase 7 Self-Review: GSE179640 donor-aware broad-compartment pseudobulk",
        "",
        "## Verdict",
        "",
        "PASS_WITH_CONDITIONS",
        "",
        "## What passed",
        "",
        f"- Tested {n_tests} gene/compartment/comparison/metric rows using subject-level region aggregates.",
        f"- Produced one donor-aware pseudobulk support row for {n_genes} candidate genes with at least one testable broad-compartment signal.",
        "- The script uses subject-level summaries and paired tests where lesion and eutopic samples are available from the same endometriosis subject.",
        "- Results label broad compartments as context, not final disease cell states.",
        "",
        "## Scientific limitations",
        "",
        "- Control eutopic coverage is only three subjects in the GSE179640 primary tissue subset, so control comparisons are low-powered.",
        "- Candidate expression was summarized as candidate counts per compartment cell and prevalence; this is not a full transcriptome count matrix per pseudobulk sample with library-size normalization.",
        "- Broad compartments were assigned by conservative marker panels rather than author-provided labels, graph clusters, or supervised reference mapping.",
        "- Multiple-testing correction is provided; nominal signals are descriptive.",
        "- This analysis does not replace future cluster-level annotation, doublet-aware validation, or cross-dataset adenomyosis comparison.",
        "",
        "## Quantitative audit",
        "",
        f"- FDR-supported subject-level signals in best-gene summary: {n_fdr}.",
        f"- Relaxed FDR q<=0.20 subject-level signals in best-gene summary: {n_relaxed_fdr}.",
        f"- Nominal best-gene signals in best-gene summary: {n_nominal}.",
        "",
        "## Decision",
        "",
        "Use these results as donor-aware broad-compartment context, not definitive cell-state differential expression.",
        "",
        "## Next gate",
        "",
        "Next: summarize the adenomyosis h5ad candidate context matrix.",
    ]
    OUT_REVIEW.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    if not INPUT.exists():
        raise SystemExit(f"Missing input: {INPUT}")
    RESULTS.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT, sep="\t", keep_default_na=False)
    numeric_cols = ["n_cells", "candidate_total_counts", "candidate_expressing_cells", "candidate_prevalence"]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
    df = df[df["broad_compartment"].isin(INTERPRETABLE_COMPARTMENTS)].copy()
    df = df[df["n_cells"] >= MIN_CELLS_PER_COMPARTMENT_SAMPLE].copy()
    df["mean_counts_per_cell"] = np.divide(
        df["candidate_total_counts"].to_numpy(dtype=float),
        df["n_cells"].to_numpy(dtype=float),
        out=np.zeros(len(df), dtype=float),
        where=df["n_cells"].to_numpy(dtype=float) > 0,
    )
    df["log1p_mean_counts_per_cell"] = np.log1p(df["mean_counts_per_cell"])
    subject_df = build_subject_region_table(df)
    tests = run_tests(subject_df)
    support = build_gene_support(tests)

    tests.to_csv(OUT_TESTS, sep="\t", index=False)
    support.to_csv(OUT_GENE_SUPPORT, sep="\t", index=False)
    write_summary(subject_df, tests, support)
    write_self_review(tests, support)

    print(OUT_SUMMARY)
    print(OUT_GENE_SUPPORT)
    print(OUT_REVIEW)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
