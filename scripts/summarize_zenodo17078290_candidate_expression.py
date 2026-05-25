#!/usr/bin/env python3
"""Summarize GWAS candidate localization in the Zenodo 17078290 adenomyosis h5ad."""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW = PROJECT_ROOT / "data" / "raw_downloads"
RESULTS = PROJECT_ROOT / "results" / "singlecell"
PY_PKG_DIR = PROJECT_ROOT / "tools" / "py_packages"

if sys.version_info[:2] != (3, 12):
    raise SystemExit(
        "Use the bundled Python 3.12 runtime: "
        "python3 "
        "scripts/summarize_zenodo17078290_candidate_expression.py"
    )

sys.path.insert(0, str(PY_PKG_DIR))

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import sparse, stats  # noqa: E402


H5AD = RAW / "Zenodo_17078290__celllable.diff_PRO.h5ad"
GWAS = PROJECT_ROOT / "results" / "gwas" / "gwas_candidate_gene_universe.tsv"
BULK = PROJECT_ROOT / "results" / "bulk" / "bulk_candidate_expression_support_scores.tsv"
COVERAGE = RESULTS / "Zenodo17078290_candidate_gene_coverage.tsv"

OUT_BY_SAMPLE_CLUSTER = RESULTS / "Zenodo17078290_candidate_expression_by_sample_cluster.tsv"
OUT_SUMMARY_GROUP_CLUSTER = RESULTS / "Zenodo17078290_candidate_expression_summary_by_group_cluster.tsv"
OUT_TESTS = RESULTS / "Zenodo17078290_candidate_group_pseudobulk_tests.tsv"
OUT_MATRIX = RESULTS / "Zenodo17078290_candidate_localization_matrix.tsv"
OUT_SUMMARY = RESULTS / "Zenodo17078290_candidate_expression_summary.md"
OUT_REVIEW = RESULTS / "phase9_zenodo17078290_candidate_expression_self_review.md"

FIBROVASCULAR_CLUSTERS = {"Fibroblasts", "Endothelial cells", "Mural cells"}
IMMUNE_CLUSTERS = {"Mononuclear phagocytes", "T and NK cells", "B cells", "Mast cells"}
EPITHELIAL_CLUSTERS = {"Epithelial cells"}
ADENOMYOSIS_GROUPS = {"AM_EcM", "AM_EuM"}
MIN_CELLS_PER_SAMPLE_CLUSTER = 20

COMPARISONS = [
    {
        "comparison": "paired_AM_EcM_vs_AM_EuM",
        "mode": "paired",
        "case_group": "AM_EcM",
        "reference_group": "AM_EuM",
    },
    {
        "comparison": "unpaired_AM_EcM_vs_C_EM",
        "mode": "unpaired",
        "case_group": "AM_EcM",
        "reference_group": "C_EM",
    },
    {
        "comparison": "unpaired_AM_EuM_vs_C_EM",
        "mode": "unpaired",
        "case_group": "AM_EuM",
        "reference_group": "C_EM",
    },
]
METRICS = ["log1p_mean_counts_per_cell", "candidate_prevalence"]


def clean_id(value: object) -> str:
    return str(value).split(".")[0]


def sample_index(sample: str) -> str:
    match = re.search(r"(\d+)$", sample)
    return match.group(1) if match else sample


def bh_adjust(values: pd.Series) -> pd.Series:
    pvals = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=values.index, dtype=float)
    valid = pvals.dropna()
    if valid.empty:
        return out
    order = valid.sort_values().index.to_list()
    n = len(order)
    running = 1.0
    adjusted = {}
    for rank_from_end, idx in enumerate(reversed(order), start=1):
        rank = n - rank_from_end + 1
        q = float(valid.loc[idx]) * n / rank
        running = min(running, q)
        adjusted[idx] = min(running, 1.0)
    for idx, q in adjusted.items():
        out.loc[idx] = q
    return out


def safe_mannwhitneyu(case_values: np.ndarray, ref_values: np.ndarray) -> float:
    if len(case_values) < 3 or len(ref_values) < 3:
        return math.nan
    if np.all(case_values == case_values[0]) and np.all(ref_values == ref_values[0]) and case_values[0] == ref_values[0]:
        return 1.0
    try:
        return float(stats.mannwhitneyu(case_values, ref_values, alternative="two-sided").pvalue)
    except Exception:
        return math.nan


def safe_wilcoxon(case_values: np.ndarray, ref_values: np.ndarray) -> float:
    if len(case_values) < 3:
        return math.nan
    diffs = case_values - ref_values
    if np.allclose(diffs, 0):
        return 1.0
    try:
        return float(stats.wilcoxon(case_values, ref_values, zero_method="wilcox", alternative="two-sided").pvalue)
    except Exception:
        return math.nan


def load_candidate_metadata() -> pd.DataFrame:
    gwas = pd.read_csv(GWAS, sep="\t", keep_default_na=False)
    bulk = pd.read_csv(BULK, sep="\t", keep_default_na=False)
    coverage = pd.read_csv(COVERAGE, sep="\t", keep_default_na=False)
    gwas["gene_id"] = gwas["gene_id"].map(clean_id)
    bulk["gene_id"] = bulk["gene_id"].map(clean_id)
    coverage["gene_id"] = coverage["gene_id"].map(clean_id)
    merged = gwas.merge(
        bulk[["gene_id", "bulk_expression_support_score_20", "bulk_support_class"]],
        on="gene_id",
        how="left",
    ).merge(
        coverage[["gene_id", "matched", "match_basis"]],
        on="gene_id",
        how="left",
    )
    merged["matched"] = merged["matched"].astype(str).str.lower().eq("true")
    return merged


def matrix_subset(raw_matrix, gene_indices: list[int]):
    subset = raw_matrix[:, gene_indices]
    if sparse.issparse(subset):
        return subset.tocsr()
    return sparse.csr_matrix(subset)


def summarize_by_sample_cluster(adata: ad.AnnData, candidates: pd.DataFrame) -> pd.DataFrame:
    matched = candidates[candidates["matched"] & candidates["gene_symbol"].astype(bool)].copy()
    var_to_idx = {str(name): idx for idx, name in enumerate(adata.var_names)}
    matched = matched[matched["gene_symbol"].isin(var_to_idx)].copy()
    matched = matched.drop_duplicates("gene_id")
    gene_symbols = matched["gene_symbol"].tolist()
    gene_indices = [var_to_idx[symbol] for symbol in gene_symbols]
    raw_matrix = adata.layers["raw"] if "raw" in adata.layers else adata.X
    mat = matrix_subset(raw_matrix, gene_indices)
    obs = adata.obs[["sample", "gname", "cluster_standard"]].copy()
    obs["sample"] = obs["sample"].astype(str)
    obs["gname"] = obs["gname"].astype(str)
    obs["cluster_standard"] = obs["cluster_standard"].astype(str)
    obs["sample_index"] = obs["sample"].map(sample_index)

    candidate_records = matched.to_dict("records")
    rows: list[dict[str, object]] = []
    for (sample, gname, cluster, idx), group in obs.groupby(["sample", "gname", "cluster_standard", "sample_index"], observed=True):
        cell_idx = group.index
        positional = adata.obs.index.get_indexer(cell_idx)
        sub = mat[positional, :]
        n_cells = int(sub.shape[0])
        total_counts = np.asarray(sub.sum(axis=0)).ravel().astype(float)
        expressing_cells = np.asarray((sub > 0).sum(axis=0)).ravel().astype(int)
        mean_counts = np.divide(total_counts, n_cells, out=np.zeros_like(total_counts), where=n_cells > 0)
        prevalence = np.divide(expressing_cells, n_cells, out=np.zeros_like(total_counts), where=n_cells > 0)
        for gene_idx, candidate in enumerate(candidate_records):
            rows.append(
                {
                    "gene_id": clean_id(candidate["gene_id"]),
                    "gene_symbol": candidate["gene_symbol"],
                    "genetic_priority": candidate.get("genetic_priority", ""),
                    "ld_neighborhood_class": candidate.get("ld_neighborhood_class", ""),
                    "module_hint_preliminary": candidate.get("module_hint_preliminary", ""),
                    "bulk_expression_support_score_20": candidate.get("bulk_expression_support_score_20", ""),
                    "bulk_support_class": candidate.get("bulk_support_class", ""),
                    "match_basis": candidate.get("match_basis", ""),
                    "sample": sample,
                    "sample_index": idx,
                    "gname": gname,
                    "cluster_standard": cluster,
                    "n_cells": n_cells,
                    "candidate_total_counts": total_counts[gene_idx],
                    "candidate_expressing_cells": int(expressing_cells[gene_idx]),
                    "candidate_prevalence": prevalence[gene_idx],
                    "mean_counts_per_cell": mean_counts[gene_idx],
                    "log1p_mean_counts_per_cell": float(np.log1p(mean_counts[gene_idx])),
                }
            )
    return pd.DataFrame(rows)


def summarize_group_cluster(sample_cluster: pd.DataFrame) -> pd.DataFrame:
    eligible = sample_cluster[sample_cluster["n_cells"] >= MIN_CELLS_PER_SAMPLE_CLUSTER].copy()
    summary = (
        eligible.groupby(
            [
                "gene_id",
                "gene_symbol",
                "genetic_priority",
                "ld_neighborhood_class",
                "module_hint_preliminary",
                "bulk_expression_support_score_20",
                "bulk_support_class",
                "gname",
                "cluster_standard",
            ],
            dropna=False,
        )
        .agg(
            n_samples_with_cluster=("sample", "nunique"),
            median_cluster_cells=("n_cells", "median"),
            min_cluster_cells=("n_cells", "min"),
            median_prevalence=("candidate_prevalence", "median"),
            mean_prevalence=("candidate_prevalence", "mean"),
            sample_detection_fraction=("candidate_prevalence", lambda x: float((x > 0).mean())),
            median_log1p_mean_counts=("log1p_mean_counts_per_cell", "median"),
        )
        .reset_index()
    )
    return summary


def run_group_tests(sample_cluster: pd.DataFrame) -> pd.DataFrame:
    eligible = sample_cluster[sample_cluster["n_cells"] >= MIN_CELLS_PER_SAMPLE_CLUSTER].copy()
    rows: list[dict[str, object]] = []
    for (gene_id, gene_symbol, cluster), sub in eligible.groupby(["gene_id", "gene_symbol", "cluster_standard"], dropna=False):
        first = sub.iloc[0]
        for comparison in COMPARISONS:
            case = sub[sub["gname"] == comparison["case_group"]].copy()
            ref = sub[sub["gname"] == comparison["reference_group"]].copy()
            for metric in METRICS:
                if comparison["mode"] == "paired":
                    merged = case[["sample_index", metric]].merge(
                        ref[["sample_index", metric]],
                        on="sample_index",
                        suffixes=("_case", "_reference"),
                    )
                    case_values = merged[f"{metric}_case"].to_numpy(dtype=float)
                    ref_values = merged[f"{metric}_reference"].to_numpy(dtype=float)
                    p_value = safe_wilcoxon(case_values, ref_values)
                    effect = float(np.median(case_values - ref_values)) if len(case_values) else math.nan
                    n_pairs = int(len(merged))
                    n_case = int(len(case_values))
                    n_ref = int(len(ref_values))
                else:
                    case_values = case[metric].to_numpy(dtype=float)
                    ref_values = ref[metric].to_numpy(dtype=float)
                    p_value = safe_mannwhitneyu(case_values, ref_values)
                    effect = (
                        float(np.median(case_values) - np.median(ref_values))
                        if len(case_values) and len(ref_values)
                        else math.nan
                    )
                    n_pairs = 0
                    n_case = int(len(case_values))
                    n_ref = int(len(ref_values))
                rows.append(
                    {
                        "gene_id": gene_id,
                        "gene_symbol": gene_symbol,
                        "genetic_priority": first.get("genetic_priority", ""),
                        "ld_neighborhood_class": first.get("ld_neighborhood_class", ""),
                        "module_hint_preliminary": first.get("module_hint_preliminary", ""),
                        "bulk_expression_support_score_20": first.get("bulk_expression_support_score_20", ""),
                        "bulk_support_class": first.get("bulk_support_class", ""),
                        "cluster_standard": cluster,
                        "comparison": comparison["comparison"],
                        "comparison_mode": comparison["mode"],
                        "case_group": comparison["case_group"],
                        "reference_group": comparison["reference_group"],
                        "metric": metric,
                        "n_case_samples": n_case,
                        "n_reference_samples": n_ref,
                        "n_paired_samples": n_pairs,
                        "case_median": float(np.median(case_values)) if len(case_values) else math.nan,
                        "reference_median": float(np.median(ref_values)) if len(ref_values) else math.nan,
                        "effect": effect,
                        "p_value": p_value,
                    }
                )
    tests = pd.DataFrame(rows)
    tests["q_value_global"] = bh_adjust(tests["p_value"])
    tests["q_value_by_comparison_metric"] = np.nan
    for _, idx in tests.groupby(["comparison", "metric"]).groups.items():
        tests.loc[idx, "q_value_by_comparison_metric"] = bh_adjust(tests.loc[idx, "p_value"])
    tests["interpretation_guardrail"] = "n3_sample_level_screen_not_definitive_DE"
    return tests


def detectability_score(value: float) -> int:
    if value >= 0.20:
        return 5
    if value >= 0.10:
        return 4
    if value >= 0.05:
        return 3
    if value > 0:
        return 1
    return 0


def robustness_score(row: pd.Series | None) -> int:
    if row is None:
        return 0
    n_samples = int(row.get("n_samples_with_cluster", 0))
    median_cells = float(row.get("median_cluster_cells", 0))
    detection_fraction = float(row.get("sample_detection_fraction", 0))
    if n_samples >= 3 and median_cells >= 50 and detection_fraction >= 0.67:
        return 5
    if n_samples >= 2 and median_cells >= 20 and detection_fraction >= 0.50:
        return 3
    if n_samples >= 1 and detection_fraction > 0:
        return 1
    return 0


def group_spread_score(summary: pd.DataFrame) -> int:
    group_prev = {
        group: float(summary[summary["gname"] == group]["median_prevalence"].max())
        if not summary[summary["gname"] == group].empty
        else 0.0
        for group in ADENOMYOSIS_GROUPS
    }
    n_groups = sum(1 for value in group_prev.values() if value >= 0.05)
    if n_groups == 2:
        return 5
    if n_groups == 1:
        return 2
    return 0


def build_localization_matrix(candidates: pd.DataFrame, summary: pd.DataFrame, tests: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    summary = summary.copy()
    summary["median_prevalence"] = pd.to_numeric(summary["median_prevalence"], errors="coerce").fillna(0)
    tests = tests.copy()
    for _, candidate in candidates.iterrows():
        gene_id = clean_id(candidate["gene_id"])
        gene_summary = summary[summary["gene_id"] == gene_id].copy()
        matched = bool(candidate.get("matched", False))
        if gene_summary.empty:
            rows.append(
                {
                    "gene_id": gene_id,
                    "gene_symbol": candidate.get("gene_symbol", ""),
                    "genetic_priority": candidate.get("genetic_priority", ""),
                    "ld_neighborhood_class": candidate.get("ld_neighborhood_class", ""),
                    "module_hint_preliminary": candidate.get("module_hint_preliminary", ""),
                    "bulk_expression_support_score_20": candidate.get("bulk_expression_support_score_20", ""),
                    "bulk_support_class": candidate.get("bulk_support_class", ""),
                    "matched_in_h5ad": matched,
                    "zenodo17078290_localization_score_25": 0,
                    "zenodo17078290_localization_class": "not_detected_or_unmatched",
                    "interpretation_guardrail": "adenomyosis_h5ad_scRNA_not_spatial",
                }
            )
            continue

        adeno = gene_summary[gene_summary["gname"].isin(ADENOMYOSIS_GROUPS)].copy()
        adeno = adeno.sort_values(
            ["median_prevalence", "median_cluster_cells", "sample_detection_fraction"],
            ascending=[False, False, False],
        )
        top = adeno.iloc[0] if not adeno.empty else None
        max_adeno_prev = float(top["median_prevalence"]) if top is not None else 0.0
        max_control_prev = (
            float(gene_summary[gene_summary["gname"] == "C_EM"]["median_prevalence"].max())
            if not gene_summary[gene_summary["gname"] == "C_EM"].empty
            else 0.0
        )
        def max_prev_for(clusters: set[str]) -> tuple[float, str]:
            eligible = adeno[adeno["cluster_standard"].isin(clusters)].copy()
            if eligible.empty:
                return 0.0, ""
            eligible = eligible.sort_values("median_prevalence", ascending=False)
            row = eligible.iloc[0]
            return float(row["median_prevalence"]), f"{row['gname']}:{row['cluster_standard']}"

        fibro_prev, fibro_label = max_prev_for(FIBROVASCULAR_CLUSTERS)
        immune_prev, immune_label = max_prev_for(IMMUNE_CLUSTERS)
        epithelial_prev, epithelial_label = max_prev_for(EPITHELIAL_CLUSTERS)

        gene_tests = tests[(tests["gene_id"] == gene_id) & tests["p_value"].notna()].copy()
        if not gene_tests.empty:
            gene_tests["abs_effect"] = gene_tests["effect"].abs()
            gene_tests = gene_tests.sort_values(["abs_effect", "p_value"], ascending=[False, True])
            best_test = gene_tests.iloc[0]
            best_p = float(best_test["p_value"])
            best_effect = float(best_test["effect"])
            best_comparison = str(best_test["comparison"])
            best_metric = str(best_test["metric"])
            best_cluster_test = str(best_test["cluster_standard"])
        else:
            best_p = math.nan
            best_effect = math.nan
            best_comparison = ""
            best_metric = ""
            best_cluster_test = ""

        loc_score = detectability_score(max_adeno_prev)
        robust = robustness_score(top)
        fibro_score = detectability_score(fibro_prev)
        spread = group_spread_score(gene_summary)
        effect_bonus = 0
        if not math.isnan(best_p) and best_p <= 0.10 and abs(best_effect) >= 0.05:
            effect_bonus = 5
        elif not math.isnan(best_effect) and abs(best_effect) >= 0.10:
            effect_bonus = 3
        elif not math.isnan(best_effect) and abs(best_effect) > 0:
            effect_bonus = 1
        score = loc_score + robust + fibro_score + spread + effect_bonus
        if score >= 15:
            klass = "high_adenomyosis_scRNA_localization"
        elif score >= 10:
            klass = "moderate_adenomyosis_scRNA_localization"
        elif score >= 5:
            klass = "limited_adenomyosis_scRNA_localization"
        else:
            klass = "minimal_adenomyosis_scRNA_localization"
        rows.append(
            {
                "gene_id": gene_id,
                "gene_symbol": candidate.get("gene_symbol", ""),
                "genetic_priority": candidate.get("genetic_priority", ""),
                "ld_neighborhood_class": candidate.get("ld_neighborhood_class", ""),
                "module_hint_preliminary": candidate.get("module_hint_preliminary", ""),
                "bulk_expression_support_score_20": candidate.get("bulk_expression_support_score_20", ""),
                "bulk_support_class": candidate.get("bulk_support_class", ""),
                "matched_in_h5ad": matched,
                "top_adeno_group": top.get("gname", "") if top is not None else "",
                "top_adeno_cluster": top.get("cluster_standard", "") if top is not None else "",
                "top_adeno_median_prevalence": max_adeno_prev,
                "control_max_median_prevalence": max_control_prev,
                "fibrovascular_max_prevalence": fibro_prev,
                "fibrovascular_top_group_cluster": fibro_label,
                "immune_max_prevalence": immune_prev,
                "immune_top_group_cluster": immune_label,
                "epithelial_max_prevalence": epithelial_prev,
                "epithelial_top_group_cluster": epithelial_label,
                "best_sample_level_comparison": best_comparison,
                "best_sample_level_metric": best_metric,
                "best_sample_level_cluster": best_cluster_test,
                "best_sample_level_effect": best_effect,
                "best_sample_level_p_value": best_p,
                "zenodo17078290_location_score_5": loc_score,
                "zenodo17078290_sample_robustness_score_5": robust,
                "zenodo17078290_fibrovascular_score_5": fibro_score,
                "zenodo17078290_adeno_group_spread_score_5": spread,
                "zenodo17078290_effect_bonus_5": effect_bonus,
                "zenodo17078290_localization_score_25": score,
                "zenodo17078290_localization_class": klass,
                "interpretation_guardrail": "adenomyosis_h5ad_scRNA_not_spatial",
            }
        )
    out = pd.DataFrame(rows)
    out = out.sort_values(
        ["zenodo17078290_localization_score_25", "fibrovascular_max_prevalence", "top_adeno_median_prevalence", "gene_symbol"],
        ascending=[False, False, False, True],
    )
    return out


def write_summary(sample_cluster: pd.DataFrame, group_summary: pd.DataFrame, tests: pd.DataFrame, matrix: pd.DataFrame) -> None:
    class_counts = matrix["zenodo17078290_localization_class"].value_counts().to_dict()
    group_counts = sample_cluster[["sample", "gname"]].drop_duplicates()["gname"].value_counts().to_dict()
    cluster_counts = (
        sample_cluster[["sample", "gname", "cluster_standard", "n_cells"]]
        .drop_duplicates()
        .groupby(["gname", "cluster_standard"])["n_cells"]
        .median()
        .reset_index()
        .sort_values(["gname", "cluster_standard"])
    )
    high_bulk = matrix[matrix["bulk_support_class"] == "high_bulk_support"].copy()
    lines = [
        "# Zenodo 17078290 adenomyosis candidate expression summary",
        "",
        "## Dataset interpretation guardrail",
        "",
        "- The h5ad contains scRNA-seq cell labels and embeddings but no detected spatial coordinates in `obsm`/`uns`; use this as adenomyosis scRNA/cell-label evidence, not spatial evidence.",
        "- Encoded groups are retained as `C_EM`, `AM_EuM`, and `AM_EcM`; biological labels should be verified against the source paper before manuscript submission.",
        "- Sample-level comparisons have only three samples per encoded group and are descriptive/supportive rather than definitive DE.",
        "",
        "## Coverage",
        "",
        f"- Candidate expression summarized for {sample_cluster['gene_id'].nunique()} matched candidates.",
        f"- Total sample-cluster candidate rows: {len(sample_cluster)}.",
        f"- Localization matrix rows: {len(matrix)}.",
        f"- Localization class counts: `{class_counts}`",
        "- Localization score range: 0-25 from detectability, sample robustness, fibrovascular localization, adenomyosis encoded-group spread and descriptive effect-size bonus.",
        f"- Encoded group sample counts: `{group_counts}`",
        "",
        "## Median cells per encoded group and broad cell label",
        "",
        "gname\tcluster_standard\tmedian_cells",
    ]
    for _, row in cluster_counts.iterrows():
        lines.append(f"{row['gname']}\t{row['cluster_standard']}\t{row['n_cells']:.1f}")
    lines.extend(
        [
            "",
            "## Top adenomyosis scRNA-localized candidates",
            "",
            "gene_symbol\tgene_id\tscore\tclass\ttop_group\ttop_cluster\ttop_prev\tfibrovascular_prev\tbest_comparison\tbest_effect\tbest_p",
        ]
    )
    for _, row in matrix.head(25).iterrows():
        lines.append(
            "\t".join(
                [
                    str(row["gene_symbol"]),
                    str(row["gene_id"]),
                    str(row["zenodo17078290_localization_score_25"]),
                    str(row["zenodo17078290_localization_class"]),
                    str(row.get("top_adeno_group", "")),
                    str(row.get("top_adeno_cluster", "")),
                    f"{float(row.get('top_adeno_median_prevalence', 0)):.4f}",
                    f"{float(row.get('fibrovascular_max_prevalence', 0)):.4f}",
                    str(row.get("best_sample_level_comparison", "")),
                    f"{float(row.get('best_sample_level_effect', 0)):.4g}" if pd.notna(row.get("best_sample_level_effect", np.nan)) else "",
                    f"{float(row.get('best_sample_level_p_value', 0)):.4g}" if pd.notna(row.get("best_sample_level_p_value", np.nan)) else "",
                ]
            )
        )
    lines.extend(
        [
            "",
            "## High-bulk-support candidates in adenomyosis h5ad",
            "",
            "gene_symbol\tgene_id\tbulk_score\tadeno_score\tadeno_class\ttop_group\ttop_cluster\ttop_prev\tfibrovascular_prev",
        ]
    )
    for _, row in high_bulk.iterrows():
        lines.append(
            "\t".join(
                [
                    str(row["gene_symbol"]),
                    str(row["gene_id"]),
                    str(row["bulk_expression_support_score_20"]),
                    str(row["zenodo17078290_localization_score_25"]),
                    str(row["zenodo17078290_localization_class"]),
                    str(row.get("top_adeno_group", "")),
                    str(row.get("top_adeno_cluster", "")),
                    f"{float(row.get('top_adeno_median_prevalence', 0)):.4f}",
                    f"{float(row.get('fibrovascular_max_prevalence', 0)):.4f}",
                ]
            )
        )
    lines.extend(
        [
            "",
            "## Output files",
            "",
            f"- By sample/cluster: `{OUT_BY_SAMPLE_CLUSTER}`",
            f"- Group/cluster summary: `{OUT_SUMMARY_GROUP_CLUSTER}`",
            f"- Sample-level tests: `{OUT_TESTS}`",
            f"- Localization matrix: `{OUT_MATRIX}`",
            f"- Self-review: `{OUT_REVIEW}`",
        ]
    )
    OUT_SUMMARY.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_self_review(sample_cluster: pd.DataFrame, tests: pd.DataFrame, matrix: pd.DataFrame) -> None:
    n_fdr = int((pd.to_numeric(tests["q_value_by_comparison_metric"], errors="coerce") <= 0.10).sum())
    lines = [
        "# Phase 9 Self-Review: Zenodo 17078290 adenomyosis candidate expression",
        "",
        "## Verdict",
        "",
        "PASS_WITH_CONDITIONS",
        "",
        "## What passed",
        "",
        f"- Summarized candidate expression for {sample_cluster['gene_id'].nunique()} matched candidate genes across sample, encoded group and author broad cell label.",
        "- Used raw count layer when available and computed prevalence plus counts per cell.",
        "- Produced an adenomyosis localization matrix for all 102 GWAS candidate records, including unmatched genes with zero support.",
        "- Explicitly marks this dataset as scRNA/cell-label evidence, not spatial evidence.",
        "",
        "## Quantitative checks",
        "",
        f"- Sample-level comparison rows: {len(tests)}.",
        f"- Tests with q<=0.10 within comparison/metric families: {n_fdr}.",
        "",
        "## Limitations",
        "",
        "- Encoded group names (`C_EM`, `AM_EuM`, `AM_EcM`) are used conservatively and require source-paper verification before final biological wording.",
        "- There are only three samples per encoded group, so p-values are low-powered and should not drive target ranking alone.",
        "- The h5ad lacks detected spatial coordinates; spatial niche conclusions still require another source or literature-level support.",
        "- Cell labels are author-provided broad labels; subcluster/state-level claims require additional marker validation.",
        "",
        "## Decision",
        "",
        "This adenomyosis layer is suitable for cross-disease candidate localization and scoring as scRNA evidence. It is not sufficient for spatial claims or definitive adenomyosis differential expression.",
    ]
    OUT_REVIEW.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    candidates = load_candidate_metadata()
    adata = ad.read_h5ad(H5AD, backed="r")
    try:
        sample_cluster = summarize_by_sample_cluster(adata, candidates)
    finally:
        adata.file.close()
    group_summary = summarize_group_cluster(sample_cluster)
    tests = run_group_tests(sample_cluster)
    matrix = build_localization_matrix(candidates, group_summary, tests)

    sample_cluster.to_csv(OUT_BY_SAMPLE_CLUSTER, sep="\t", index=False)
    group_summary.to_csv(OUT_SUMMARY_GROUP_CLUSTER, sep="\t", index=False)
    tests.to_csv(OUT_TESTS, sep="\t", index=False)
    matrix.to_csv(OUT_MATRIX, sep="\t", index=False)
    write_summary(sample_cluster, group_summary, tests, matrix)
    write_self_review(sample_cluster, tests, matrix)
    print(OUT_SUMMARY)
    print(OUT_MATRIX)
    print(OUT_REVIEW)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
