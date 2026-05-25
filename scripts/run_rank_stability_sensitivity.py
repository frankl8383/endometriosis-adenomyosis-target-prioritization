#!/usr/bin/env python3
"""Run rank-stability sensitivity analysis for target prioritisation."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT_ROOT / "results"
INTEGRATION = RESULTS / "integration"
PY_PKG_DIR = PROJECT_ROOT / "tools" / "py_packages"

if sys.version_info[:2] != (3, 12):
    raise SystemExit(
        "Use the bundled Python 3.12 runtime: "
        "python3 "
        "scripts/run_rank_stability_sensitivity.py"
    )

sys.path.insert(0, str(PY_PKG_DIR))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


INPUT = INTEGRATION / "preliminary_target_priority_matrix.tsv"
OUT_LAYER = INTEGRATION / "rank_stability_leave_one_layer_out.tsv"
OUT_BOOT = INTEGRATION / "rank_stability_weight_bootstrap.tsv"
OUT_PERM = INTEGRATION / "rank_stability_gene_label_permutation.tsv"
OUT_COMBINED = INTEGRATION / "rank_stability_combined_matrix.tsv"
OUT_SUMMARY = INTEGRATION / "rank_stability_summary.md"
OUT_REVIEW = INTEGRATION / "phase13_rank_stability_self_review.md"

LAYERS = {
    "genetics": ("genetics_support_score_20", 20.0),
    "bulk": ("bulk_expression_support_score_20", 20.0),
    "singlecell": ("cross_disease_scRNA_score_40", 40.0),
    "druggability": ("druggability_score_20", 20.0),
}


def to_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").fillna(0)


def add_layer_fractions(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for layer, (col, weight) in LAYERS.items():
        out[col] = to_num(out[col])
        out[f"{layer}_fraction"] = (out[col] / weight).clip(0, 1)
    return out


def rank_from_score(df: pd.DataFrame, score_col: str, rank_col: str) -> pd.DataFrame:
    out = df.sort_values([score_col, "druggability_score_20"], ascending=False).copy()
    out[rank_col] = np.arange(1, len(out) + 1)
    return out[["gene_id", "gene_symbol", score_col, rank_col]]


def leave_one_layer_out(df: pd.DataFrame) -> pd.DataFrame:
    base_cols = ["gene_id", "gene_symbol", "pre_rank_stability_rank"]
    rows = df[base_cols].copy()
    total_weight = sum(weight for _, weight in LAYERS.values())
    for held_out, (_, held_weight) in LAYERS.items():
        included = [layer for layer in LAYERS if layer != held_out]
        score = np.zeros(len(df))
        included_weight = total_weight - held_weight
        for layer in included:
            col, weight = LAYERS[layer]
            score += to_num(df[col])
        score = score / included_weight * 100.0
        tmp = df[["gene_id", "gene_symbol", "druggability_score_20"]].copy()
        tmp[f"score_without_{held_out}_rescaled_100"] = score.round(6)
        ranked = rank_from_score(tmp, f"score_without_{held_out}_rescaled_100", f"rank_without_{held_out}")
        rows = rows.merge(ranked.drop(columns=["gene_symbol"]), on="gene_id", how="left")
    rank_cols = [col for col in rows.columns if col.startswith("rank_without_")]
    rows["leave_one_layer_top10_count"] = rows[rank_cols].le(10).sum(axis=1)
    rows["leave_one_layer_top20_count"] = rows[rank_cols].le(20).sum(axis=1)
    rows["leave_one_layer_best_rank"] = rows[rank_cols].min(axis=1)
    rows["leave_one_layer_worst_rank"] = rows[rank_cols].max(axis=1)
    rows["leave_one_layer_rank_range"] = rows["leave_one_layer_worst_rank"] - rows["leave_one_layer_best_rank"]
    return rows.sort_values(["leave_one_layer_top10_count", "leave_one_layer_top20_count", "leave_one_layer_best_rank"], ascending=[False, False, True])


def bootstrap_weights(df: pd.DataFrame, n_boot: int, top_n: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    layer_names = list(LAYERS)
    base_weights = np.array([LAYERS[layer][1] for layer in layer_names], dtype=float)
    fractions = np.column_stack([df[f"{layer}_fraction"].to_numpy(dtype=float) for layer in layer_names])
    top_counts = np.zeros(len(df), dtype=int)
    top20_counts = np.zeros(len(df), dtype=int)
    ranks = np.empty((n_boot, len(df)), dtype=np.int16)

    for i in range(n_boot):
        jitter = rng.uniform(0.8, 1.2, size=len(layer_names))
        weights = base_weights * jitter
        weights = weights / weights.sum() * 100.0
        scores = fractions @ weights
        order = np.argsort(-scores, kind="mergesort")
        rank = np.empty(len(df), dtype=np.int16)
        rank[order] = np.arange(1, len(df) + 1)
        ranks[i] = rank
        top_counts[rank <= top_n] += 1
        top20_counts[rank <= 20] += 1

    out = df[["gene_id", "gene_symbol"]].copy()
    out["bootstrap_top10_frequency"] = top_counts / n_boot
    out["bootstrap_top20_frequency"] = top20_counts / n_boot
    out["bootstrap_median_rank"] = np.median(ranks, axis=0)
    out["bootstrap_q10_rank"] = np.quantile(ranks, 0.10, axis=0)
    out["bootstrap_q90_rank"] = np.quantile(ranks, 0.90, axis=0)
    out["bootstrap_rank_iqr_like_width"] = out["bootstrap_q90_rank"] - out["bootstrap_q10_rank"]
    return out.sort_values(["bootstrap_top10_frequency", "bootstrap_median_rank"], ascending=[False, True])


def gene_label_permutation(df: pd.DataFrame, n_perm: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed + 100_000)
    layer_cols = [LAYERS[layer][0] for layer in LAYERS]
    observed = to_num(df["final_target_priority_score_100_pre_rank_stability"]).to_numpy(dtype=float)
    exceed = np.zeros(len(df), dtype=int)
    max_scores = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        perm_score = np.zeros(len(df), dtype=float)
        for col in layer_cols:
            values = to_num(df[col]).to_numpy(dtype=float).copy()
            rng.shuffle(values)
            perm_score += values
        exceed += perm_score >= observed
        max_scores[i] = perm_score.max()
    out = df[["gene_id", "gene_symbol"]].copy()
    out["gene_label_permutation_p_ge_observed"] = (exceed + 1) / (n_perm + 1)
    out["observed_exceeds_95pct_null_max"] = observed > np.quantile(max_scores, 0.95)
    out["observed_total_score"] = observed
    return out.sort_values(["gene_label_permutation_p_ge_observed", "observed_total_score"], ascending=[True, False])


def stability_class(row: pd.Series) -> str:
    tier = str(row.get("priority_tier", ""))
    if "not_direct_target" in tier:
        return "stable_marker_context_not_direct_target" if row["bootstrap_top20_frequency"] >= 0.5 else "context_unstable"
    if (
        row["bootstrap_top10_frequency"] >= 0.70
        and row["leave_one_layer_top20_count"] >= 3
        and row["gene_label_permutation_p_ge_observed"] <= 0.10
    ):
        return "rank_stable_target_hypothesis"
    if row["bootstrap_top10_frequency"] >= 0.40 and row["leave_one_layer_top20_count"] >= 2:
        return "moderately_stable_target_hypothesis"
    if row["bootstrap_top20_frequency"] >= 0.40:
        return "rank_sensitive_secondary_candidate"
    return "rank_unstable_or_context_only"


def combine(df: pd.DataFrame, layer: pd.DataFrame, boot: pd.DataFrame, perm: pd.DataFrame) -> pd.DataFrame:
    keep_cols = [
        "pre_rank_stability_rank",
        "gene_id",
        "gene_symbol",
        "gene_description",
        "final_target_priority_score_100_pre_rank_stability",
        "pre_druggability_biologic_evidence_score_80",
        "druggability_score_20",
        "priority_tier",
        "manuscript_role",
        "suggested_shortlist_label",
        "review_flags",
        "genetic_priority",
        "bulk_support_class",
        "cross_disease_scRNA_class",
        "dominant_cross_disease_axis",
        "actionability_class",
        "drug_direction_evidence",
        "druggability_penalty_reasons",
    ]
    out = df[keep_cols].merge(layer.drop(columns=["gene_symbol", "pre_rank_stability_rank"]), on="gene_id", how="left")
    out = out.merge(boot.drop(columns=["gene_symbol"]), on="gene_id", how="left")
    out = out.merge(perm.drop(columns=["gene_symbol"]), on="gene_id", how="left")
    out["rank_stability_class"] = out.apply(stability_class, axis=1)
    class_order = {
        "rank_stable_target_hypothesis": 1,
        "moderately_stable_target_hypothesis": 2,
        "rank_sensitive_secondary_candidate": 3,
        "stable_marker_context_not_direct_target": 4,
        "rank_unstable_or_context_only": 5,
        "context_unstable": 6,
    }
    out["rank_stability_order"] = out["rank_stability_class"].map(class_order).fillna(99).astype(int)
    out = out.sort_values(
        [
            "rank_stability_order",
            "bootstrap_top10_frequency",
            "final_target_priority_score_100_pre_rank_stability",
        ],
        ascending=[True, False, False],
    )
    return out


def write_summary(combined: pd.DataFrame, args: argparse.Namespace) -> None:
    class_counts = combined["rank_stability_class"].value_counts().to_dict()
    stable = combined[
        combined["rank_stability_class"].isin(["rank_stable_target_hypothesis", "moderately_stable_target_hypothesis"])
    ].sort_values(
        ["rank_stability_order", "bootstrap_top10_frequency", "final_target_priority_score_100_pre_rank_stability"],
        ascending=[True, False, False],
    )
    lines = [
        "# Rank-stability sensitivity summary",
        "",
        f"- Candidate rows: {len(combined)}",
        f"- Bootstrap iterations: {args.n_boot}",
        f"- Gene-label permutation iterations: {args.n_perm}",
        f"- Stability classes: {class_counts}",
        "",
        "## Stable or moderately stable target hypotheses",
        "",
        "| pre-rank | gene | score | boot top10 | leave-one-layer top20 | permutation p | stability | role |",
        "| ---: | --- | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    for _, row in stable.head(25).iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    str(int(row["pre_rank_stability_rank"])),
                    str(row["gene_symbol"] or row["gene_id"]),
                    f"{float(row['final_target_priority_score_100_pre_rank_stability']):.3f}",
                    f"{float(row['bootstrap_top10_frequency']):.3f}",
                    str(int(row["leave_one_layer_top20_count"])),
                    f"{float(row['gene_label_permutation_p_ge_observed']):.4f}",
                    str(row["rank_stability_class"]),
                    str(row["manuscript_role"]),
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "## Method notes",
            "",
            "- Leave-one-layer-out ranks remove genetics, bulk, single-cell or druggability one at a time and rescale remaining layers to 100.",
            "- Bootstrap ranks perturb layer weights by +/-20% and renormalize to 100.",
            "- Gene-label permutation is an internal concordance null; it is not a genome-wide matched-background null.",
            "- Stability results are for prioritisation robustness, not causal inference.",
        ]
    )
    OUT_SUMMARY.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_review(df: pd.DataFrame, layer: pd.DataFrame, boot: pd.DataFrame, perm: pd.DataFrame, combined: pd.DataFrame, args: argparse.Namespace) -> None:
    failures: list[str] = []
    if len(combined) != 102:
        failures.append(f"Expected 102 combined rows, observed {len(combined)}.")
    if len(layer) != 102 or len(boot) != 102 or len(perm) != 102:
        failures.append("One or more sensitivity tables do not contain 102 rows.")
    if combined["bootstrap_top10_frequency"].min() < 0 or combined["bootstrap_top10_frequency"].max() > 1:
        failures.append("Bootstrap top10 frequency is outside [0, 1].")
    if combined["gene_label_permutation_p_ge_observed"].min() <= 0 or combined["gene_label_permutation_p_ge_observed"].max() > 1:
        failures.append("Permutation p-values are outside (0, 1].")
    if int((combined["rank_stability_class"] == "rank_stable_target_hypothesis").sum()) == 0:
        failures.append("No rank-stable target hypotheses found; inspect whether thresholds are too strict or evidence is unstable.")
    hspg2 = combined[combined["gene_symbol"] == "HSPG2"]
    if not hspg2.empty and "direct_target" not in str(hspg2.iloc[0]["rank_stability_class"]):
        pass
    if not hspg2.empty and hspg2.iloc[0]["rank_stability_class"] == "rank_stable_target_hypothesis":
        failures.append("HSPG2 is incorrectly marked as a rank-stable direct target hypothesis.")

    status = "PASS" if not failures else "FAIL"
    lines = [
        "# Phase 13 self-review: rank-stability sensitivity",
        "",
        f"Status: {status}",
        "",
        "## Checks",
        "",
        f"- Candidate rows: {len(combined)}",
        f"- Bootstrap iterations: {args.n_boot}",
        f"- Gene-label permutation iterations: {args.n_perm}",
        f"- Rank-stable target hypotheses: {int((combined['rank_stability_class'] == 'rank_stable_target_hypothesis').sum())}",
        f"- Moderately stable target hypotheses: {int((combined['rank_stability_class'] == 'moderately_stable_target_hypothesis').sum())}",
        f"- Stable/context-not-direct genes: {int((combined['rank_stability_class'] == 'stable_marker_context_not_direct_target').sum())}",
        "",
        "## Guardrails",
        "",
        "- Permutation p-values test internal layer concordance only and should not be reported as genome-wide significance.",
        "- Stable rank does not override druggability penalties or systemic safety cautions.",
        "- Final shortlist still needs manual mechanism-direction review.",
        "",
    ]
    if failures:
        lines.extend(["## Failures", ""])
        lines.extend(f"- {failure}" for failure in failures)
    else:
        lines.extend(["## Decision", "", "Rank-stability analysis passes automated checks and can be used for manual shortlist review."])
    OUT_REVIEW.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-boot", type=int, default=2000)
    parser.add_argument("--n-perm", type=int, default=2000)
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260523)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    df = pd.read_csv(INPUT, sep="\t", keep_default_na=False)
    df = add_layer_fractions(df)

    layer = leave_one_layer_out(df)
    boot = bootstrap_weights(df, args.n_boot, args.top_n, args.seed)
    perm = gene_label_permutation(df, args.n_perm, args.seed)
    combined = combine(df, layer, boot, perm)

    layer.to_csv(OUT_LAYER, sep="\t", index=False)
    boot.to_csv(OUT_BOOT, sep="\t", index=False)
    perm.to_csv(OUT_PERM, sep="\t", index=False)
    combined.to_csv(OUT_COMBINED, sep="\t", index=False)
    write_summary(combined, args)
    write_review(df, layer, boot, perm, combined, args)

    print(f"Wrote {OUT_LAYER}")
    print(f"Wrote {OUT_BOOT}")
    print(f"Wrote {OUT_PERM}")
    print(f"Wrote {OUT_COMBINED}")
    print(f"Wrote {OUT_SUMMARY}")
    print(f"Wrote {OUT_REVIEW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
