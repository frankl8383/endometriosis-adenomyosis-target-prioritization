#!/usr/bin/env python3
"""Subject-level GSE203191 cell-state proportion validation."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy.stats import mannwhitneyu


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW = PROJECT_ROOT / "data" / "raw_downloads"
OUT = PROJECT_ROOT / "results" / "singlecell"


CELLSTATE_RULES = {
    "uNK_total": lambda s: s.startswith("uNK"),
    "uNK1": lambda s: s == "uNK1",
    "uNK2": lambda s: s == "uNK2",
    "stromal_total": lambda s: s.startswith("Stromal"),
    "myeloid_total": lambda s: s.startswith("Myeloid"),
    "b_cell": lambda s: s == "B",
    "epithelial_total": lambda s: s.startswith("Epithelial"),
    "cd4_t": lambda s: s == "CD4T",
    "cd8_t_total": lambda s: s.startswith("CD8T"),
}


def bh_adjust(pvalues: list[float]) -> list[float]:
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i])
    adjusted = [1.0] * n
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        original_rank = n - rank + 1
        value = min(prev, pvalues[idx] * n / original_rank)
        adjusted[idx] = value
        prev = value
    return adjusted


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(RAW / "GSE203191__GSE203191_Shih_endo_meta_frame.tsv.gz", sep="\t")
    total = df.groupby(["subjectID", "pheno"]).size().rename("total_cells").reset_index()
    rows: list[dict[str, object]] = []
    for _, subject in total.iterrows():
        sub = df[df["subjectID"] == subject["subjectID"]]
        record: dict[str, object] = {
            "subjectID": subject["subjectID"],
            "pheno": subject["pheno"],
            "total_cells": subject["total_cells"],
        }
        for state, rule in CELLSTATE_RULES.items():
            count = int(sub["clusterID"].map(rule).sum())
            record[f"{state}_cells"] = count
            record[f"{state}_fraction"] = count / subject["total_cells"]
        rows.append(record)
    prop = pd.DataFrame(rows)
    prop.to_csv(OUT / "GSE203191_cellstate_proportions_by_subject.tsv", sep="\t", index=False)

    tests: list[dict[str, object]] = []
    for state in CELLSTATE_RULES:
        col = f"{state}_fraction"
        control = prop.loc[prop["pheno"] == "Control", col].dropna()
        diagnosed = prop.loc[prop["pheno"] == "Diagnosed", col].dropna()
        symptomatic = prop.loc[prop["pheno"] == "Symptomatic", col].dropna()
        for comparison, a, b in [
            ("Diagnosed_vs_Control", diagnosed, control),
            ("Symptomatic_vs_Control", symptomatic, control),
            ("Diagnosed_vs_Symptomatic", diagnosed, symptomatic),
        ]:
            if len(a) < 2 or len(b) < 2:
                p = 1.0
                stat = float("nan")
            else:
                res = mannwhitneyu(a, b, alternative="two-sided")
                p = float(res.pvalue)
                stat = float(res.statistic)
            tests.append(
                {
                    "cell_state": state,
                    "comparison": comparison,
                    "n_group_a": len(a),
                    "n_group_b": len(b),
                    "median_group_a_fraction": float(a.median()) if len(a) else float("nan"),
                    "median_group_b_fraction": float(b.median()) if len(b) else float("nan"),
                    "delta_median_group_a_minus_b": float(a.median() - b.median()) if len(a) and len(b) else float("nan"),
                    "mannwhitney_u": stat,
                    "p_value": p,
                }
            )
    test_df = pd.DataFrame(tests)
    test_df["adj_p_value"] = bh_adjust(test_df["p_value"].tolist())
    test_df.to_csv(OUT / "GSE203191_cellstate_proportion_tests.tsv", sep="\t", index=False)

    focused = test_df[test_df["comparison"] == "Diagnosed_vs_Control"].sort_values("p_value")
    lines = [
        "# GSE203191 cell-state proportion validation",
        "",
        f"Subjects: {prop['subjectID'].nunique()}",
        f"Cells: {int(prop['total_cells'].sum())}",
        "",
        "## Diagnosed versus control",
        "",
        dataframe_to_markdown(
            focused[
                [
                    "cell_state",
                    "n_group_a",
                    "n_group_b",
                    "median_group_a_fraction",
                    "median_group_b_fraction",
                    "delta_median_group_a_minus_b",
                    "p_value",
                    "adj_p_value",
                ]
            ]
        ),
        "",
        "Interpretation: this validates cell-state composition signals only; expression-level candidate localization still requires matrix data.",
        "",
    ]
    (OUT / "GSE203191_cellstate_proportion_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print(OUT / "GSE203191_cellstate_proportion_summary.md")
    return 0


def dataframe_to_markdown(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for _, row in df.iterrows():
        vals = []
        for col in cols:
            value = row[col]
            if isinstance(value, float):
                value = f"{value:.4g}"
            vals.append(str(value).replace("|", "\\|"))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


if __name__ == "__main__":
    raise SystemExit(main())
