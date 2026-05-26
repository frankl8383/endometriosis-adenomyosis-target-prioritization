#!/usr/bin/env python3
"""Create shortlist tables and Figure 6."""

from __future__ import annotations

import os
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT_ROOT / "results"
INTEGRATION = RESULTS / "integration"
TABLES = RESULTS / "tables"
FIGURES = RESULTS / "figures"
GWAS = RESULTS / "gwas"
PY_PKG_DIR = PROJECT_ROOT / "tools" / "py_packages"

if sys.version_info[:2] != (3, 12):
    raise SystemExit(
        "Use the bundled Python 3.12 runtime: "
        "/Users/doctorliu/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 "
        "scripts/make_phase14_tables_and_figure6.py"
    )

sys.path.insert(0, str(PY_PKG_DIR))
os.environ.setdefault("MPLCONFIGDIR", str(FIGURES / ".mplconfig"))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402


REVIEW = INTEGRATION / "phase14_manual_target_directionality_review.tsv"
PRELIM = INTEGRATION / "preliminary_target_priority_matrix.tsv"
OUT_SHORTLIST_TSV = TABLES / "table_final_target_shortlist.tsv"
OUT_SHORTLIST_MD = TABLES / "table_final_target_shortlist.md"
OUT_CONTEXT_TSV = TABLES / "table_context_and_exploratory_candidates.tsv"
OUT_FIG_DATA = FIGURES / "figure6_experimental_hypothesis_triage_data.tsv"
OUT_FIG_SVG = FIGURES / "figure6_experimental_hypothesis_triage.svg"
OUT_FIG_PNG = FIGURES / "figure6_experimental_hypothesis_triage.png"
OUT_REVIEW = INTEGRATION / "phase15_figure6_table_self_review.md"

PALETTE = {
    "genetics": "#5B7FCA",
    "bulk": "#77B7A5",
    "singlecell": "#D69C4E",
    "druggability": "#B75F5F",
    "neutral_light": "#E7E7E7",
    "neutral_mid": "#8A8A8A",
    "neutral_dark": "#3F3F3F",
    "primary": "#315A9E",
    "safety": "#B75F5F",
    "secondary": "#6C8F5C",
    "context": "#8A8A8A",
}


def apply_publication_style() -> None:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["font.size"] = 8
    plt.rcParams["axes.spines.right"] = False
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.linewidth"] = 0.8
    plt.rcParams["legend.frameon"] = False


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.08, 1.04, label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def manuscript_use(row: pd.Series) -> str:
    tier = str(row["manual_evidence_tier"])
    category = str(row.get("manual_category", ""))
    gene = str(row.get("gene_symbol", ""))
    if tier == "A_primary" and gene == "LY96":
        return "pathway-level"
    if tier == "A_primary":
        return "primary"
    if tier == "B_primary_safety_limited":
        return "safety-limited"
    if category == "safety_limited_secondary_candidate":
        return "safety-limited"
    if tier == "C_secondary":
        return "secondary"
    if tier == "D_context_or_anchor":
        return "context"
    if tier == "E_marker_not_target":
        return "marker"
    return "exploratory"


def compact_citations(text: str, limit: int = 2) -> str:
    pmids = [item for item in str(text).split("|") if item]
    return ", ".join(pmids[:limit]) + ("..." if len(pmids) > limit else "")


def category_label(row: pd.Series) -> str:
    category = str(row.get("manual_category", ""))
    tier = str(row.get("manual_evidence_tier", ""))
    if category == "safety_limited_secondary_candidate":
        return "Safety-limited secondary"
    if tier == "A_primary" and row.get("gene_symbol") == "LY96":
        return "Primary pathway-level"
    if tier == "A_primary":
        return "Primary"
    if tier == "B_primary_safety_limited":
        return "Safety-limited fibrovascular signal"
    if tier == "C_secondary":
        return "Secondary"
    return tier or category


def clean_phrase(value: object) -> str:
    text = str(value)
    text = text.replace("_", " ")
    text = text.replace("scRNA", "scRNA-seq")
    text = text.replace("MD2", "MD-2")
    text = text.replace("VEGFR2", "VEGFR2")
    return text


def stability_label(value: object) -> str:
    mapping = {
        "rank_stable_target_hypothesis": "rank-stable prioritisation candidate",
        "moderately_stable_target_hypothesis": "moderately stable prioritisation candidate",
        "rank_sensitive_secondary_candidate": "rank-sensitive secondary candidate",
        "stable_marker_context_not_direct_target": "stable context marker, not direct target",
    }
    text = str(value)
    return mapping.get(text, clean_phrase(text).replace("target hypothesis", "prioritisation candidate"))


def soften_table_text(df: pd.DataFrame) -> pd.DataFrame:
    """Return a reviewer-facing copy with over-strong internal labels softened."""
    out = df.copy()
    replacements = {
        "rank_stable_target_hypothesis": "rank_stable_prioritisation_candidate",
        "moderately_stable_target_hypothesis": "moderately_stable_prioritisation_candidate",
        "therapeutic direction": "perturbation direction",
        "fibrovascular localization": "fibrovascular cell-context",
        "immune localization": "immune cell-context",
        "cell-state localization": "cell-context evidence",
        "target hypothesis": "model-system hypothesis",
    }
    for col in out.select_dtypes(include="object").columns:
        series = out[col].astype(str)
        for old, new in replacements.items():
            series = series.str.replace(old, new, regex=False)
        out[col] = series
    return out


def mapping_note(distance: object) -> str:
    value = pd.to_numeric(pd.Series([distance]), errors="coerce").iloc[0]
    if pd.isna(value):
        return "distance unavailable"
    dist = int(round(float(value)))
    kb = dist / 1000.0
    if dist <= 50_000:
        return f"lead-proximal ({kb:.1f} kb)"
    if dist <= 250_000:
        return f"<=250 kb ({kb:.0f} kb)"
    if dist <= 500_000:
        return f"<=500 kb ({kb:.0f} kb)"
    return f"distal; ±1 Mb only ({kb:.0f} kb)"


def add_mapping_context(df: pd.DataFrame) -> pd.DataFrame:
    universe_path = GWAS / "gwas_candidate_gene_universe.tsv"
    out = df.copy()
    if not universe_path.exists():
        out["distance_to_best_lead_bp"] = ""
        out["retained_at_250kb"] = ""
        out["retained_at_500kb"] = ""
        out["retained_at_1mb"] = ""
        out["gene_mapping_note"] = "mapping table unavailable"
        return out

    universe = pd.read_csv(universe_path, sep="\t", keep_default_na=False)
    cols = [
        "gene_id",
        "best_lead_snp",
        "best_phenotype",
        "distance_to_best_lead_bp",
        "genetic_priority",
        "ld_neighborhood_class",
    ]
    cols = [col for col in cols if col in universe.columns]
    out = out.merge(universe[cols].drop_duplicates("gene_id"), on="gene_id", how="left")
    distances = pd.to_numeric(out["distance_to_best_lead_bp"], errors="coerce")
    out["retained_at_250kb"] = distances.le(250_000).fillna(False)
    out["retained_at_500kb"] = distances.le(500_000).fillna(False)
    out["retained_at_1mb"] = distances.le(1_000_000).fillna(False)
    out["gene_mapping_note"] = out["distance_to_best_lead_bp"].apply(mapping_note)
    return out


def layer_summary(row: pd.Series) -> str:
    pieces = [
        f"genetics {float(row.get('genetics_support_score_20', 0) or 0):.0f}/20",
        f"bulk {float(row.get('bulk_expression_support_score_20', 0) or 0):.0f}/20",
        f"scRNA {float(row.get('cross_disease_scRNA_score_40', 0) or 0):.1f}/40",
        f"targetability {float(row.get('druggability_score_20', 0) or 0):.0f}/20",
    ]
    return "; ".join(pieces)


def interpretation(row: pd.Series) -> str:
    gene = str(row["gene_symbol"])
    if gene == "KDM1A":
        return "Epithelial-epigenetic model-system hypothesis; test KDM1A/LSD1 perturbation, including inhibition."
    if gene == "LY96":
        return "MD-2/TLR4 pathway-level immune hypothesis; LY96-specific evidence remains limited."
    if gene == "PDGFRA":
        return "Stromal-fibrovascular remodelling hypothesis; test receptor modulation in lesion-relevant systems."
    if gene == "KDR":
        return "Fibrovascular angiogenic signal; not a recommendation for systemic VEGFR2 inhibition."
    if gene == "KIT":
        return "Safety-limited immune/mast-cell axis hypothesis; perturbation evidence is not KIT-specific."
    if gene == "ECE1":
        return "Secondary endothelin vascular-pain axis hypothesis with limited ECE1-specific evidence."
    if gene == "C1QA":
        return "Secondary complement-axis hypothesis; context may determine protective versus pathogenic effects."
    return clean_phrase(row.get("proposed_direction", ""))


def caution(row: pd.Series) -> str:
    gene = str(row["gene_symbol"])
    if gene == "KDM1A":
        return "Broad epigenetic regulator; local, time-limited and fertility-aware validation needed."
    if gene == "LY96":
        return "Innate-immunity modulation may affect infection response and reproductive immune balance."
    if gene == "PDGFRA":
        return "Systemic growth-factor receptor inhibition can affect stromal repair and vascular homeostasis."
    if gene == "KDR":
        return "High vascular, wound-healing and reproductive safety concern for systemic anti-angiogenesis."
    if gene == "KIT":
        return "Haematological and mast-cell liabilities; multi-target CSF1R/KIT evidence only."
    if gene == "ECE1":
        return "Endothelin biology is vascular and pain-linked; direction needs functional testing."
    if gene == "C1QA":
        return "Complement modulation is systemic and stage/context dependent."
    return clean_phrase(row.get("manual_safety", ""))


def selected_evidence(row: pd.Series) -> str:
    pmids = [item for item in str(row.get("selected_pmids", "")).split("|") if item]
    if not pmids:
        return "Targeted literature sparse"
    return f"Targeted literature support ({len(pmids)} records; full identifiers in supplementary table)"


def claim_label(row: pd.Series) -> str:
    gene = str(row["gene_symbol"])
    if gene == "KDM1A":
        return "Primary model-system hypothesis; not causal proof"
    if gene == "PDGFRA":
        return "Primary fibrovascular-stromal follow-up hypothesis"
    if gene == "LY96":
        return "Pathway-level hypothesis"
    if gene == "KDR":
        return "Safety-limited fibrovascular signal"
    if gene == "KIT":
        return "Safety-limited secondary hypothesis"
    if gene in {"ECE1", "C1QA"}:
        return "Secondary/exploratory axis"
    raw = str(row.get("claim_strength", "")).strip()
    if raw:
        return clean_phrase(raw).replace("target hypothesis", "experimental hypothesis")
    if gene in {"KDM1A", "PDGFRA"}:
        return "Primary experimental hypothesis"
    return "Hypothesis for functional testing"


def write_tables(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    TABLES.mkdir(parents=True, exist_ok=True)
    df = add_mapping_context(df)
    shortlist = df[df["manual_evidence_tier"].isin(["A_primary", "B_primary_safety_limited", "C_secondary"])].copy()
    shortlist = shortlist.sort_values(["manual_evidence_tier", "final_target_priority_score_100_pre_rank_stability"], ascending=[True, False])
    shortlist_cols = [
        "manual_evidence_tier",
        "manual_category",
        "gene_symbol",
        "best_lead_snp",
        "best_phenotype",
        "distance_to_best_lead_bp",
        "retained_at_250kb",
        "retained_at_500kb",
        "retained_at_1mb",
        "gene_mapping_note",
        "proposed_direction",
        "dominant_cross_disease_axis",
        "rank_stability_class",
        "final_target_priority_score_100_pre_rank_stability",
        "genetics_support_score_20",
        "bulk_expression_support_score_20",
        "cross_disease_scRNA_score_40",
        "druggability_score_20",
        "bootstrap_top10_frequency",
        "leave_one_layer_top20_count",
        "gene_label_permutation_p_ge_observed",
        "manual_safety",
        "literature_level",
        "selected_pmids",
        "manual_rationale",
        "claim_strength",
    ]
    soften_table_text(shortlist[shortlist_cols]).to_csv(OUT_SHORTLIST_TSV, sep="\t", index=False)

    md_lines = [
        "| category | gene | mapping note | evidence/stability | permitted claim | main caution |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for _, row in shortlist.iterrows():
        score_text = (
            f"{float(row['final_target_priority_score_100_pre_rank_stability']):.1f}; "
            f"{stability_label(row['rank_stability_class'])}"
        )
        permitted_claim = f"{claim_label(row)}; {interpretation(row)}"
        md_lines.append(
            "| "
            + " | ".join(
                [
                    category_label(row),
                    str(row["gene_symbol"]),
                    str(row["gene_mapping_note"]),
                    score_text,
                    permitted_claim,
                    caution(row),
                ]
            )
            + " |"
        )
    OUT_SHORTLIST_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    context = df[~df["manual_evidence_tier"].isin(["A_primary", "B_primary_safety_limited", "C_secondary"])].copy()
    context_cols = [
        "manual_evidence_tier",
        "manual_category",
        "gene_symbol",
        "best_lead_snp",
        "best_phenotype",
        "distance_to_best_lead_bp",
        "retained_at_250kb",
        "retained_at_500kb",
        "retained_at_1mb",
        "gene_mapping_note",
        "dominant_cross_disease_axis",
        "literature_level",
        "selected_pmids",
        "manual_rationale",
        "claim_strength",
    ]
    soften_table_text(context[context_cols]).to_csv(OUT_CONTEXT_TSV, sep="\t", index=False)
    return shortlist, context


def plot_figure(df: pd.DataFrame) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    apply_publication_style()
    df = add_mapping_context(df)
    plot_df = df[df["manual_evidence_tier"].isin(["A_primary", "B_primary_safety_limited", "C_secondary", "D_context_or_anchor"])].copy()
    plot_df = plot_df.sort_values(["manual_evidence_tier", "final_target_priority_score_100_pre_rank_stability"], ascending=[False, True])
    soften_table_text(plot_df).to_csv(OUT_FIG_DATA, sep="\t", index=False)

    genes = plot_df["gene_symbol"].tolist()
    y = np.arange(len(plot_df))

    fig = plt.figure(figsize=(7.2, 5.8), constrained_layout=False)
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.35, 1.0],
        height_ratios=[1.2, 0.95],
        left=0.10,
        right=0.94,
        top=0.93,
        bottom=0.12,
        hspace=0.42,
        wspace=0.34,
    )
    ax_bar = fig.add_subplot(gs[:, 0])
    ax_heat = fig.add_subplot(gs[0, 1])
    ax_cat = fig.add_subplot(gs[1, 1])

    layers = [
        ("genetics_support_score_20", "Genetics", PALETTE["genetics"]),
        ("bulk_expression_support_score_20", "Bulk", PALETTE["bulk"]),
        ("cross_disease_scRNA_score_40", "scRNA", PALETTE["singlecell"]),
        ("druggability_score_20", "Targetability", PALETTE["druggability"]),
    ]
    left = np.zeros(len(plot_df))
    for col, label, color in layers:
        values = pd.to_numeric(plot_df[col], errors="coerce").fillna(0).to_numpy()
        ax_bar.barh(y, values, left=left, height=0.64, color=color, edgecolor="white", linewidth=0.5, label=label)
        left += values
    ax_bar.set_yticks(y)
    ax_bar.set_yticklabels(genes)
    ax_bar.set_xlabel("Evidence score")
    ax_bar.set_xlim(0, 100)
    ax_bar.grid(axis="x", color="#D8D8D8", linewidth=0.5)
    ax_bar.legend(loc="lower right", fontsize=7, ncols=1)
    ax_bar.set_title("Evidence layers")
    add_panel_label(ax_bar, "a")

    heat_values = np.column_stack(
        [
            pd.to_numeric(plot_df["bootstrap_top10_frequency"], errors="coerce").fillna(0).to_numpy(),
            pd.to_numeric(plot_df["leave_one_layer_top20_count"], errors="coerce").fillna(0).to_numpy() / 4.0,
            np.clip(
                -np.log10(pd.to_numeric(plot_df["gene_label_permutation_p_ge_observed"], errors="coerce").fillna(1).to_numpy()),
                0,
                2,
            )
            / 2.0,
        ]
    )
    im = ax_heat.imshow(heat_values, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax_heat.set_yticks(y)
    ax_heat.set_yticklabels(genes, fontsize=7)
    ax_heat.set_xticks([0, 1, 2])
    ax_heat.set_xticklabels(["Bootstrap", "Layer\nrobust", "Null\nscore"], fontsize=7)
    ax_heat.tick_params(length=0)
    ax_heat.set_title("Rank stability")
    add_panel_label(ax_heat, "b")
    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.02)
    cbar.ax.tick_params(labelsize=7, length=2)
    cbar.set_label("Scaled", fontsize=7)

    use_colors = {
        "primary": PALETTE["primary"],
        "pathway-level": "#4D76B8",
        "safety-limited": PALETTE["safety"],
        "secondary": PALETTE["secondary"],
        "context": PALETTE["neutral_mid"],
        "marker": PALETTE["neutral_light"],
        "exploratory": PALETTE["neutral_light"],
    }
    uses = plot_df.apply(manuscript_use, axis=1).tolist()
    cat_colors = [use_colors[item] for item in uses]
    ax_cat.barh(y, [1] * len(y), height=0.58, color=cat_colors, edgecolor="white", linewidth=0.5)
    ax_cat.set_yticks(y)
    ax_cat.set_yticklabels(genes, fontsize=7)
    ax_cat.set_xticks([])
    ax_cat.set_xlim(0, 1)
    ax_cat.set_title("Follow-up class")
    add_panel_label(ax_cat, "c")
    for yi, use, color in zip(y, uses, cat_colors):
        text_color = "white" if use in {"primary", "pathway-level", "safety-limited", "secondary", "context"} else PALETTE["neutral_dark"]
        ax_cat.text(0.97, yi, use, ha="right", va="center", fontsize=7, color=text_color)

    for ax in [ax_heat, ax_cat]:
        for spine in ax.spines.values():
            spine.set_visible(False)

    fig.savefig(OUT_FIG_SVG)
    fig.savefig(OUT_FIG_PNG, dpi=1200)
    plt.close(fig)


def write_self_review(shortlist: pd.DataFrame, context: pd.DataFrame) -> None:
    failures: list[str] = []
    for path in [OUT_SHORTLIST_TSV, OUT_SHORTLIST_MD, OUT_CONTEXT_TSV, OUT_FIG_DATA, OUT_FIG_SVG, OUT_FIG_PNG]:
        if not path.exists() or path.stat().st_size == 0:
            failures.append(f"Missing or empty output: {path}")
    if len(shortlist) < 5:
        failures.append("Shortlist table has fewer than five candidates.")
    if "HSPG2" in set(shortlist["gene_symbol"]):
        failures.append("HSPG2 incorrectly appears in the direct shortlist.")
    if "KDR" in set(shortlist["gene_symbol"]):
        safety = shortlist.loc[shortlist["gene_symbol"] == "KDR", "manual_safety"].iloc[0]
        if "caution" not in str(safety).lower():
            failures.append("KDR appears without explicit safety caution.")
    if "KIT" in set(shortlist["gene_symbol"]):
        tier = shortlist.loc[shortlist["gene_symbol"] == "KIT", "manual_evidence_tier"].iloc[0]
        category = shortlist.loc[shortlist["gene_symbol"] == "KIT", "manual_category"].iloc[0]
        if tier != "C_secondary" or "safety_limited_secondary" not in str(category):
            failures.append("KIT should appear as a safety-limited secondary candidate.")
    if not OUT_FIG_SVG.read_text(encoding="utf-8", errors="ignore").strip().startswith("<?xml"):
        failures.append("Figure 6 SVG does not look like a valid SVG/XML file.")

    status = "PASS" if not failures else "FAIL"
    lines = [
        "# Phase 15 self-review: shortlist tables and Figure 6 draft",
        "",
        f"Status: {status}",
        "",
        "## Checks",
        "",
        f"- Direct/secondary shortlist rows: {len(shortlist)}",
        f"- Context/exploratory rows: {len(context)}",
        f"- Figure SVG size: {OUT_FIG_SVG.stat().st_size if OUT_FIG_SVG.exists() else 0} bytes",
        f"- Figure PNG size: {OUT_FIG_PNG.stat().st_size if OUT_FIG_PNG.exists() else 0} bytes",
        "",
        "## Guardrails",
        "",
        "- The figure separates manual use class from numerical evidence score.",
        "- Safety-limited genes are not visually indistinguishable from unrestricted primary candidates.",
        "- Context genes are preserved in tables but excluded from the direct target shortlist.",
        "",
    ]
    if failures:
        lines.extend(["## Failures", ""])
        lines.extend(f"- {failure}" for failure in failures)
    else:
        lines.extend(["## Decision", "", "The shortlist tables and Figure 6 draft pass automated checks and are ready for visual/manual refinement."])
    OUT_REVIEW.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    FIGURES.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(REVIEW, sep="\t", keep_default_na=False)
    prelim = pd.read_csv(PRELIM, sep="\t", keep_default_na=False)
    needed = [
        "gene_id",
        "genetics_support_score_20",
        "bulk_expression_support_score_20",
        "cross_disease_scRNA_score_40",
        "druggability_score_20",
    ]
    available = [col for col in needed if col in prelim.columns]
    df = df.merge(prelim[available], on="gene_id", how="left", suffixes=("", "_from_prelim"))
    for col in needed[1:]:
        alt = f"{col}_from_prelim"
        if col not in df.columns and alt in df.columns:
            df[col] = df[alt]
        elif col in df.columns and alt in df.columns:
            df[col] = df[col].where(df[col].astype(str).str.len() > 0, df[alt])
    shortlist, context = write_tables(df)
    plot_figure(df)
    write_self_review(shortlist, context)
    print(f"Wrote {OUT_SHORTLIST_TSV}")
    print(f"Wrote {OUT_SHORTLIST_MD}")
    print(f"Wrote {OUT_CONTEXT_TSV}")
    print(f"Wrote {OUT_FIG_SVG}")
    print(f"Wrote {OUT_FIG_PNG}")
    print(f"Wrote {OUT_REVIEW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
