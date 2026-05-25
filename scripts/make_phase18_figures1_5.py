#!/usr/bin/env python3
"""Create draft manuscript Figures 1-5 from audited project outputs."""

from __future__ import annotations

import ast
import os
import re
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT_ROOT / "results"
FIGURES = RESULTS / "figures"
GWAS = RESULTS / "gwas"
BULK = RESULTS / "bulk"
SC = RESULTS / "singlecell"
INTEGRATION = RESULTS / "integration"
PY_PKG_DIR = PROJECT_ROOT / "tools" / "py_packages"

if sys.version_info[:2] != (3, 12):
    raise SystemExit(
        "Use the bundled Python 3.12 runtime: "
        "python3 "
        "scripts/make_phase18_figures1_5.py"
    )

sys.path.insert(0, str(PY_PKG_DIR))
os.environ.setdefault("MPLCONFIGDIR", str(FIGURES / ".mplconfig"))

import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch  # noqa: E402


PALETTE = {
    "genetics": "#5B7FCA",
    "bulk": "#77B7A5",
    "singlecell": "#D69C4E",
    "action": "#B75F5F",
    "neutral_0": "#F2F2F2",
    "neutral_1": "#D8D8D8",
    "neutral_2": "#8A8A8A",
    "neutral_3": "#3F3F3F",
    "fibro": "#4C78A8",
    "immune": "#B279A2",
    "epithelial": "#72B7B2",
    "mixed": "#E3B35A",
    "limited": "#BEBEBE",
    "warning": "#C46A46",
}

FIGURE_PATHS = {
    "figure1": (FIGURES / "figure1_study_design_audit.svg", FIGURES / "figure1_study_design_audit.png"),
    "figure2": (FIGURES / "figure2_genetic_candidate_universe.svg", FIGURES / "figure2_genetic_candidate_universe.png"),
    "figure3": (FIGURES / "figure3_bulk_expression_support.svg", FIGURES / "figure3_bulk_expression_support.png"),
    "figure4": (FIGURES / "figure4_singlecell_celllabel_context.svg", FIGURES / "figure4_singlecell_celllabel_context.png"),
    "figure5": (FIGURES / "figure5_cross_disease_convergence.svg", FIGURES / "figure5_cross_disease_convergence.png"),
}
CAPTIONS = FIGURES / "figure1_5_captions.md"
SELF_REVIEW = INTEGRATION / "phase18_figures1_5_self_review.md"


def apply_publication_style() -> None:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["font.size"] = 7.5
    plt.rcParams["axes.spines.right"] = False
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.linewidth"] = 0.8
    plt.rcParams["legend.frameon"] = False
    plt.rcParams["figure.facecolor"] = "white"


def add_panel_label(ax: plt.Axes, label: str, x: float = -0.08, y: float = 1.03) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def save_figure(fig: plt.Figure, name: str) -> None:
    svg, png = FIGURE_PATHS[name]
    fig.savefig(svg)
    fig.savefig(png, dpi=300)
    plt.close(fig)


def clean_label(value: str, max_len: int = 28) -> str:
    value = str(value).replace("_", " ")
    if len(value) <= max_len:
        return value
    return value[: max_len - 1] + "…"


def parse_md_int(pattern: str, text: str) -> int:
    match = re.search(pattern, text)
    if not match:
        raise ValueError(f"Could not parse integer with pattern: {pattern}")
    return int(match.group(1).replace(",", ""))


def parse_literal_counts(label: str, text: str) -> dict[str, int]:
    match = re.search(re.escape(label) + r":\s*`({.*?})`", text)
    if not match:
        raise ValueError(f"Could not parse count dictionary for {label}")
    return {str(k): int(v) for k, v in ast.literal_eval(match.group(1)).items()}


def row_for(df: pd.DataFrame, column: str, value: str) -> pd.Series:
    hit = df[df[column].astype(str) == value]
    if hit.empty:
        raise ValueError(f"No row where {column} == {value}")
    return hit.iloc[0]


def horizontal_bar(ax: plt.Axes, labels: list[str], values: list[float], colors: list[str], xlabel: str) -> None:
    y = np.arange(len(labels))
    ax.barh(y, values, color=colors, edgecolor="white", linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.grid(axis="x", color="#E4E4E4", linewidth=0.5)


def draw_flow(ax: plt.Axes, nodes: list[tuple[str, str, str]]) -> None:
    ax.set_axis_off()
    xs = np.linspace(0.06, 0.94, len(nodes))
    y = 0.55
    w = 0.145
    h = 0.34
    for i, (title, subtitle, color) in enumerate(nodes):
        x = xs[i]
        box = FancyBboxPatch(
            (x - w / 2, y - h / 2),
            w,
            h,
            boxstyle="round,pad=0.012,rounding_size=0.018",
            facecolor=color,
            edgecolor="white",
            linewidth=1.0,
        )
        ax.add_patch(box)
        text_color = "white" if color in {PALETTE["genetics"], PALETTE["action"], PALETTE["fibro"]} else PALETTE["neutral_3"]
        ax.text(x, y + 0.045, title, ha="center", va="center", fontsize=8, fontweight="bold", color=text_color)
        ax.text(x, y - 0.075, subtitle, ha="center", va="center", fontsize=6.5, color=text_color)
        if i < len(nodes) - 1:
            arrow = FancyArrowPatch(
                (x + w / 2 + 0.006, y),
                (xs[i + 1] - w / 2 - 0.006, y),
                arrowstyle="-|>",
                mutation_scale=10,
                linewidth=0.9,
                color=PALETTE["neutral_2"],
            )
            ax.add_patch(arrow)


def load_sources() -> dict[str, object]:
    sources: dict[str, object] = {}
    sources["gwas_signal"] = pd.read_csv(GWAS / "gwas_signal_summary.tsv", sep="\t")
    sources["coord_md"] = (GWAS / "coordinate_mapping_audit.md").read_text(encoding="utf-8")
    sources["neighborhoods"] = pd.read_csv(GWAS / "coordinate_window_neighborhoods.tsv", sep="\t")
    sources["ld"] = pd.read_csv(GWAS / "ld_supported_shared_neighborhoods.tsv", sep="\t")
    sources["universe"] = pd.read_csv(GWAS / "gwas_candidate_gene_universe.tsv", sep="\t")
    sources["bulk_meta"] = pd.read_csv(BULK / "bulk_metadata_audit.tsv", sep="\t")
    sources["bulk_extract"] = pd.read_csv(BULK / "candidate_expression_extraction_audit.tsv", sep="\t")
    sources["bulk_models"] = pd.read_csv(BULK / "models" / "bulk_candidate_model_summary.tsv", sep="\t")
    sources["bulk_scores"] = pd.read_csv(BULK / "bulk_candidate_expression_support_scores.tsv", sep="\t")
    sources["gse179640_md"] = (SC / "GSE179640_broad_compartment_annotation_summary.md").read_text(encoding="utf-8")
    sources["gse179640_v2_md"] = (SC / "GSE179640_singlecell_candidate_evidence_matrix_v2_summary.md").read_text(encoding="utf-8")
    sources["gse203191"] = pd.read_csv(SC / "GSE203191_cellstate_proportion_tests.tsv", sep="\t")
    sources["zenodo_md"] = (SC / "Zenodo17078290_candidate_expression_summary.md").read_text(encoding="utf-8")
    sources["zenodo_audit_md"] = (SC / "Zenodo17078290_h5ad_audit.md").read_text(encoding="utf-8")
    sources["cross"] = pd.read_csv(INTEGRATION / "cross_disease_singlecell_localization_matrix.tsv", sep="\t")
    sources["cross_md"] = (INTEGRATION / "cross_disease_singlecell_localization_summary.md").read_text(encoding="utf-8")
    return sources


def make_figure1(s: dict[str, object]) -> None:
    gwas = s["gwas_signal"]
    bulk_meta = s["bulk_meta"]
    coord_md = s["coord_md"]
    gse179640_md = s["gse179640_md"]
    zenodo_audit_md = s["zenodo_audit_md"]
    gse203191 = s["gse203191"]

    gse179640_cells = parse_md_int(r"QC-passing cells annotated:\s*([0-9,]+)", gse179640_md)
    zenodo_shape = re.search(r"Shape:\s*([0-9,]+) observations x ([0-9,]+) variables", zenodo_audit_md)
    if not zenodo_shape:
        raise ValueError("Could not parse Zenodo h5ad shape")
    zenodo_cells = int(zenodo_shape.group(1).replace(",", ""))
    menstrual_cells = int(43054)
    coord_rows = parse_md_int(r"SNP evidence rows:\s*([0-9,]+)", coord_md)

    fig = plt.figure(figsize=(8.0, 5.5))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 1.0], left=0.17, right=0.96, top=0.87, bottom=0.11, hspace=0.64, wspace=0.46)
    ax_flow = fig.add_subplot(gs[0, :])
    nodes = [
        ("GWAS", "3 EUR files", PALETTE["genetics"]),
        ("Mapped\nrsIDs", f"{coord_rows:,} rows", "#B8C8EA"),
        ("Candidate\nuniverse", "102 records", "#DDE8F7"),
        ("Bulk", "4 datasets", PALETTE["bulk"]),
        ("scRNA", "3 resources", PALETTE["singlecell"]),
        ("Hypothesis\ntriage", "rank + safety", PALETTE["action"]),
    ]
    draw_flow(ax_flow, nodes)
    ax_flow.text(0.01, 0.97, "a", transform=ax_flow.transAxes, fontsize=10, fontweight="bold", va="top")
    ax_flow.text(0.5, 0.98, "Genetics-first ordered evidence gates", ha="center", va="top", transform=ax_flow.transAxes, fontsize=10)

    ax_gwas = fig.add_subplot(gs[1, 0])
    labels = ["Adenomyosis", "Endometriosis", "Endometriosis\nwithout adenomyosis"]
    values = [float(v) / 1e6 for v in gwas["rows"]]
    colors = [PALETTE["genetics"], "#7895D2", "#A5B6E3"]
    horizontal_bar(ax_gwas, labels, values, colors, "SNP rows (millions)")
    ax_gwas.set_title("GWAS inputs")
    add_panel_label(ax_gwas, "b")

    ax_bulk = fig.add_subplot(gs[1, 1])
    labels = bulk_meta["dataset"].astype(str).tolist()
    values = pd.to_numeric(bulk_meta["matched_samples"], errors="coerce").tolist()
    horizontal_bar(ax_bulk, labels, values, [PALETTE["bulk"]] * len(labels), "Matched samples")
    ax_bulk.set_title("Bulk transcriptome inputs")
    add_panel_label(ax_bulk, "c")

    ax_sc = fig.add_subplot(gs[2, 0])
    labels = ["GSE179640\nendo scRNA", "Zenodo 17078290\nadeno scRNA", "GSE203191\nmenstrual effluent"]
    values = [gse179640_cells / 1000, zenodo_cells / 1000, menstrual_cells / 1000]
    horizontal_bar(ax_sc, labels, values, [PALETTE["singlecell"], "#E4B965", "#ECCE90"], "Cells (thousands)")
    ax_sc.set_title("Single-cell resources")
    add_panel_label(ax_sc, "d")

    ax_guard = fig.add_subplot(gs[2, 1])
    ax_guard.set_axis_off()
    guard_items = [
        ("GWAS", "coordinate-window\ncandidates"),
        ("Bulk", "candidate-level\nmodels"),
        ("scRNA", "broad labels /\ndonor-aware support"),
        ("Adenomyosis\nh5ad", "scRNA only;\nno spatial coordinates"),
    ]
    y0 = 0.84
    for i, (left, right) in enumerate(guard_items):
        y = y0 - i * 0.19
        ax_guard.text(0.02, y, left, ha="left", va="center", fontsize=7.3, fontweight="bold", color=PALETTE["neutral_3"])
        ax_guard.text(0.46, y, right, ha="left", va="center", fontsize=6.9, color=PALETTE["neutral_3"], linespacing=0.9)
        ax_guard.plot([0.02, 0.95], [y - 0.08, y - 0.08], color="#E1E1E1", lw=0.6)
    ax_guard.text(0.0, 1.03, "e", transform=ax_guard.transAxes, fontsize=10, fontweight="bold", va="bottom")

    data_rows = []
    for _, row in gwas.iterrows():
        data_rows.append({"panel": "gwas_rows", "resource": row["phenotype"], "value": int(row["rows"]), "unit": "SNP rows"})
    for _, row in bulk_meta.iterrows():
        data_rows.append({"panel": "bulk_samples", "resource": row["dataset"], "value": int(row["matched_samples"]), "unit": "matched samples"})
    for resource, value in [
        ("GSE179640 endometriosis scRNA", gse179640_cells),
        ("Zenodo 17078290 adenomyosis scRNA", zenodo_cells),
        ("GSE203191 menstrual effluent", menstrual_cells),
    ]:
        data_rows.append({"panel": "singlecell_cells", "resource": resource, "value": int(value), "unit": "cells"})
    data_rows.append({"panel": "coordinate_mapping", "resource": "rsID evidence rows", "value": int(coord_rows), "unit": "rows"})
    data = pd.DataFrame(data_rows)
    data.to_csv(FIGURES / "figure1_study_design_audit_data.tsv", sep="\t", index=False)
    fig.suptitle("Audited public resources and evidence gates", fontsize=12, y=0.975)
    save_figure(fig, "figure1")


def make_figure2(s: dict[str, object]) -> None:
    gwas = s["gwas_signal"]
    neighborhoods = s["neighborhoods"]
    ld = s["ld"]
    universe = s["universe"]
    gw = neighborhoods[neighborhoods["threshold"] == "genome_wide"]
    ld_gw = ld[ld["threshold"] == "genome_wide"]

    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(2, 2, left=0.10, right=0.96, top=0.91, bottom=0.12, hspace=0.48, wspace=0.38)

    ax1 = fig.add_subplot(gs[0, 0])
    labels = ["Adenomyosis", "Endometriosis", "Endometriosis\nwithout adeno"]
    values = pd.to_numeric(gwas["genome_wide_snps_p_lt_5e_8"], errors="coerce").tolist()
    horizontal_bar(ax1, labels, values, [PALETTE["genetics"], "#7895D2", "#A5B6E3"], "SNPs at P < 5 x 10^-8")
    ax1.set_xscale("log")
    ax1.set_title("Genome-wide SNP signals")
    add_panel_label(ax1, "a")

    ax2 = fig.add_subplot(gs[0, 1])
    cats = ["All", "Shared by\n>=2 phenotypes", "All 3"]
    vals = [len(gw), int((gw["n_phenotypes"] >= 2).sum()), int((gw["n_phenotypes"] == 3).sum())]
    ax2.bar(cats, vals, color=[PALETTE["neutral_2"], PALETTE["genetics"], PALETTE["action"]], edgecolor="white")
    ax2.set_ylabel("Coordinate-window neighborhoods")
    ax2.set_title("Shared region construction")
    for i, v in enumerate(vals):
        ax2.text(i, v + max(vals) * 0.03, str(v), ha="center", va="bottom", fontsize=8)
    add_panel_label(ax2, "b")

    ax3 = fig.add_subplot(gs[1, 0])
    class_order = ["same_lead_snp_supported", "high_ld_supported", "moderate_ld_supported", "weak_or_distinct_ld"]
    vals = [int((ld_gw["ld_neighborhood_class"] == c).sum()) for c in class_order]
    labels = ["Same lead", "High LD", "Moderate LD", "Weak/distinct"]
    ax3.bar(labels, vals, color=[PALETTE["action"], PALETTE["genetics"], PALETTE["singlecell"], PALETTE["neutral_2"]], edgecolor="white")
    ax3.set_ylabel("Shared neighborhoods")
    ax3.tick_params(axis="x", rotation=25)
    ax3.set_title("LD sensitivity classes")
    for i, v in enumerate(vals):
        ax3.text(i, v + 0.08, str(v), ha="center", va="bottom", fontsize=8)
    add_panel_label(ax3, "c")

    ax4 = fig.add_subplot(gs[1, 1])
    priority_order = [
        "tier1_genome_wide_shared_ld_supported",
        "tier2_genome_wide_shared_moderate_ld",
        "tier3_genome_wide_shared_weak_or_unresolved_ld",
    ]
    counts = universe["genetic_priority"].value_counts().to_dict()
    vals = [counts.get(k, 0) for k in priority_order]
    labels = ["Tier 1\nLD-supported", "Tier 2\nmoderate LD", "Tier 3\nweak/unresolved"]
    ax4.bar(labels, vals, color=[PALETTE["genetics"], "#8CA4D9", PALETTE["neutral_1"]], edgecolor="white")
    ax4.set_ylabel("Candidate records")
    ax4.set_title("Genetics-derived universe")
    for i, v in enumerate(vals):
        ax4.text(i, v + 1.5, str(v), ha="center", va="bottom", fontsize=8)
    add_panel_label(ax4, "d")

    pd.DataFrame(
        {
            "metric": ["gw_snp_adenomyosis", "gw_snp_endometriosis", "gw_snp_endo_wo_adeno", "neighborhood_all", "neighborhood_shared2", "neighborhood_all3", "tier1", "tier2", "tier3"],
            "value": values + [len(gw), int((gw["n_phenotypes"] >= 2).sum()), int((gw["n_phenotypes"] == 3).sum())] + vals,
        }
    ).to_csv(FIGURES / "figure2_genetic_candidate_universe_data.tsv", sep="\t", index=False)
    fig.suptitle("Coordinate-window GWAS mapping defines a constrained candidate universe", fontsize=12)
    save_figure(fig, "figure2")


def make_figure3(s: dict[str, object]) -> None:
    models = s["bulk_models"].copy()
    scores = s["bulk_scores"].copy()
    selected = scores.sort_values("bulk_expression_support_score_20", ascending=False, kind="mergesort").head(15).copy()

    fig = plt.figure(figsize=(8.2, 6.0))
    gs = fig.add_gridspec(2, 2, left=0.31, right=0.96, top=0.90, bottom=0.13, hspace=0.55, wspace=0.38)

    ax1 = fig.add_subplot(gs[0, 0])
    model_keep = models[models["analysis"].str.contains("GSE141549|GSE51981|GSE313775|GSE234354", regex=True)].copy()
    label_map = {
        "GSE234354_cycle_stage_F_test": "GSE234354 cycle stage",
        "GSE313775_endometriosis_vs_control_Th1": "GSE313775 Th1 disease",
        "GSE313775_endometriosis_vs_control_Th1_17": "GSE313775 Th1/17 disease",
        "GSE313775_endometriosis_vs_control_Th17": "GSE313775 Th17 disease",
        "GSE313775_Th17_interaction_vs_Th1": "GSE313775 Th17 interaction",
        "GSE313775_Th1_17_interaction_vs_Th1": "GSE313775 Th1/17 interaction",
        "GSE141549_tissue_contrast_lesion-control_endometrium": "GSE141549 lesion vs control endometrium",
        "GSE141549_tissue_contrast_lesion-patient_eutopic_endometrium": "GSE141549 lesion vs eutopic endometrium",
        "GSE141549_tissue_contrast_patient_eutopic_endometrium-control_endometrium": "GSE141549 eutopic vs control endometrium",
        "GSE141549_tissue_contrast_patient_peritoneum-control_peritoneum": "GSE141549 peritoneum disease",
        "GSE141549_tissue_contrast_lesion-patient_peritoneum": "GSE141549 lesion vs peritoneum",
        "GSE51981_Endometriosis_vs_Non_Endometriosis_adjusted_cycle_all_mapped_probes": "GSE51981 disease, all probes",
        "GSE51981_Endometriosis_vs_Non_Endometriosis_adjusted_cycle_single_gene_probes": "GSE51981 disease, single-gene probes",
    }
    model_keep["short"] = model_keep["analysis"].map(label_map).fillna(model_keep["analysis"].astype(str))
    model_keep = model_keep.sort_values("fdr_lt_0_10", ascending=True).tail(9)
    horizontal_bar(
        ax1,
        model_keep["short"].tolist(),
        pd.to_numeric(model_keep["fdr_lt_0_10"], errors="coerce").tolist(),
        [PALETTE["bulk"]] * len(model_keep),
        "Candidates with FDR < 0.10",
    )
    ax1.set_title("Candidate-level bulk model support")
    add_panel_label(ax1, "a")

    ax2 = fig.add_subplot(gs[0, 1])
    class_order = ["high_bulk_support", "moderate_bulk_support", "limited_bulk_support"]
    counts = scores["bulk_support_class"].value_counts().to_dict()
    vals = [counts.get(k, 0) for k in class_order]
    ax2.bar(["High", "Moderate", "Limited"], vals, color=[PALETTE["action"], PALETTE["bulk"], PALETTE["neutral_1"]], edgecolor="white")
    ax2.set_ylabel("Genes")
    ax2.set_title("Bulk evidence classes")
    for i, v in enumerate(vals):
        ax2.text(i, v + 1.2, str(v), ha="center", fontsize=8)
    add_panel_label(ax2, "b")

    ax3 = fig.add_subplot(gs[1, 0])
    heat_cols = ["lesion_support_score_6", "independent_endometrium_validation_score_4", "immune_support_score_4", "cycle_nonconfounded_score_4"]
    heat = selected[heat_cols].to_numpy(dtype=float)
    denom = np.array([6, 4, 4, 4], dtype=float)
    im = ax3.imshow(heat / denom, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
    ax3.set_yticks(np.arange(len(selected)))
    ax3.set_yticklabels(selected["gene_symbol"].tolist(), fontsize=6.5)
    ax3.set_xticks(np.arange(len(heat_cols)))
    ax3.set_xticklabels(["Lesion", "GSE51981", "Immune", "Low cycle\nconfounding"], fontsize=6.5)
    ax3.tick_params(length=0)
    ax3.set_title("Top bulk-supported candidates")
    add_panel_label(ax3, "c")
    cbar = fig.colorbar(im, ax=ax3, fraction=0.046, pad=0.02)
    cbar.set_label("Scaled support", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    ax4 = fig.add_subplot(gs[1, 1])
    x = pd.to_numeric(selected["lesion_logFC"], errors="coerce")
    y = pd.to_numeric(selected["immune_logFC"], errors="coerce")
    cycle = selected["cycle_driven_flag"].astype(str) == "True"
    colors = np.where(cycle, PALETTE["warning"], PALETTE["bulk"])
    ax4.scatter(x, y, s=45 + 4 * pd.to_numeric(selected["bulk_expression_support_score_20"], errors="coerce"), c=colors, edgecolor="white", linewidth=0.5)
    for _, row in selected.head(10).iterrows():
        ax4.text(float(row["lesion_logFC"]) + 0.025, float(row["immune_logFC"]) + 0.035, str(row["gene_symbol"]), fontsize=6.1, ha="left", va="bottom")
    ax4.axhline(0, color="#CFCFCF", lw=0.7)
    ax4.axvline(0, color="#CFCFCF", lw=0.7)
    ax4.set_xlabel("Lesion logFC")
    ax4.set_ylabel("Immune-cell logFC")
    ax4.set_title("Direction and cycle sensitivity")
    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="", color=PALETTE["warning"], label="Cycle-sensitive"),
        Line2D([0], [0], marker="o", linestyle="", color=PALETTE["bulk"], label="Lower cycle flag"),
    ]
    ax4.legend(handles=legend_handles, fontsize=6.5, loc="best")
    add_panel_label(ax4, "d")

    selected.to_csv(FIGURES / "figure3_bulk_expression_support_data.tsv", sep="\t", index=False)
    fig.suptitle("Candidate-level bulk expression support across lesion, immune and cycle-control layers", fontsize=12)
    save_figure(fig, "figure3")


def make_figure4(s: dict[str, object]) -> None:
    gse179640_md = s["gse179640_md"]
    gse179640_v2_md = s["gse179640_v2_md"]
    zenodo_md = s["zenodo_md"]
    cross = s["cross"].copy()

    comp_counts = parse_literal_counts("Broad-compartment totals", gse179640_md)
    main_comps = {
        "Stromal/\nfibroblast": comp_counts.get("stromal_fibroblast", 0),
        "T/NK": comp_counts.get("t_nk", 0),
        "Epithelial": comp_counts.get("epithelial", 0),
        "Myeloid/\nmacrophage": comp_counts.get("myeloid_macrophage", 0),
        "Endothelial": comp_counts.get("endothelial", 0),
        "Mural/\nsmooth muscle": comp_counts.get("mural_smooth_muscle", 0),
    }
    v2_counts = parse_literal_counts("Integrated class counts", gse179640_v2_md)
    z_counts = parse_literal_counts("Localization class counts", zenodo_md)

    fig = plt.figure(figsize=(8.0, 6.0))
    gs = fig.add_gridspec(2, 2, left=0.17, right=0.96, top=0.90, bottom=0.12, hspace=0.50, wspace=0.42)
    ax1 = fig.add_subplot(gs[0, 0])
    horizontal_bar(ax1, list(main_comps.keys()), [v / 1000 for v in main_comps.values()], [PALETTE["singlecell"]] * len(main_comps), "Cells (thousands)")
    ax1.set_title("GSE179640 broad compartments")
    add_panel_label(ax1, "a")

    ax2 = fig.add_subplot(gs[0, 1])
    labels = ["Moderate\nrelaxed", "Moderate\ndonor-aware", "Suggestive", "Limited", "Minimal"]
    vals = [
        v2_counts.get("moderate_relaxed_fdr_singlecell_support", 0),
        v2_counts.get("moderate_donor_aware_singlecell_support", 0),
        v2_counts.get("suggestive_donor_aware_singlecell_support", 0),
        v2_counts.get("limited_singlecell_support", 0),
        v2_counts.get("minimal_singlecell_support", 0),
    ]
    ax2.bar(labels, vals, color=[PALETTE["singlecell"], "#E6BD72", "#F0D49A", PALETTE["neutral_1"], PALETTE["neutral_2"]], edgecolor="white")
    ax2.set_ylabel("Genes")
    ax2.set_title("Endometriosis scRNA support")
    ax2.tick_params(axis="x", rotation=20)
    for i, v in enumerate(vals):
        ax2.text(i, v + 0.8, str(v), ha="center", fontsize=8)
    add_panel_label(ax2, "b")

    ax3 = fig.add_subplot(gs[1, 0])
    labels = ["High", "Moderate", "Limited", "Minimal", "Unmatched"]
    vals = [
        z_counts.get("high_adenomyosis_scRNA_localization", 0),
        z_counts.get("moderate_adenomyosis_scRNA_localization", 0),
        z_counts.get("limited_adenomyosis_scRNA_localization", 0),
        z_counts.get("minimal_adenomyosis_scRNA_localization", 0),
        z_counts.get("not_detected_or_unmatched", 0),
    ]
    ax3.bar(labels, vals, color=[PALETTE["singlecell"], "#E6BD72", "#F0D49A", PALETTE["neutral_1"], PALETTE["neutral_2"]], edgecolor="white")
    ax3.set_ylabel("Genes")
    ax3.set_title("Adenomyosis h5ad cell-label context\nn=3/group; no spatial coordinates")
    ax3.tick_params(axis="x", rotation=20)
    for i, v in enumerate(vals):
        ax3.text(i, v + 1.0, str(v), ha="center", fontsize=8)
    add_panel_label(ax3, "c")

    ax4 = fig.add_subplot(gs[1, 1])
    selected_genes = ["HSPG2", "SSPN", "KDR", "KDM1A", "LY96", "PDGFRA", "ECE1", "C1QA", "KIT"]
    sub = cross[cross["gene_symbol"].isin(selected_genes)].copy()
    sub["axis_color"] = sub["dominant_cross_disease_axis"].map(
        {"fibrovascular": PALETTE["fibro"], "immune": PALETTE["immune"], "epithelial": PALETTE["epithelial"], "mixed": PALETTE["mixed"]}
    ).fillna(PALETTE["neutral_2"])
    ax4.scatter(
        pd.to_numeric(sub["endo_scRNA_normalized_score_20"], errors="coerce"),
        pd.to_numeric(sub["adeno_scRNA_normalized_score_20"], errors="coerce"),
        s=55 + pd.to_numeric(sub["pre_druggability_biologic_evidence_score_80"], errors="coerce"),
        c=sub["axis_color"],
        edgecolor="white",
        linewidth=0.6,
    )
    callout_positions = {
        "KDR": (20.15, 19.35),
        "ECE1": (20.35, 18.95),
        "SSPN": (20.35, 18.55),
        "LY96": (20.35, 18.15),
        "KDM1A": (20.10, 17.55),
    }
    label_offsets = {
        "HSPG2": (0.16, 0.12),
        "PDGFRA": (0.16, 0.10),
        "C1QA": (0.16, 0.13),
        "KIT": (0.16, -0.05),
    }
    for _, row in sub.iterrows():
        x0 = float(row["endo_scRNA_normalized_score_20"])
        y0 = float(row["adeno_scRNA_normalized_score_20"])
        gene = str(row["gene_symbol"])
        if gene in callout_positions:
            ax4.annotate(
                gene,
                xy=(x0, y0),
                xytext=callout_positions[gene],
                textcoords="data",
                fontsize=6.2,
                arrowprops=dict(arrowstyle="-", color="#9A9A9A", lw=0.5, shrinkA=2, shrinkB=2),
            )
        else:
            dx, dy = label_offsets.get(gene, (0.12, 0.08))
            ax4.text(x0 + dx, y0 + dy, gene, fontsize=6.4)
    ax4.set_xlim(10, 22.2)
    ax4.set_ylim(10, 20.8)
    ax4.set_xlabel("Endometriosis scRNA context score")
    ax4.set_ylabel("Adenomyosis cell-label context score")
    ax4.set_title("Cross-disease candidate context")
    ax4.grid(color="#E4E4E4", lw=0.5)
    add_panel_label(ax4, "d")

    figure4_rows = []
    for label, value in main_comps.items():
        figure4_rows.append({"panel": "gse179640_compartment", "label": label, "value": int(value)})
    for label, value in zip(
        ["moderate_relaxed", "moderate_donor", "suggestive", "limited", "minimal"],
        [
            v2_counts.get("moderate_relaxed_fdr_singlecell_support", 0),
            v2_counts.get("moderate_donor_aware_singlecell_support", 0),
            v2_counts.get("suggestive_donor_aware_singlecell_support", 0),
            v2_counts.get("limited_singlecell_support", 0),
            v2_counts.get("minimal_singlecell_support", 0),
        ],
    ):
        figure4_rows.append({"panel": "endo_scRNA_class", "label": label, "value": int(value)})
    for label, value in zip(
        ["high", "moderate", "limited", "minimal", "unmatched"],
        [
            z_counts.get("high_adenomyosis_scRNA_localization", 0),
            z_counts.get("moderate_adenomyosis_scRNA_localization", 0),
            z_counts.get("limited_adenomyosis_scRNA_localization", 0),
            z_counts.get("minimal_adenomyosis_scRNA_localization", 0),
            z_counts.get("not_detected_or_unmatched", 0),
        ],
    ):
        figure4_rows.append({"panel": "adeno_scRNA_class", "label": label, "value": int(value)})
    pd.DataFrame(figure4_rows).to_csv(FIGURES / "figure4_singlecell_celllabel_context_data.tsv", sep="\t", index=False)
    sub[
        [
            "gene_symbol",
            "dominant_cross_disease_axis",
            "endo_scRNA_normalized_score_20",
            "adeno_scRNA_normalized_score_20",
            "pre_druggability_biologic_evidence_score_80",
        ]
    ].to_csv(FIGURES / "figure4_singlecell_candidate_scores.tsv", sep="\t", index=False)
    fig.suptitle("Broad-compartment and h5ad cell-label context of candidate genes", fontsize=12)
    save_figure(fig, "figure4")


def make_figure5(s: dict[str, object]) -> None:
    cross = s["cross"].copy()
    cross_md = s["cross_md"]
    class_counts = parse_literal_counts("Cross-disease scRNA classes", cross_md)
    axis_counts = parse_literal_counts("Dominant axis counts", cross_md)
    top = cross.sort_values("pre_druggability_biologic_evidence_score_80", ascending=False).head(16).copy()

    fig = plt.figure(figsize=(7.2, 6.0))
    gs = fig.add_gridspec(2, 2, left=0.10, right=0.97, top=0.91, bottom=0.13, hspace=0.50, wspace=0.38)

    ax1 = fig.add_subplot(gs[0, 0])
    class_order = [
        "shared_scRNA_localization",
        "adenomyosis_dominant_with_endometriosis_signal",
        "endometriosis_dominant_with_adenomyosis_signal",
        "limited_or_single_layer_scRNA_localization",
        "minimal_scRNA_localization",
    ]
    labels = ["Shared", "Adeno-\ndominant", "Endo-\ndominant", "Limited/\nsingle", "Minimal"]
    vals = [class_counts.get(k, 0) for k in class_order]
    ax1.bar(labels, vals, color=[PALETTE["genetics"], "#87A1DD", "#C4A2C2", PALETTE["neutral_1"], PALETTE["neutral_2"]], edgecolor="white")
    ax1.set_ylabel("Genes")
    ax1.set_title("Cross-disease scRNA context classes")
    for i, v in enumerate(vals):
        ax1.text(i, v + 1.0, str(v), ha="center", fontsize=8)
    add_panel_label(ax1, "a")

    ax2 = fig.add_subplot(gs[0, 1])
    axis_order = ["fibrovascular", "epithelial", "mixed", "immune", "none_detected"]
    vals = [axis_counts.get(k, 0) for k in axis_order]
    colors = [PALETTE["fibro"], PALETTE["epithelial"], PALETTE["mixed"], PALETTE["immune"], PALETTE["neutral_1"]]
    ax2.bar(["Fibro-\nvascular", "Epithelial", "Mixed", "Immune", "None"], vals, color=colors, edgecolor="white")
    ax2.set_ylabel("Genes")
    ax2.set_title("Dominant context axes")
    for i, v in enumerate(vals):
        ax2.text(i, v + 0.8, str(v), ha="center", fontsize=8)
    add_panel_label(ax2, "b")

    ax3 = fig.add_subplot(gs[1, 0])
    matrix = np.column_stack(
        [
            pd.to_numeric(top["genetics_support_score_20"], errors="coerce") / 20,
            pd.to_numeric(top["bulk_expression_support_score_20"], errors="coerce") / 20,
            pd.to_numeric(top["cross_disease_scRNA_score_40"], errors="coerce") / 40,
            pd.to_numeric(top["pre_druggability_biologic_evidence_score_80"], errors="coerce") / 80,
        ]
    )
    im = ax3.imshow(matrix, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
    ax3.set_yticks(np.arange(len(top)))
    ax3.set_yticklabels(top["gene_symbol"].tolist(), fontsize=6.2)
    ax3.set_xticks(np.arange(4))
    ax3.set_xticklabels(["Genetics", "Bulk", "scRNA", "Biological\nevidence"], fontsize=6.5)
    ax3.tick_params(length=0)
    ax3.set_title("Top biological evidence genes")
    add_panel_label(ax3, "c")
    cbar = fig.colorbar(im, ax=ax3, fraction=0.046, pad=0.02)
    cbar.set_label("Scaled evidence", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    ax4 = fig.add_subplot(gs[1, 1])
    flag_cols = ["shared_fibrovascular_flag", "shared_immune_flag", "shared_epithelial_flag"]
    flag_matrix = top[flag_cols].astype(bool).astype(int).to_numpy()
    cmap_colors = np.array([[0.93, 0.93, 0.93, 1.0], [0.30, 0.49, 0.66, 1.0]])
    ax4.imshow(flag_matrix, aspect="auto", cmap=plt.matplotlib.colors.ListedColormap(cmap_colors), vmin=0, vmax=1)
    ax4.set_yticks(np.arange(len(top)))
    ax4.set_yticklabels(top["gene_symbol"].tolist(), fontsize=6.2)
    ax4.set_xticks(np.arange(3))
    ax4.set_xticklabels(["Fibro-\nvascular", "Immune", "Epithelial"], fontsize=6.5)
    ax4.tick_params(length=0)
    ax4.set_title("Shared-axis flags")
    ax4.text(
        0.5,
        -0.18,
        "Filled cells indicate assigned axis support",
        transform=ax4.transAxes,
        ha="center",
        va="top",
        fontsize=6.4,
        color=PALETTE["neutral_3"],
    )
    add_panel_label(ax4, "d")

    top.to_csv(FIGURES / "figure5_cross_disease_convergence_data.tsv", sep="\t", index=False)
    fig.suptitle("Cross-disease cell-context overlap before actionability scoring", fontsize=12)
    save_figure(fig, "figure5")


def write_captions() -> None:
    text = """# Draft captions for Figures 1-5

**Figure 1. Audited public-resource workflow for genetics-guided experimental-hypothesis prioritisation.** The analysis began with three EUR GWAS summary-statistic files and proceeded through coordinate-window candidate construction, candidate-level bulk expression support, single-cell or cell-label context evidence and model-system hypothesis triage. The adenomyosis h5ad was used as scRNA/cell-label evidence only because no usable spatial coordinates were detected; it is not spatial evidence.

**Figure 2. Cross-phenotype GWAS mapping defines the candidate universe.** SNP-level signals were converted into coordinate-window neighbourhoods and then stratified by pairwise EUR LD sensitivity. The resulting 102 Ensembl gene records are coordinate-neighbourhood candidates, not locus-resolved causal genes or therapeutic targets.

**Figure 3. Candidate-level bulk expression support across lesion, immune and cycle-control layers.** Bulk models were restricted to the genetics-derived candidate universe. GSE234354 was used as a menstrual-cycle confounding-control layer rather than disease evidence.

**Figure 4. Broad-compartment and h5ad cell-label context evidence across endometriosis and adenomyosis resources.** GSE179640 supports broad-compartment context evidence in endometriosis, while Zenodo 17078290 supports adenomyosis scRNA/cell-label context evidence. The adenomyosis h5ad layer has n=3 samples per group and no usable spatial coordinates in the local audit. These panels do not claim final cluster-level differential expression, and the adenomyosis h5ad is not spatial evidence.

**Figure 5. Cross-resource cell-context overlap before targetability filtering.** The integrated cell-context matrix identifies broad fibrovascular, epithelial and immune candidate patterns that feed into Figure 6 scoring, where targetability, safety and directionality are added.
"""
    CAPTIONS.write_text(text, encoding="utf-8")


def image_nonblank_fraction(path: Path) -> float:
    arr = mpimg.imread(path)
    if arr.ndim == 2:
        rgb = arr
    else:
        rgb = arr[..., :3]
    nonwhite = np.any(rgb < 0.985, axis=-1)
    return float(nonwhite.mean())


def write_self_review() -> tuple[str, list[str]]:
    failures: list[str] = []
    for name, (svg, png) in FIGURE_PATHS.items():
        if not svg.exists() or svg.stat().st_size < 10_000:
            failures.append(f"{name} SVG missing or too small: {svg}")
        if not png.exists() or png.stat().st_size < 10_000:
            failures.append(f"{name} PNG missing or too small: {png}")
        if svg.exists():
            svg_text = svg.read_text(encoding="utf-8", errors="ignore")
            if "<text" not in svg_text:
                failures.append(f"{name} SVG lacks editable text nodes.")
            bad_spatial = "spatial evidence" in svg_text.lower() and "not spatial evidence" not in svg_text.lower()
            if bad_spatial:
                failures.append(f"{name} may imply spatial evidence without the not-spatial guardrail.")
        if png.exists():
            frac = image_nonblank_fraction(png)
            if frac < 0.02:
                failures.append(f"{name} PNG appears nearly blank; non-white fraction={frac:.4f}.")
    if not CAPTIONS.exists() or "not spatial evidence" not in CAPTIONS.read_text(encoding="utf-8").lower():
        failures.append("Captions missing or lack adenomyosis not-spatial guardrail.")
    required_data = [
        FIGURES / "figure1_study_design_audit_data.tsv",
        FIGURES / "figure2_genetic_candidate_universe_data.tsv",
        FIGURES / "figure3_bulk_expression_support_data.tsv",
        FIGURES / "figure4_singlecell_celllabel_context_data.tsv",
        FIGURES / "figure4_singlecell_candidate_scores.tsv",
        FIGURES / "figure5_cross_disease_convergence_data.tsv",
    ]
    for path in required_data:
        if not path.exists() or path.stat().st_size == 0:
            failures.append(f"Missing figure data table: {path}")

    status = "PASS" if not failures else "FAIL"
    lines = [
        "# Phase 18 self-review: Figures 1-5 draft panels",
        "",
        f"Status: {status}",
        "",
        "## Checks",
        "",
    ]
    for name, (svg, png) in FIGURE_PATHS.items():
        svg_size = svg.stat().st_size if svg.exists() else 0
        png_size = png.stat().st_size if png.exists() else 0
        nonblank = image_nonblank_fraction(png) if png.exists() else 0
        lines.append(f"- {name}: SVG {svg_size:,} bytes; PNG {png_size:,} bytes; non-white fraction {nonblank:.4f}.")
    lines.extend(
        [
            "",
            "## Scientific guardrails",
            "",
            "- Figure 1 separates evidence resources from interpretation gates.",
            "- Figure 2 labels outputs as coordinate-window candidates rather than causal loci.",
            "- Figure 3 is candidate-level bulk expression support rather than a genome-wide DEG figure.",
            "- Figure 4 explicitly states adenomyosis h5ad is scRNA/cell-label evidence, not spatial evidence.",
            "- Figure 5 feeds cell-context overlap into Figure 6 and does not present intervention candidates before targetability review.",
            "",
        ]
    )
    if failures:
        lines.extend(["## Failures", ""])
        lines.extend(f"- {failure}" for failure in failures)
    else:
        lines.extend(["## Decision", "", "The draft Figures 1-5 pass automated output, readability and overclaim checks."])
    SELF_REVIEW.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return status, failures


def main() -> int:
    FIGURES.mkdir(parents=True, exist_ok=True)
    INTEGRATION.mkdir(parents=True, exist_ok=True)
    apply_publication_style()
    sources = load_sources()
    make_figure1(sources)
    make_figure2(sources)
    make_figure3(sources)
    make_figure4(sources)
    make_figure5(sources)
    write_captions()
    status, failures = write_self_review()
    for name, (svg, png) in FIGURE_PATHS.items():
        print(f"Wrote {svg}")
        print(f"Wrote {png}")
    print(f"Wrote {CAPTIONS}")
    print(f"Wrote {SELF_REVIEW}")
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
