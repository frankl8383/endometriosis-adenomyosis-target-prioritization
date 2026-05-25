#!/usr/bin/env python3
"""Build enhanced reviewer-facing sensitivity and boundary tables."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
INTEGRATION = RESULTS / "integration"
GWAS = RESULTS / "gwas"
BULK = RESULTS / "bulk"
SC = RESULTS / "singlecell"
DRUG = RESULTS / "druggability"
TABLES = RESULTS / "tables"

OUT_GWAS_SYMBOL_AUDIT = INTEGRATION / "gwas_candidate_gene_symbol_audit.tsv"
OUT_WINDOW_SENSITIVITY = INTEGRATION / "gwas_window_size_sensitivity.tsv"
OUT_FINAL_EVIDENCE = INTEGRATION / "final_candidate_evidence_summary.tsv"
OUT_TARGETABILITY_SPLIT = INTEGRATION / "targetability_safety_split.tsv"
OUT_EXTENDED_SENSITIVITY = INTEGRATION / "priority_sensitivity_extended.tsv"
OUT_INTERNAL_MATCHED_NULL = INTEGRATION / "internal_matched_candidate_null.tsv"
OUT_HASHES = INTEGRATION / "expected_outputs_sha256.tsv"
OUT_REVIEW = INTEGRATION / "phase30_submission_enhancement_self_review.md"

FINAL_AND_CONTEXT_GENES = ["KDM1A", "PDGFRA", "LY96", "KDR", "KIT", "ECE1", "C1QA", "HSPG2", "SSPN"]


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def soften_export_text(df: pd.DataFrame) -> pd.DataFrame:
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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_gwas_symbol_audit(universe: pd.DataFrame) -> pd.DataFrame:
    df = universe.copy()
    df["has_exported_gene_symbol"] = df["gene_symbol"].astype(str).str.strip().ne("")
    df["gene_record_label"] = np.where(df["has_exported_gene_symbol"], df["gene_symbol"], df["gene_id"])
    df["interpretation"] = np.where(
        df["has_exported_gene_symbol"],
        "named protein-coding gene record",
        "protein-coding Ensembl gene record without exported HGNC symbol in local annotation",
    )
    cols = [
        "gene_record_label",
        "gene_symbol",
        "gene_id",
        "gene_biotype",
        "genetic_priority",
        "ld_neighborhood_class",
        "best_lead_snp",
        "best_lead_p",
        "distance_to_best_lead_bp",
        "has_exported_gene_symbol",
        "interpretation",
    ]
    out = df[[c for c in cols if c in df.columns]].copy()
    out.to_csv(OUT_GWAS_SYMBOL_AUDIT, sep="\t", index=False)
    return out


def build_window_sensitivity(universe: pd.DataFrame) -> pd.DataFrame:
    df = universe.copy()
    df["gene_record_label"] = np.where(df["gene_symbol"].astype(str).str.strip().ne(""), df["gene_symbol"], df["gene_id"])
    distances = pd.to_numeric(df["distance_to_best_lead_bp"], errors="coerce").abs()
    windows = [250_000, 500_000, 1_000_000, 2_000_000]
    rows = []
    for _, row in df.iterrows():
        distance = abs(float(row.get("distance_to_best_lead_bp", np.nan))) if str(row.get("distance_to_best_lead_bp", "")) else np.nan
        for window in windows:
            rows.append(
                {
                    "gene_record_label": row["gene_record_label"],
                    "gene_symbol": row.get("gene_symbol", ""),
                    "gene_id": row.get("gene_id", ""),
                    "window_bp": window,
                    "distance_to_best_lead_bp": distance,
                    "retained_by_distance_to_best_lead": bool(pd.notna(distance) and distance <= window),
                    "genetic_priority": row.get("genetic_priority", ""),
                    "ld_neighborhood_class": row.get("ld_neighborhood_class", ""),
                    "best_lead_snp": row.get("best_lead_snp", ""),
                    "best_phenotype": row.get("best_phenotype", ""),
                    "interpretation": "Distance-to-best-lead sensitivity only; this is not LD clumping or fine mapping.",
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(OUT_WINDOW_SENSITIVITY, sep="\t", index=False)
    return out


def build_final_evidence_summary(
    universe: pd.DataFrame,
    bulk: pd.DataFrame,
    cross: pd.DataFrame,
    action: pd.DataFrame,
    rank: pd.DataFrame,
    manual: pd.DataFrame,
) -> pd.DataFrame:
    base = pd.DataFrame({"gene_symbol": FINAL_AND_CONTEXT_GENES})
    for df, cols in [
        (universe, ["gene_symbol", "gene_id", "genetic_priority", "ld_neighborhood_class", "best_lead_snp", "best_lead_p", "distance_to_best_lead_bp"]),
        (bulk, ["gene_symbol", "bulk_expression_support_score_20", "bulk_support_class", "cycle_adj_p", "cycle_driven_flag", "cycle_nonconfounded_score_4", "lesion_top_contrast", "lesion_logFC", "lesion_adj_p", "immune_top_contrast", "immune_logFC", "immune_adj_p"]),
        (cross, ["gene_symbol", "cross_disease_scRNA_score_40", "cross_disease_scRNA_class", "dominant_cross_disease_axis", "endo_top_broad_compartment", "top_adeno_cluster", "adeno_guardrail", "cross_disease_guardrail"]),
        (action, ["gene_symbol", "actionability_class", "druggability_score_20", "drug_direction_evidence", "druggability_penalty_reasons"]),
        (rank, ["gene_symbol", "final_target_priority_score_100_pre_rank_stability", "rank_stability_class", "bootstrap_top10_frequency", "leave_one_layer_top20_count", "gene_label_permutation_p_ge_observed", "rank_without_druggability", "rank_without_singlecell"]),
        (manual, ["gene_symbol", "manual_evidence_tier", "manual_category", "claim_strength", "manual_safety", "manual_rationale"]),
    ]:
        keep = [c for c in cols if c in df.columns]
        base = base.merge(df[keep], on="gene_symbol", how="left")
    base["cycle_boundary"] = np.where(
        base["cycle_driven_flag"].astype(str).str.lower().eq("true"),
        "cycle-sensitive; expression treated as context/directionality support, not disease-specific validation",
        "not flagged as cycle-driven in reference model",
    )
    base["claim_boundary"] = base["manual_evidence_tier"].map(
        {
            "A_primary": "experimental model-system hypothesis only",
            "B_primary_safety_limited": "safety-limited fibrovascular signal; not systemic treatment recommendation",
            "C_secondary": "secondary/exploratory model-system hypothesis",
            "D_context_or_anchor": "mechanistic anchor, not direct intervention candidate",
            "E_marker_not_target": "lesion-context marker, not direct intervention candidate",
        }
    ).fillna("context/exploratory candidate")
    base = soften_export_text(base)
    base.to_csv(OUT_FINAL_EVIDENCE, sep="\t", index=False)
    return base


def build_targetability_split(prelim: pd.DataFrame) -> pd.DataFrame:
    df = prelim.copy()
    for col in [
        "target_class_score_5",
        "drug_precedent_score_5",
        "evidence_convergence_score_4",
        "safety_acceptability_score_4",
        "local_repurposing_feasibility_score_2",
        "druggability_penalty_points",
    ]:
        df[col] = pd.to_numeric(df.get(col, 0), errors="coerce").fillna(0)
    df["technical_targetability_score_14"] = (
        df["target_class_score_5"] + df["drug_precedent_score_5"] + df["evidence_convergence_score_4"]
    )
    df["disease_context_safety_plausibility_score_6_minus_penalty"] = (
        df["safety_acceptability_score_4"] + df["local_repurposing_feasibility_score_2"] - df["druggability_penalty_points"]
    )
    df["split_interpretation"] = "Technical druggability is separated from disease-context safety/plausibility."
    keep = [
        "gene_symbol",
        "gene_id",
        "technical_targetability_score_14",
        "disease_context_safety_plausibility_score_6_minus_penalty",
        "actionability_class",
        "drug_direction_evidence",
        "druggability_penalty_reasons",
        "split_interpretation",
    ]
    out = df[keep].copy()
    out = soften_export_text(out)
    out.to_csv(OUT_TARGETABILITY_SPLIT, sep="\t", index=False)
    return out


def build_extended_sensitivity(prelim: pd.DataFrame, rank: pd.DataFrame) -> pd.DataFrame:
    df = prelim.copy()
    for col in [
        "genetics_support_score_20",
        "bulk_expression_support_score_20",
        "cross_disease_scRNA_score_40",
        "druggability_score_20",
        "final_target_priority_score_100_pre_rank_stability",
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
    df["score_without_targetability_rescaled_100"] = (
        (df["genetics_support_score_20"] + df["bulk_expression_support_score_20"] + df["cross_disease_scRNA_score_40"]) / 80 * 100
    )
    df["score_without_expression_rescaled_100"] = (
        (df["genetics_support_score_20"] + df["cross_disease_scRNA_score_40"] + df["druggability_score_20"]) / 80 * 100
    )
    df["score_without_cell_context_rescaled_100"] = (
        (df["genetics_support_score_20"] + df["bulk_expression_support_score_20"] + df["druggability_score_20"]) / 60 * 100
    )
    df["rank_without_targetability_recomputed"] = df["score_without_targetability_rescaled_100"].rank(ascending=False, method="min").astype(int)
    df["rank_without_expression_recomputed"] = df["score_without_expression_rescaled_100"].rank(ascending=False, method="min").astype(int)
    df["rank_without_cell_context_recomputed"] = df["score_without_cell_context_rescaled_100"].rank(ascending=False, method="min").astype(int)
    merged = df.merge(
        rank[["gene_symbol", "rank_without_druggability", "rank_without_bulk", "rank_without_singlecell"]],
        on="gene_symbol",
        how="left",
        suffixes=("", "_rank_table"),
    )
    merged["manual_literature_removed_status"] = (
        "No formal literature-free score was calculated because targeted literature was used only after computational ranking for manual classification."
    )
    keep = [
        "gene_symbol",
        "gene_id",
        "final_target_priority_score_100_pre_rank_stability",
        "pre_rank_stability_rank",
        "score_without_targetability_rescaled_100",
        "rank_without_targetability_recomputed",
        "rank_without_druggability",
        "score_without_expression_rescaled_100",
        "rank_without_expression_recomputed",
        "rank_without_bulk",
        "score_without_cell_context_rescaled_100",
        "rank_without_cell_context_recomputed",
        "rank_without_singlecell",
        "manual_literature_removed_status",
    ]
    out = merged[[c for c in keep if c in merged.columns]].copy()
    out = soften_export_text(out)
    out.to_csv(OUT_EXTENDED_SENSITIVITY, sep="\t", index=False)
    return out


def build_internal_matched_null(prelim: pd.DataFrame) -> pd.DataFrame:
    df = prelim.copy()
    score_col = "final_target_priority_score_100_pre_rank_stability"
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce").fillna(0)
    rows = []
    for gene in FINAL_AND_CONTEXT_GENES:
        row = df[df["gene_symbol"] == gene]
        if row.empty:
            continue
        r = row.iloc[0]
        matched = df[
            (df["bulk_support_class"].astype(str) == str(r["bulk_support_class"]))
            & (df["actionability_class"].astype(str) == str(r["actionability_class"]))
        ].copy()
        if matched.empty:
            matched = df.copy()
        rows.append(
            {
                "gene_symbol": gene,
                "observed_score": float(r[score_col]),
                "matched_by": "bulk_support_class + actionability_class",
                "matched_candidate_count": int(len(matched)),
                "matched_candidates_with_score_ge_observed": int((matched[score_col] >= float(r[score_col])).sum()),
                "internal_matched_empirical_p_ge": float(((matched[score_col] >= float(r[score_col])).sum()) / len(matched)),
                "interpretation": "Internal 102-candidate matched check; not a genome-wide matched-background permutation.",
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(OUT_INTERNAL_MATCHED_NULL, sep="\t", index=False)
    return out


def build_hashes(paths: list[Path]) -> pd.DataFrame:
    rows = []
    for path in paths:
        if path.exists():
            rows.append({"file": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": sha256(path)})
    out = pd.DataFrame(rows)
    out.to_csv(OUT_HASHES, sep="\t", index=False)
    return out


def write_review(outputs: dict[str, pd.DataFrame]) -> None:
    failures: list[str] = []
    symbol = outputs["symbol"]
    final = outputs["final"]
    if len(symbol) != 102:
        failures.append(f"Expected 102 candidate records in symbol audit, observed {len(symbol)}.")
    if int((~symbol["has_exported_gene_symbol"]).sum()) != 10:
        failures.append("Expected 10 records without exported gene symbols.")
    if not final["cycle_boundary"].astype(str).str.contains("cycle-sensitive").all():
        failures.append("Expected all final/context genes to carry cycle-sensitive boundary text.")
    for path in [
        OUT_GWAS_SYMBOL_AUDIT,
        OUT_WINDOW_SENSITIVITY,
        OUT_FINAL_EVIDENCE,
        OUT_TARGETABILITY_SPLIT,
        OUT_EXTENDED_SENSITIVITY,
        OUT_INTERNAL_MATCHED_NULL,
        OUT_HASHES,
    ]:
        if not path.exists() or path.stat().st_size == 0:
            failures.append(f"Missing or empty output: {path}")
    status = "PASS" if not failures else "FAIL"
    lines = [
        "# Phase 30 self-review: submission enhancement tables",
        "",
        f"Status: {status}",
        "",
        "## Outputs",
        "",
        f"- `{OUT_GWAS_SYMBOL_AUDIT.relative_to(ROOT)}`",
        f"- `{OUT_WINDOW_SENSITIVITY.relative_to(ROOT)}`",
        f"- `{OUT_FINAL_EVIDENCE.relative_to(ROOT)}`",
        f"- `{OUT_TARGETABILITY_SPLIT.relative_to(ROOT)}`",
        f"- `{OUT_EXTENDED_SENSITIVITY.relative_to(ROOT)}`",
        f"- `{OUT_INTERNAL_MATCHED_NULL.relative_to(ROOT)}`",
        f"- `{OUT_HASHES.relative_to(ROOT)}`",
        "",
        "## Guardrails",
        "",
        "- Window sensitivity is distance-to-best-lead only and is not described as LD clumping or fine mapping.",
        "- Internal matched null checks use only the 102-candidate universe and are not genome-wide background tests.",
        "- Literature-free scoring is explicitly not recalculated because targeted literature was used after computational ranking.",
        "- Cycle sensitivity is foregrounded for all final/context genes.",
        "",
    ]
    if failures:
        lines.extend(["## Failures", "", *[f"- {f}" for f in failures]])
    else:
        lines.extend(["## Decision", "", "Enhanced supplementary tables are ready for package assembly."])
    OUT_REVIEW.write_text("\n".join(lines) + "\n", encoding="utf-8")
    if failures:
        raise SystemExit(1)


def main() -> int:
    INTEGRATION.mkdir(parents=True, exist_ok=True)
    universe = read_tsv(GWAS / "gwas_candidate_gene_universe.tsv")
    bulk = read_tsv(BULK / "bulk_candidate_expression_support_scores.tsv")
    cross = read_tsv(INTEGRATION / "cross_disease_singlecell_localization_matrix.tsv")
    action = read_tsv(DRUG / "target_actionability_scores.tsv")
    rank = read_tsv(INTEGRATION / "rank_stability_combined_matrix.tsv")
    manual = read_tsv(INTEGRATION / "phase14_manual_target_directionality_review.tsv")
    prelim = read_tsv(INTEGRATION / "preliminary_target_priority_matrix.tsv")

    outputs = {
        "symbol": build_gwas_symbol_audit(universe),
        "window": build_window_sensitivity(universe),
        "final": build_final_evidence_summary(universe, bulk, cross, action, rank, manual),
        "split": build_targetability_split(prelim),
        "extended": build_extended_sensitivity(prelim, rank),
        "matched": build_internal_matched_null(prelim),
    }
    build_hashes(
        [
            OUT_GWAS_SYMBOL_AUDIT,
            OUT_WINDOW_SENSITIVITY,
            OUT_FINAL_EVIDENCE,
            OUT_TARGETABILITY_SPLIT,
            OUT_EXTENDED_SENSITIVITY,
            OUT_INTERNAL_MATCHED_NULL,
            RESULTS / "figures" / "figure1_study_design_audit.svg",
            RESULTS / "figures" / "figure6_experimental_hypothesis_triage.svg",
            TABLES / "table_final_target_shortlist.tsv",
        ]
    )
    outputs["hashes"] = pd.read_csv(OUT_HASHES, sep="\t")
    write_review(outputs)
    print(OUT_REVIEW)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
