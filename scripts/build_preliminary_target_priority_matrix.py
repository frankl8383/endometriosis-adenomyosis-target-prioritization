#!/usr/bin/env python3
"""Build a manuscript-facing preliminary target priority matrix."""

from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT_ROOT / "results"
DRUG = RESULTS / "druggability"
INTEGRATION = RESULTS / "integration"
PY_PKG_DIR = PROJECT_ROOT / "tools" / "py_packages"

if sys.version_info[:2] != (3, 12):
    raise SystemExit(
        "Use the bundled Python 3.12 runtime: "
        "python3 "
        "scripts/build_preliminary_target_priority_matrix.py"
    )

sys.path.insert(0, str(PY_PKG_DIR))

import pandas as pd  # noqa: E402


INPUT = DRUG / "target_actionability_scores.tsv"
OUT_MATRIX = INTEGRATION / "preliminary_target_priority_matrix.tsv"
OUT_SUMMARY = INTEGRATION / "preliminary_target_priority_summary.md"
OUT_REVIEW = INTEGRATION / "phase12_preliminary_target_priority_self_review.md"


DIRECT_EXCLUSION_PENALTIES = [
    "structural_or_ecm_readout_without_direct_drug_precedent",
    "no_external_actionability_record",
]


def is_truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def priority_tier(row: pd.Series) -> str:
    total = float(row["final_target_priority_score_100_pre_rank_stability"])
    drug = float(row["druggability_score_20"])
    penalties = str(row.get("druggability_penalty_reasons", ""))
    actionability = str(row.get("actionability_class", ""))

    if any(token in penalties for token in DIRECT_EXCLUSION_PENALTIES):
        if total >= 70:
            return "tier_1_biologic_marker_or_context_gene_not_direct_target"
        return "tier_3_context_gene_not_direct_target"
    if actionability == "minimal_actionability":
        if total >= 70:
            return "tier_1_biologic_marker_or_context_gene_not_direct_target"
        return "tier_3_context_gene_not_direct_target"
    if total >= 75 and drug >= 12:
        return "tier_1_reviewable_target_hypothesis"
    if total >= 65 and drug >= 8:
        return "tier_2_supporting_target_hypothesis"
    if total >= 55:
        return "tier_3_biologic_context_gene"
    return "tier_4_low_priority_for_current_manuscript"


def manuscript_role(row: pd.Series) -> str:
    tier = str(row["priority_tier"])
    cross_class = str(row.get("cross_disease_scRNA_class", ""))
    actionability = str(row.get("actionability_class", ""))
    axis = str(row.get("dominant_cross_disease_axis", ""))
    if "not_direct_target" in tier:
        return "mechanistic_marker_or_niche_context"
    if actionability == "minimal_actionability":
        return "supporting_biologic_context"
    if "shared" in cross_class:
        return f"shared_{axis}_candidate_target"
    if "endometriosis" in cross_class:
        return f"endometriosis_prioritized_{axis}_candidate"
    if "adenomyosis" in cross_class:
        return f"adenomyosis_prioritized_{axis}_candidate"
    return f"single_layer_{axis}_candidate"


def review_flags(row: pd.Series) -> str:
    flags: list[str] = []
    if not str(row.get("gene_symbol", "")).strip():
        flags.append("missing_hgnc_symbol")
    penalties = str(row.get("druggability_penalty_reasons", ""))
    if penalties:
        flags.extend(penalties.split("|"))
    if "multiple_safety_liabilities" in penalties or "systemic_modulation_safety_caution" in penalties:
        flags.append("requires_safety_section_attention")
    if str(row.get("actionability_class", "")) == "minimal_actionability":
        flags.append("poor_drug_precedent")
    if float(row.get("endo_pseudobulk_score_10", 0) or 0) < 5:
        flags.append("weak_endometriosis_donor_aware_de")
    if str(row.get("adeno_guardrail", "")).strip():
        flags.append(str(row.get("adeno_guardrail", "")))
    return "|".join(dict.fromkeys([flag for flag in flags if flag]))


def suggested_shortlist_label(row: pd.Series) -> str:
    tier = str(row["priority_tier"])
    penalties = str(row.get("druggability_penalty_reasons", ""))
    if tier == "tier_1_reviewable_target_hypothesis" and "multiple_safety_liabilities" not in penalties:
        return "shortlist_review_candidate"
    if tier == "tier_1_reviewable_target_hypothesis":
        return "shortlist_candidate_with_safety_caution"
    if "not_direct_target" in tier:
        return "main_text_marker_or_module_context"
    if tier == "tier_2_supporting_target_hypothesis":
        return "supplement_or_secondary_candidate"
    return "background_or_sensitivity_candidate"


def write_summary(df: pd.DataFrame) -> None:
    tier_counts = df["priority_tier"].value_counts().to_dict()
    role_counts = df["manuscript_role"].value_counts().head(15).to_dict()
    shortlist = df[df["suggested_shortlist_label"].str.startswith("shortlist")].head(20)
    context = df[df["suggested_shortlist_label"] == "main_text_marker_or_module_context"].head(12)

    lines = [
        "# Preliminary target priority matrix summary",
        "",
        f"- Candidate rows: {len(df)}",
        f"- Priority tier counts: {tier_counts}",
        f"- Top manuscript-role counts: {role_counts}",
        "",
        "## Reviewable target hypotheses",
        "",
        "| rank | gene | score | druggability | role | actionability | cautions |",
        "| ---: | --- | ---: | ---: | --- | --- | --- |",
    ]
    for _, row in shortlist.iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    str(int(row["pre_rank_stability_rank"])),
                    str(row["gene_symbol"] or row["gene_id"]),
                    f"{float(row['final_target_priority_score_100_pre_rank_stability']):.3f}",
                    f"{float(row['druggability_score_20']):.3f}",
                    str(row["manuscript_role"]),
                    str(row["actionability_class"]),
                    str(row["review_flags"]),
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Main-text marker/context genes not to overclaim as direct targets",
            "",
            "| rank | gene | score | reason | role |",
            "| ---: | --- | ---: | --- | --- |",
        ]
    )
    for _, row in context.iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    str(int(row["pre_rank_stability_rank"])),
                    str(row["gene_symbol"] or row["gene_id"]),
                    f"{float(row['final_target_priority_score_100_pre_rank_stability']):.3f}",
                    str(row["review_flags"]),
                    str(row["manuscript_role"]),
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- This matrix is preliminary because rank-stability, directionality and manual target-safety review are still pending.",
            "- Genes with high biological evidence but minimal actionability are preserved as mechanistic context rather than removed.",
            "- Safety-cautioned targets remain eligible for discussion if local delivery, short-term intervention or disease-stage specificity can be justified.",
        ]
    )
    OUT_SUMMARY.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_review(df: pd.DataFrame) -> None:
    failures: list[str] = []
    if len(df) != 102:
        failures.append(f"Expected 102 candidates, observed {len(df)}.")
    if float(df["final_target_priority_score_100_pre_rank_stability"].max()) > 100:
        failures.append("A preliminary total score exceeds 100.")
    if df["pre_rank_stability_rank"].duplicated().any():
        failures.append("Ranks are not unique.")
    hspg2 = df[df["gene_symbol"] == "HSPG2"]
    if not hspg2.empty and not str(hspg2.iloc[0]["priority_tier"]).endswith("not_direct_target"):
        failures.append("HSPG2 is not marked as a non-direct-target ECM/context gene.")
    missing_symbol_high = df[
        (df["gene_symbol"].astype(str).str.len() == 0)
        & (df["priority_tier"] == "tier_1_reviewable_target_hypothesis")
    ]
    if not missing_symbol_high.empty:
        failures.append("One or more missing-symbol genes are tier-1 target hypotheses.")
    kdr = df[df["gene_symbol"] == "KDR"]
    if not kdr.empty and "systemic_modulation_safety_caution" not in str(kdr.iloc[0]["review_flags"]):
        failures.append("KDR lacks systemic safety caution.")

    status = "PASS" if not failures else "FAIL"
    lines = [
        "# Phase 12 self-review: preliminary target priority matrix",
        "",
        f"Status: {status}",
        "",
        "## Checks",
        "",
        f"- Candidate rows: {len(df)}",
        f"- Maximum score: {float(df['final_target_priority_score_100_pre_rank_stability'].max()):.3f}",
        f"- Tier-1 reviewable target hypotheses: {int((df['priority_tier'] == 'tier_1_reviewable_target_hypothesis').sum())}",
        f"- Tier-1 biologic marker/context genes: {int((df['priority_tier'] == 'tier_1_biologic_marker_or_context_gene_not_direct_target').sum())}",
        f"- Shortlist review candidates: {int(df['suggested_shortlist_label'].str.startswith('shortlist').sum())}",
        "",
        "## Guardrails",
        "",
        "- Preliminary target rank is not a final therapeutic recommendation.",
        "- Direct-target tiering requires druggability support; biology-only markers are explicitly separated.",
        "- This matrix should feed rank-stability and manual directionality review before a final manuscript shortlist is declared.",
        "",
    ]
    if failures:
        lines.extend(["## Failures", ""])
        lines.extend(f"- {failure}" for failure in failures)
    else:
        lines.extend(["## Decision", "", "The matrix passes automated scientific guardrail checks for use in downstream rank-stability analysis."])
    OUT_REVIEW.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    INTEGRATION.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT, sep="\t", keep_default_na=False)
    df = df.sort_values(
        ["final_target_priority_score_100_pre_rank_stability", "druggability_score_20"],
        ascending=False,
    ).reset_index(drop=True)
    df["pre_rank_stability_rank"] = range(1, len(df) + 1)
    df["priority_tier"] = df.apply(priority_tier, axis=1)
    df["manuscript_role"] = df.apply(manuscript_role, axis=1)
    df["review_flags"] = df.apply(review_flags, axis=1)
    df["suggested_shortlist_label"] = df.apply(suggested_shortlist_label, axis=1)

    front_cols = [
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
        "actionability_guardrail",
    ]
    ordered_cols = front_cols + [col for col in df.columns if col not in front_cols]
    df[ordered_cols].to_csv(OUT_MATRIX, sep="\t", index=False)
    write_summary(df)
    write_review(df)
    print(f"Wrote {OUT_MATRIX}")
    print(f"Wrote {OUT_SUMMARY}")
    print(f"Wrote {OUT_REVIEW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
