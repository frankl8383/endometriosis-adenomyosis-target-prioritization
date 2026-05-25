#!/usr/bin/env python3
"""Build manual directionality, safety and literature review for target shortlist."""

from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT_ROOT / "results"
INTEGRATION = RESULTS / "integration"
LIT = RESULTS / "literature"
PY_PKG_DIR = PROJECT_ROOT / "tools" / "py_packages"

if sys.version_info[:2] != (3, 12):
    raise SystemExit(
        "Use the bundled Python 3.12 runtime: "
        "python3 "
        "scripts/build_phase14_manual_target_review.py"
    )

sys.path.insert(0, str(PY_PKG_DIR))

import pandas as pd  # noqa: E402


PRELIM = INTEGRATION / "preliminary_target_priority_matrix.tsv"
STABILITY = INTEGRATION / "rank_stability_combined_matrix.tsv"
LITERATURE = LIT / "phase14_pubmed_literature_hits.tsv"

OUT_TSV = INTEGRATION / "phase14_manual_target_directionality_review.tsv"
OUT_SUMMARY = INTEGRATION / "phase14_manual_target_directionality_summary.md"
OUT_SELF_REVIEW = INTEGRATION / "phase14_manual_target_directionality_self_review.md"


CURATION = {
    "KDM1A": {
        "manual_category": "primary_shortlist_candidate",
        "proposed_direction": "inhibit_LSD1_KDM1A",
        "literature_level": "direct_functional_support_in_endometriosis_and_adenomyosis",
        "selected_pmids": ["39737119", "27062244", "24388204"],
        "manual_safety": "epigenetic_systemic_caution_but_drug_precedent_exists",
        "manual_rationale": "Direct studies report KDM1A/LSD1 inhibition suppressing adenomyosis stromal-cell phenotypes and endometriosis lesion growth; high rank stability and druggability support a primary epigenetic model-system hypothesis.",
        "claim_strength": "strong_translational_hypothesis_not_causal_proof",
    },
    "KDR": {
        "manual_category": "safety_limited_shortlist_candidate",
        "proposed_direction": "inhibit_or_locally_modulate_VEGFR2_KDR_axis",
        "literature_level": "genetic_and_angiogenesis_axis_support",
        "selected_pmids": ["27453397", "23635398", "26773192", "37210898", "39595579"],
        "manual_safety": "high_systemic_antiangiogenic_caution",
        "manual_rationale": "KDR has genetic and fibrovascular cell-context support, but systemic VEGFR2 inhibition has major vascular, wound-healing and reproductive safety concerns; manuscript claims should emphasize local/short-term modulation or use KDR as a fibrovascular signal rather than a simple repurposing target.",
        "claim_strength": "strong_biology_high_safety_caution",
    },
    "LY96": {
        "manual_category": "primary_shortlist_candidate",
        "proposed_direction": "inhibit_or_dampen_LY96_TLR4_MD2_inflammatory_signaling",
        "literature_level": "pathway_support_stronger_than_gene_specific_support",
        "selected_pmids": ["23855795", "19365133", "18596029", "23252918"],
        "manual_safety": "innate_immune_modulation_caution",
        "manual_rationale": "LY96 is best framed as the MD-2/TLR4 axis rather than a stand-alone gene claim. Rank stability and immune cell-context evidence are strong, while direct LY96-named endometriosis literature is limited; TLR4/endometriosis literature supports a cautious immune-axis hypothesis.",
        "claim_strength": "moderate_strong_axis_hypothesis",
    },
    "PDGFRA": {
        "manual_category": "primary_shortlist_candidate",
        "proposed_direction": "inhibit_or_rebalance_PDGFRA_stromal_fibrovascular_signaling",
        "literature_level": "endometrial_stromal_function_and_recent_endometriosis_multiomics_support",
        "selected_pmids": ["41721175", "16815388", "23349855", "22588000"],
        "manual_safety": "growth_factor_receptor_caution_moderate",
        "manual_rationale": "PDGFRA combines fibrovascular/stromal biology, high tractability and literature linking PDGF signaling to endometrial stromal proliferation/motility. It is a plausible stromal-remodeling target but not yet disease-specific causal proof.",
        "claim_strength": "moderate_strong_translational_hypothesis",
    },
    "ECE1": {
        "manual_category": "secondary_shortlist_candidate",
        "proposed_direction": "inhibit_or_modulate_endothelin_ECE1_axis",
        "literature_level": "endothelin_axis_support_limited_gene_specific_ECE1_support",
        "selected_pmids": ["19626996", "29034546", "41533156", "17712175", "33470475"],
        "manual_safety": "vascular_pain_axis_caution",
        "manual_rationale": "ECE1 is rank stable and druggable, with a plausible endothelin pain/vascular axis in endometriosis. Direct ECE1 disease evidence is thinner than KDM1A/KDR/KIT, so it should be secondary unless directionality strengthens.",
        "claim_strength": "moderate_axis_hypothesis",
    },
    "STK38L": {
        "manual_category": "exploratory_discovery_candidate",
        "proposed_direction": "uncertain",
        "literature_level": "little_or_no_direct_disease_literature",
        "selected_pmids": [],
        "manual_safety": "unknown_context_requires_validation",
        "manual_rationale": "STK38L is rank stable in the integrated matrix but lacks direct endometriosis/adenomyosis literature and clear perturbation direction; keep as a novel discovery candidate in supplement or discussion, not as a main target.",
        "claim_strength": "exploratory_only",
    },
    "SRD5A3": {
        "manual_category": "exploratory_discovery_candidate",
        "proposed_direction": "uncertain_steroid_glycosylation_context",
        "literature_level": "broad_steroid_endometrium_literature_no_direct_SRD5A3_support",
        "selected_pmids": [],
        "manual_safety": "unknown_context_requires_validation",
        "manual_rationale": "SRD5A3 is stable and druggable by generic criteria but has no direct disease literature in the targeted PubMed search; its role should be treated as exploratory.",
        "claim_strength": "exploratory_only",
    },
    "KIT": {
        "manual_category": "safety_limited_secondary_candidate",
        "proposed_direction": "inhibit_KIT_CSF1R_related_inflammatory_signaling",
        "literature_level": "direct_drug_perturbation_support_but_not_KIT_specific",
        "selected_pmids": ["38227801", "40472665", "40508002", "31424502"],
        "manual_safety": "hematologic_and_mast_cell_receptor_safety_caution",
        "manual_rationale": "Pexidartinib provides direct endometriosis perturbation evidence for CSF1R/KIT inhibition, but the evidence is not KIT-specific and KIT has notable haematological and mast-cell safety liabilities. Treat KIT as a safety-limited secondary immune-axis hypothesis rather than a primary KIT-specific target.",
        "claim_strength": "secondary_axis_hypothesis_safety_limited_not_KIT_specific",
    },
    "WNT4": {
        "manual_category": "genetic_mechanism_anchor_not_near_term_drug_target",
        "proposed_direction": "mechanism_context_direction_uncertain",
        "literature_level": "strong_genetic_and_reproductive_biology_support",
        "selected_pmids": ["31821471", "38354602", "26363035", "27453397"],
        "manual_safety": "developmental_pathway_druggability_caution",
        "manual_rationale": "WNT4 is a strong genetic and reproductive-development anchor, but drug direction and safety are unclear. It should support genetic-to-cell-state interpretation rather than be sold as a near-term repurposing target.",
        "claim_strength": "strong_mechanistic_anchor_not_direct_drug_target",
    },
    "ITPR2": {
        "manual_category": "exploratory_discovery_candidate",
        "proposed_direction": "uncertain_calcium_signaling",
        "literature_level": "little_or_no_direct_disease_literature",
        "selected_pmids": [],
        "manual_safety": "calcium_signaling_broad_systemic_caution",
        "manual_rationale": "ITPR2 is rank sensitive, has limited direct disease literature and no clear drug direction; keep as secondary/exploratory.",
        "claim_strength": "exploratory_only",
    },
    "C1QA": {
        "manual_category": "secondary_shortlist_candidate",
        "proposed_direction": "modulate_complement_C1q_axis_context_dependent",
        "literature_level": "complement_axis_support_with_recent_C1q_endometriosis_evidence",
        "selected_pmids": ["38983846", "29572748", "34099740", "36624343"],
        "manual_safety": "complement_systemic_immune_caution",
        "manual_rationale": "C1QA/C1q has complement and proangiogenic endometriosis evidence, but complement biology is context-dependent and systemic modulation is risky. Use as secondary immune-vascular target hypothesis.",
        "claim_strength": "moderate_axis_hypothesis",
    },
    "HSPG2": {
        "manual_category": "marker_context_not_direct_target",
        "proposed_direction": "not_direct_target",
        "literature_level": "uterine_ECM_marker_support_little_direct_target_support",
        "selected_pmids": ["15214943", "27619726", "42175707"],
        "manual_safety": "basement_membrane_structural_marker_not_actionable",
        "manual_rationale": "HSPG2/perlecan is a strong fibrovascular/ECM niche marker in this project, but lacks direct drug precedent and should not be overclaimed as a direct therapeutic target.",
        "claim_strength": "strong_context_marker",
    },
    "SSPN": {
        "manual_category": "marker_context_not_direct_target",
        "proposed_direction": "not_direct_target",
        "literature_level": "little_or_no_direct_disease_literature",
        "selected_pmids": [],
        "manual_safety": "limited_actionability",
        "manual_rationale": "SSPN is biologically strong in the matrix but has poor drug precedent and little disease literature; keep as fibrovascular/stromal context marker.",
        "claim_strength": "context_marker_only",
    },
    "ESR1": {
        "manual_category": "known_hormonal_axis_control_not_novel_shortlist",
        "proposed_direction": "context_dependent_hormonal_modulation",
        "literature_level": "large_known_hormone_axis_literature",
        "selected_pmids": [],
        "manual_safety": "high_systemic_hormonal_safety_and_cycle_confounding_caution",
        "manual_rationale": "ESR1 is biologically important but not a novel genetics-to-target discovery here; extensive safety, hormonal and cycle confounding issues make it better as a known-axis comparator.",
        "claim_strength": "known_axis_context",
    },
}


def selected_citations(pmids: list[str], lit: pd.DataFrame) -> tuple[str, str]:
    if not pmids:
        return "", ""
    lit = lit.copy()
    lit["pmid"] = lit["pmid"].astype(str)
    rows = []
    urls = []
    for pmid in pmids:
        hit = lit[lit["pmid"] == str(pmid)]
        if hit.empty:
            rows.append(f"{pmid}: not_fetched")
            urls.append(f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/")
            continue
        first = hit.iloc[0]
        rows.append(f"{pmid}: {first.get('year', '')}; {first.get('title', '')}")
        urls.append(str(first.get("pubmed_url", f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/")))
    return " || ".join(rows), "|".join(urls)


def direct_pubmed_count(gene: str, lit: pd.DataFrame) -> int:
    rows = lit[lit["query_key"] == f"gene_disease_{gene}"]
    if rows.empty:
        return 0
    return int(pd.to_numeric(rows["pubmed_total_count"], errors="coerce").fillna(0).max())


def evidence_tier(row: pd.Series) -> str:
    category = str(row["manual_category"])
    if category == "primary_shortlist_candidate":
        return "A_primary"
    if category == "safety_limited_shortlist_candidate":
        return "B_primary_safety_limited"
    if category in {"secondary_shortlist_candidate", "safety_limited_secondary_candidate"}:
        return "C_secondary"
    if "anchor" in category or "control" in category:
        return "D_context_or_anchor"
    if category == "marker_context_not_direct_target":
        return "E_marker_not_target"
    return "F_exploratory"


def main() -> int:
    prelim = pd.read_csv(PRELIM, sep="\t", keep_default_na=False)
    stability = pd.read_csv(STABILITY, sep="\t", keep_default_na=False)
    lit = pd.read_csv(LITERATURE, sep="\t", keep_default_na=False)

    stability_cols = [
        "gene_id",
        "rank_stability_class",
        "bootstrap_top10_frequency",
        "leave_one_layer_top20_count",
        "gene_label_permutation_p_ge_observed",
    ]
    merged = prelim.merge(stability[stability_cols], on="gene_id", how="left", suffixes=("", "_stability"))
    rows = []
    for gene, curation in CURATION.items():
        hit = merged[merged["gene_symbol"] == gene]
        if hit.empty:
            base = {"gene_symbol": gene, "gene_id": ""}
        else:
            base = hit.iloc[0].to_dict()
        pmids = list(curation["selected_pmids"])
        citation_text, citation_urls = selected_citations(pmids, lit)
        row = {
            "gene_symbol": gene,
            "gene_id": base.get("gene_id", ""),
            "manual_evidence_tier": "",
            "manual_category": curation["manual_category"],
            "proposed_direction": curation["proposed_direction"],
            "literature_level": curation["literature_level"],
            "selected_pmids": "|".join(pmids),
            "selected_citations": citation_text,
            "selected_pubmed_urls": citation_urls,
            "direct_gene_pubmed_total_count": direct_pubmed_count(gene, lit),
            "manual_safety": curation["manual_safety"],
            "manual_rationale": curation["manual_rationale"],
            "claim_strength": curation["claim_strength"],
            "rank_stability_class": base.get("rank_stability_class", ""),
            "bootstrap_top10_frequency": base.get("bootstrap_top10_frequency", ""),
            "leave_one_layer_top20_count": base.get("leave_one_layer_top20_count", ""),
            "gene_label_permutation_p_ge_observed": base.get("gene_label_permutation_p_ge_observed", ""),
            "pre_rank_stability_rank": base.get("pre_rank_stability_rank", ""),
            "final_target_priority_score_100_pre_rank_stability": base.get("final_target_priority_score_100_pre_rank_stability", ""),
            "dominant_cross_disease_axis": base.get("dominant_cross_disease_axis", ""),
            "bulk_support_class": base.get("bulk_support_class", ""),
            "cycle_driven_flag": "",
            "lesion_logFC": "",
            "lesion_adj_p": "",
            "endo_best_pseudobulk_effect": base.get("endo_best_pseudobulk_effect", ""),
            "endo_best_pseudobulk_p_value": base.get("endo_best_pseudobulk_p_value", ""),
            "adeno_best_sample_level_effect": base.get("adeno_best_sample_level_effect", ""),
            "adeno_best_sample_level_p_value": base.get("adeno_best_sample_level_p_value", ""),
            "actionability_class": base.get("actionability_class", ""),
            "drug_direction_evidence": base.get("drug_direction_evidence", ""),
            "druggability_penalty_reasons": base.get("druggability_penalty_reasons", ""),
            "review_flags": base.get("review_flags", ""),
        }
        rows.append(row)

    out = pd.DataFrame(rows)
    out["manual_evidence_tier"] = out.apply(evidence_tier, axis=1)
    out = out.sort_values(["manual_evidence_tier", "final_target_priority_score_100_pre_rank_stability"], ascending=[True, False])
    out.to_csv(OUT_TSV, sep="\t", index=False)

    primary = out[out["manual_evidence_tier"].isin(["A_primary", "B_primary_safety_limited"])]
    secondary = out[out["manual_evidence_tier"].isin(["C_secondary", "D_context_or_anchor"])]
    lines = [
        "# Phase 14 manual target directionality and safety review",
        "",
        f"- Reviewed genes: {len(out)}",
        f"- Primary or primary safety-limited shortlist candidates: {len(primary)}",
        f"- Secondary/context candidates: {len(secondary)}",
        "",
        "## Primary Manuscript Shortlist Input",
        "",
        "| tier | gene | direction | axis | stability | safety | selected PMIDs |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for _, row in primary.iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["manual_evidence_tier"]),
                    str(row["gene_symbol"]),
                    str(row["proposed_direction"]),
                    str(row["dominant_cross_disease_axis"]),
                    str(row["rank_stability_class"]),
                    str(row["manual_safety"]),
                    str(row["selected_pmids"]),
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "## Secondary or Context Use",
            "",
            "| tier | gene | manuscript use | reason |",
            "| --- | --- | --- | --- |",
        ]
    )
    for _, row in secondary.iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["manual_evidence_tier"]),
                    str(row["gene_symbol"]),
                    str(row["manual_category"]),
                    str(row["manual_rationale"]),
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "## Interpretation Rules",
            "",
            "- Primary shortlist genes still require careful wording as candidate therapeutic hypotheses, not validated therapies.",
            "- Safety-limited genes can appear in the main shortlist only with explicit systemic-toxicity and local-delivery caveats.",
            "- Stable computational rank is not sufficient for main-text target status when direct literature or directionality is weak.",
            "- Marker/context genes are retained for mechanism figures and module interpretation but not promoted as direct drug targets.",
        ]
    )
    OUT_SUMMARY.write_text("\n".join(lines) + "\n", encoding="utf-8")

    failures = []
    stable_genes = set(
        stability[
            stability["rank_stability_class"].isin(["rank_stable_target_hypothesis", "moderately_stable_target_hypothesis"])
        ]["gene_symbol"]
    )
    reviewed_genes = set(out["gene_symbol"])
    missing_stable = sorted(stable_genes - reviewed_genes)
    if missing_stable:
        failures.append(f"Stable genes missing manual review: {', '.join(missing_stable)}")
    forbidden_primary = {"HSPG2", "SSPN", "STK38L", "SRD5A3", "ITPR2"}
    bad_primary = sorted(
        set(out[out["manual_evidence_tier"].isin(["A_primary", "B_primary_safety_limited"])]["gene_symbol"]) & forbidden_primary
    )
    if bad_primary:
        failures.append(f"Genes with weak direct actionability were incorrectly primary: {', '.join(bad_primary)}")
    for gene in ["KDR", "KIT", "ESR1"]:
        row = out[out["gene_symbol"] == gene]
        if not row.empty and "caution" not in str(row.iloc[0]["manual_safety"]).lower():
            failures.append(f"{gene} lacks explicit safety caution.")
    kit_row = out[out["gene_symbol"] == "KIT"]
    if not kit_row.empty and str(kit_row.iloc[0]["manual_evidence_tier"]) != "C_secondary":
        failures.append("KIT should be downgraded to a safety-limited secondary hypothesis.")
    primary_no_pmids = primary[primary["selected_pmids"].astype(str).str.len() == 0]["gene_symbol"].tolist()
    if primary_no_pmids:
        failures.append(f"Primary/safety-limited candidates without selected PMIDs: {', '.join(primary_no_pmids)}")
    if int((out["manual_evidence_tier"] == "A_primary").sum()) < 3:
        failures.append("Fewer than three primary shortlist candidates after manual review.")

    status = "PASS" if not failures else "FAIL"
    review_lines = [
        "# Phase 14 self-review: manual target directionality and safety",
        "",
        f"Status: {status}",
        "",
        "## Checks",
        "",
        f"- Reviewed genes: {len(out)}",
        f"- Stable genes reviewed: {len(stable_genes - (stable_genes - reviewed_genes))}/{len(stable_genes)}",
        f"- A_primary candidates: {int((out['manual_evidence_tier'] == 'A_primary').sum())}",
        f"- B_primary_safety_limited candidates: {int((out['manual_evidence_tier'] == 'B_primary_safety_limited').sum())}",
        f"- C_secondary candidates: {int((out['manual_evidence_tier'] == 'C_secondary').sum())}",
        f"- Marker/context genes: {int((out['manual_evidence_tier'] == 'E_marker_not_target').sum())}",
        "",
        "## Guardrails",
        "",
        "- No therapeutic efficacy claim is made from public-data prioritisation alone.",
        "- Directionality is considered provisional unless drug mechanism, disease expression and functional literature agree.",
        "- Cycle-driven expression flags weaken expression-direction claims and should be discussed in limitations.",
        "- Adenomyosis h5ad evidence remains single-cell/cell-label evidence, not spatial evidence.",
        "",
    ]
    if failures:
        review_lines.extend(["## Failures", ""])
        review_lines.extend(f"- {failure}" for failure in failures)
    else:
        review_lines.extend(
            [
                "## Decision",
                "",
                "Manual directionality/safety review passes the predefined guardrails. Use this table as the input for final Figure 6/table shortlist drafting, while preserving explicit caveats.",
            ]
        )
    OUT_SELF_REVIEW.write_text("\n".join(review_lines) + "\n", encoding="utf-8")

    print(f"Wrote {OUT_TSV}")
    print(f"Wrote {OUT_SUMMARY}")
    print(f"Wrote {OUT_SELF_REVIEW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
