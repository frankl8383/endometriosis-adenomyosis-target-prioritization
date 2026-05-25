#!/usr/bin/env python3
"""Build a GWAS-derived candidate gene universe for expression validation.

This is intentionally a candidate universe, not a target shortlist. It combines
shared coordinate-window neighborhoods, Ensembl pairwise LD sensitivity classes
and positional protein-coding gene overlap.
"""

from __future__ import annotations

import csv
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GWAS_DIR = PROJECT_ROOT / "results" / "gwas"
LD_NEIGHBORHOODS = GWAS_DIR / "ld_supported_shared_neighborhoods.tsv"
GENES = GWAS_DIR / "coordinate_window_neighborhood_genes.tsv"
OUT_TSV = GWAS_DIR / "gwas_candidate_gene_universe.tsv"
OUT_MD = GWAS_DIR / "gwas_candidate_gene_universe_summary.md"

MODULE_HINTS = {
    "WNT4": "hormone_development",
    "ESR1": "hormone_response",
    "CCDC170": "hormone_response",
    "FSHB": "endocrine_axis",
    "KDR": "angiogenesis",
    "KIT": "stromal_immune_development",
    "LY96": "immune_inflammation",
    "HRH1": "immune_inflammation",
    "HSPG2": "ecm_basement_membrane",
    "ITPR2": "calcium_signaling",
    "SSPN": "membrane_ecm",
    "JPH1": "calcium_signaling",
    "CDC42": "cell_migration_remodeling",
    "ATG7": "autophagy_stress",
}


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def priority_from_class(threshold: str, ld_class: str) -> str:
    if threshold == "genome_wide" and ld_class in {"same_lead_snp_supported", "high_ld_supported"}:
        return "tier1_genome_wide_shared_ld_supported"
    if threshold == "genome_wide" and ld_class == "moderate_ld_supported":
        return "tier2_genome_wide_shared_moderate_ld"
    if threshold == "genome_wide":
        return "tier3_genome_wide_shared_weak_or_unresolved_ld"
    if threshold == "suggestive" and ld_class in {"same_lead_snp_supported", "high_ld_supported"}:
        return "sensitivity_suggestive_shared_ld_supported"
    return "sensitivity_suggestive_or_unresolved"


def build_rows() -> list[dict[str, str]]:
    neighborhoods = read_table(LD_NEIGHBORHOODS)
    genes = read_table(GENES)
    neighborhood_by_id = {row["neighborhood_id"]: row for row in neighborhoods}
    rows: list[dict[str, str]] = []

    for gene in genes:
        neighborhood = neighborhood_by_id.get(gene["neighborhood_id"])
        if not neighborhood:
            continue
        priority = priority_from_class(neighborhood["threshold"], neighborhood["ld_neighborhood_class"])
        rows.append(
            {
                "gene_symbol": gene["gene_symbol"],
                "gene_id": gene["gene_id"],
                "gene_biotype": gene["gene_biotype"],
                "neighborhood_id": gene["neighborhood_id"],
                "threshold": neighborhood["threshold"],
                "genetic_priority": priority,
                "ld_neighborhood_class": neighborhood["ld_neighborhood_class"],
                "max_r2": neighborhood["max_r2"],
                "chrom": neighborhood["chrom"],
                "neighborhood_start": neighborhood["start"],
                "neighborhood_end": neighborhood["end"],
                "phenotypes_present": neighborhood["phenotypes_present"],
                "n_phenotypes": neighborhood["n_phenotypes"],
                "best_lead_snp": neighborhood["best_lead_snp"],
                "best_lead_p": neighborhood["best_lead_p"],
                "best_phenotype": neighborhood["best_phenotype"],
                "gene_start": gene["gene_start"],
                "gene_end": gene["gene_end"],
                "distance_to_best_lead_bp": gene["distance_to_best_lead_bp"],
                "module_hint_preliminary": MODULE_HINTS.get(gene["gene_symbol"], ""),
                "gene_description": gene["gene_description"],
            }
        )

    rows.sort(
        key=lambda row: (
            {
                "tier1_genome_wide_shared_ld_supported": 1,
                "tier2_genome_wide_shared_moderate_ld": 2,
                "tier3_genome_wide_shared_weak_or_unresolved_ld": 3,
                "sensitivity_suggestive_shared_ld_supported": 4,
                "sensitivity_suggestive_or_unresolved": 5,
            }.get(row["genetic_priority"], 99),
            float(row["best_lead_p"]),
            int(row["distance_to_best_lead_bp"]) if row["distance_to_best_lead_bp"].isdigit() else 10**12,
            row["gene_symbol"],
        )
    )
    return rows


def write_outputs(rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "gene_symbol",
        "gene_id",
        "gene_biotype",
        "neighborhood_id",
        "threshold",
        "genetic_priority",
        "ld_neighborhood_class",
        "max_r2",
        "chrom",
        "neighborhood_start",
        "neighborhood_end",
        "phenotypes_present",
        "n_phenotypes",
        "best_lead_snp",
        "best_lead_p",
        "best_phenotype",
        "gene_start",
        "gene_end",
        "distance_to_best_lead_bp",
        "module_hint_preliminary",
        "gene_description",
    ]
    with OUT_TSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    unique_genes = {row["gene_symbol"] or row["gene_id"] for row in rows}
    lines = [
        "# GWAS Candidate Gene Universe",
        "",
        "## Scope",
        "",
        "This table is the genetics-derived candidate gene universe for downstream expression, single-cell/spatial and druggability evidence integration. It is not a target shortlist.",
        "",
        "## Counts",
        "",
        f"- Neighborhood-gene records: {len(rows):,}",
        f"- Unique gene symbols/IDs: {len(unique_genes):,}",
        "",
        "| Genetic priority | Records | Unique genes |",
        "|---|---:|---:|",
    ]
    for priority in sorted({row["genetic_priority"] for row in rows}):
        subset = [row for row in rows if row["genetic_priority"] == priority]
        subset_genes = {row["gene_symbol"] or row["gene_id"] for row in subset}
        lines.append(f"| {priority} | {len(subset):,} | {len(subset_genes):,} |")

    lines.extend(
        [
            "",
            "## Tier 1/2 Near-Lead Candidates",
            "",
            "| Priority | Gene | Neighborhood | LD class | Best lead | Distance bp | Module hint |",
            "|---|---|---|---|---|---:|---|",
        ]
    )
    for row in [
        row
        for row in rows
        if row["genetic_priority"]
        in {"tier1_genome_wide_shared_ld_supported", "tier2_genome_wide_shared_moderate_ld"}
        and (row["distance_to_best_lead_bp"].isdigit() and int(row["distance_to_best_lead_bp"]) <= 750_000)
    ][:60]:
        best = f"{row['best_lead_snp']}:{row['best_lead_p']}"
        lines.append(
            f"| {row['genetic_priority']} | {row['gene_symbol'] or row['gene_id']} | {row['neighborhood_id']} | {row['ld_neighborhood_class']} | {best} | {row['distance_to_best_lead_bp']} | {row['module_hint_preliminary']} |"
        )

    lines.extend(
        [
            "",
            "## Interpretation Gate",
            "",
            "Genes in this table enter bulk-expression validation and cell-state localization. They cannot be called therapeutic targets until expression direction, cell/spatial localization, druggability and safety evidence are integrated.",
            "",
            "## Output",
            "",
            f"- Candidate universe: `{OUT_TSV.relative_to(PROJECT_ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    rows = build_rows()
    write_outputs(rows)
    print(f"Wrote {OUT_TSV}")
    print(f"Wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
