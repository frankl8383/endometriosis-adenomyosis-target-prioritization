#!/usr/bin/env python3
"""Summarize LD support for shared coordinate-window neighborhoods."""

from __future__ import annotations

import csv
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GWAS_DIR = PROJECT_ROOT / "results" / "gwas"
NEIGHBORHOODS = GWAS_DIR / "coordinate_window_neighborhoods.tsv"
LD_TABLE = GWAS_DIR / "ensembl_pairwise_ld_sensitivity.tsv"
GENES = GWAS_DIR / "coordinate_window_neighborhood_genes.tsv"
OUT_TSV = GWAS_DIR / "ld_supported_shared_neighborhoods.tsv"
OUT_MD = GWAS_DIR / "ld_supported_shared_neighborhoods_summary.md"


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_locus_ids(value: str) -> set[str]:
    return {item for item in value.split(",") if item}


def gene_candidates(gene_rows: list[dict[str, str]], neighborhood_id: str, n: int = 8) -> str:
    subset = [row for row in gene_rows if row["neighborhood_id"] == neighborhood_id and row["gene_symbol"]]
    subset.sort(key=lambda row: (int(row["distance_to_best_lead_bp"]), row["gene_symbol"]))
    return ";".join(f"{row['gene_symbol']}({row['distance_to_best_lead_bp']}bp)" for row in subset[:n])


def ld_rows_for_neighborhood(ld_rows: list[dict[str, str]], locus_ids: set[str], threshold: str) -> list[dict[str, str]]:
    return [
        row
        for row in ld_rows
        if row["threshold"] == threshold and row["locus_a"] in locus_ids and row["locus_b"] in locus_ids
    ]


def max_float(values: list[str]) -> str:
    parsed = [float(value) for value in values if value]
    if not parsed:
        return ""
    return f"{max(parsed):.6f}"


def classify_neighborhood(rows: list[dict[str, str]]) -> str:
    classes = {row["ld_class"] for row in rows}
    if "same_snp" in classes:
        return "same_lead_snp_supported"
    if "high_ld_r2_ge_0.8" in classes:
        return "high_ld_supported"
    if "moderate_ld_r2_0.2_to_0.8" in classes:
        return "moderate_ld_supported"
    if "low_ld_r2_lt_0.2" in classes and len(classes - {"low_ld_r2_lt_0.2", "not_available"}) == 0:
        return "weak_or_distinct_ld"
    if classes == {"not_available"}:
        return "ld_not_available"
    return "mixed_unresolved"


def build_summary_rows() -> list[dict[str, str]]:
    neighborhoods = read_table(NEIGHBORHOODS)
    ld_rows = read_table(LD_TABLE)
    gene_rows = read_table(GENES)
    out_rows: list[dict[str, str]] = []

    for neighborhood in neighborhoods:
        if int(neighborhood["n_phenotypes"]) < 2:
            continue
        locus_ids = parse_locus_ids(neighborhood["locus_ids"])
        rows = ld_rows_for_neighborhood(ld_rows, locus_ids, neighborhood["threshold"])
        counts: dict[str, int] = {}
        for row in rows:
            counts[row["ld_class"]] = counts.get(row["ld_class"], 0) + 1
        out_rows.append(
            {
                "neighborhood_id": neighborhood["neighborhood_id"],
                "threshold": neighborhood["threshold"],
                "chrom": neighborhood["chrom"],
                "start": neighborhood["start"],
                "end": neighborhood["end"],
                "span_bp": neighborhood["span_bp"],
                "phenotypes_present": neighborhood["phenotypes_present"],
                "n_phenotypes": neighborhood["n_phenotypes"],
                "n_loci": neighborhood["n_loci"],
                "best_lead_snp": neighborhood["best_lead_snp"],
                "best_lead_p": neighborhood["best_lead_p"],
                "best_phenotype": neighborhood["best_phenotype"],
                "best_by_phenotype": neighborhood["best_by_phenotype"],
                "n_ld_pairs": str(len(rows)),
                "n_same_snp": str(counts.get("same_snp", 0)),
                "n_high_ld": str(counts.get("high_ld_r2_ge_0.8", 0)),
                "n_moderate_ld": str(counts.get("moderate_ld_r2_0.2_to_0.8", 0)),
                "n_low_ld": str(counts.get("low_ld_r2_lt_0.2", 0)),
                "n_ld_not_available": str(counts.get("not_available", 0)),
                "max_r2": max_float([row["r2"] for row in rows]),
                "ld_neighborhood_class": classify_neighborhood(rows) if rows else "no_pairwise_ld_test",
                "nearest_positional_genes": gene_candidates(gene_rows, neighborhood["neighborhood_id"]),
            }
        )

    out_rows.sort(key=lambda row: (row["threshold"], -int(row["n_phenotypes"]), float(row["best_lead_p"])))
    return out_rows


def write_outputs(rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "neighborhood_id",
        "threshold",
        "chrom",
        "start",
        "end",
        "span_bp",
        "phenotypes_present",
        "n_phenotypes",
        "n_loci",
        "best_lead_snp",
        "best_lead_p",
        "best_phenotype",
        "best_by_phenotype",
        "n_ld_pairs",
        "n_same_snp",
        "n_high_ld",
        "n_moderate_ld",
        "n_low_ld",
        "n_ld_not_available",
        "max_r2",
        "ld_neighborhood_class",
        "nearest_positional_genes",
    ]
    with OUT_TSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# LD-Supported Shared Coordinate-Window Neighborhoods",
        "",
        "## Scope",
        "",
        "This table summarizes pairwise Ensembl EUR LD among cross-phenotype lead SNPs within shared coordinate-window neighborhoods. It ranks evidence for shared genetic neighborhoods, while preserving the limitation that full PLINK-based LD clumping is still required.",
        "",
        "## Class Counts",
        "",
        "| Threshold | LD neighborhood class | Count |",
        "|---|---|---:|",
    ]
    for threshold in sorted({row["threshold"] for row in rows}):
        subset = [row for row in rows if row["threshold"] == threshold]
        for cls in sorted({row["ld_neighborhood_class"] for row in subset}):
            count = sum(row["ld_neighborhood_class"] == cls for row in subset)
            lines.append(f"| {threshold} | {cls} | {count:,} |")

    lines.extend(
        [
            "",
            "## Genome-Wide Shared Neighborhoods",
            "",
            "| Neighborhood | Chr:span | Phenotypes | Best lead | LD class | max r2 | Nearest positional genes |",
            "|---|---|---|---|---|---:|---|",
        ]
    )
    for row in [row for row in rows if row["threshold"] == "genome_wide"]:
        span = f"{row['chrom']}:{row['start']}-{row['end']}"
        best = f"{row['best_phenotype']}:{row['best_lead_snp']}:{row['best_lead_p']}"
        lines.append(
            f"| {row['neighborhood_id']} | {span} | {row['phenotypes_present']} | {best} | {row['ld_neighborhood_class']} | {row['max_r2'] or 'NA'} | {row['nearest_positional_genes']} |"
        )

    lines.extend(
        [
            "",
            "## Interpretation Gate",
            "",
            "Neighborhoods with same-lead, high-LD or moderate-LD support can be prioritized for gene mapping, expression localization and target scoring. Weak, mixed or unavailable LD neighborhoods remain candidates but should not be interpreted as shared genetic signals without PLINK/fine-mapping confirmation.",
            "",
            "## Output",
            "",
            f"- Summary table: `{OUT_TSV.relative_to(PROJECT_ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    rows = build_summary_rows()
    write_outputs(rows)
    print(f"Wrote {OUT_TSV}")
    print(f"Wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
