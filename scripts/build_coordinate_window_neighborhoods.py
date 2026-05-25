#!/usr/bin/env python3
"""Merge overlapping coordinate-window loci into cross-phenotype neighborhoods."""

from __future__ import annotations

import csv
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GWAS_DIR = PROJECT_ROOT / "results" / "gwas"
LOCI = GWAS_DIR / "coordinate_window_loci.tsv"
OUT_TSV = GWAS_DIR / "coordinate_window_neighborhoods.tsv"
OUT_MD = GWAS_DIR / "coordinate_window_neighborhoods_summary.md"


def read_loci() -> list[dict[str, str]]:
    with LOCI.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def chrom_sort_value(chrom: str) -> int:
    order = {str(i): i for i in range(1, 23)} | {"X": 23, "Y": 24, "MT": 25}
    return order.get(chrom, 99)


def merge_for_threshold(loci: list[dict[str, str]], threshold: str) -> list[dict[str, str]]:
    selected = [row for row in loci if row["threshold"] == threshold]
    selected.sort(key=lambda row: (chrom_sort_value(row["chrom"]), int(row["window_start"]), int(row["window_end"])))
    neighborhoods: list[dict[str, object]] = []

    for locus in selected:
        chrom = locus["chrom"]
        start = int(locus["window_start"])
        end = int(locus["window_end"])
        if not neighborhoods or neighborhoods[-1]["chrom"] != chrom or start > int(neighborhoods[-1]["end"]):
            neighborhoods.append({"chrom": chrom, "start": start, "end": end, "loci": [locus]})
        else:
            neighborhoods[-1]["end"] = max(int(neighborhoods[-1]["end"]), end)
            neighborhoods[-1]["loci"].append(locus)  # type: ignore[index]

    rows: list[dict[str, str]] = []
    for index, neighborhood in enumerate(neighborhoods, start=1):
        neighborhood_loci: list[dict[str, str]] = neighborhood["loci"]  # type: ignore[assignment]
        phenotypes = sorted({locus["phenotype"] for locus in neighborhood_loci})
        best_by_phenotype = []
        for phenotype in phenotypes:
            phenotype_loci = [locus for locus in neighborhood_loci if locus["phenotype"] == phenotype]
            best = min(phenotype_loci, key=lambda locus: float(locus["lead_p"]))
            best_by_phenotype.append(
                f"{phenotype}:{best['lead_snp']}:{best['lead_p']}:pos{best['lead_pos']}"
            )
        best_overall = min(neighborhood_loci, key=lambda locus: float(locus["lead_p"]))
        rows.append(
            {
                "neighborhood_id": f"{threshold}.n{index:03d}",
                "threshold": threshold,
                "chrom": str(neighborhood["chrom"]),
                "start": str(neighborhood["start"]),
                "end": str(neighborhood["end"]),
                "span_bp": str(int(neighborhood["end"]) - int(neighborhood["start"]) + 1),
                "phenotypes_present": ",".join(phenotypes),
                "n_phenotypes": str(len(phenotypes)),
                "n_loci": str(len(neighborhood_loci)),
                "best_lead_snp": best_overall["lead_snp"],
                "best_lead_p": best_overall["lead_p"],
                "best_phenotype": best_overall["phenotype"],
                "best_by_phenotype": ";".join(best_by_phenotype),
                "locus_ids": ",".join(locus["locus_id"] for locus in neighborhood_loci),
            }
        )
    return rows


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
        "locus_ids",
    ]
    with OUT_TSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# Coordinate-Window Genomic Neighborhoods",
        "",
        "## Scope",
        "",
        "Overlapping lead-SNP +/- 1 Mb coordinate-window loci were merged into connected genomic neighborhoods. This reduces pairwise duplicate counting but still does not replace LD-aware locus definition.",
        "",
        "## Counts",
        "",
        "| Threshold | Neighborhoods | Shared by >=2 phenotypes | Supported by all 3 phenotypes |",
        "|---|---:|---:|---:|",
    ]
    for threshold in sorted({row["threshold"] for row in rows}):
        subset = [row for row in rows if row["threshold"] == threshold]
        shared = sum(int(row["n_phenotypes"]) >= 2 for row in subset)
        all_three = sum(int(row["n_phenotypes"]) == 3 for row in subset)
        lines.append(f"| {threshold} | {len(subset):,} | {shared:,} | {all_three:,} |")

    lines.extend(
        [
            "",
            "## Genome-Wide Shared Neighborhoods",
            "",
            "| Neighborhood | Chr:span | Phenotypes | Best signal | Per-phenotype best leads |",
            "|---|---|---|---|---|",
        ]
    )
    genome_shared = sorted(
        [row for row in rows if row["threshold"] == "genome_wide" and int(row["n_phenotypes"]) >= 2],
        key=lambda row: float(row["best_lead_p"]),
    )
    for row in genome_shared:
        best = f"{row['best_phenotype']}:{row['best_lead_snp']}:{row['best_lead_p']}"
        lines.append(
            f"| {row['neighborhood_id']} | {row['chrom']}:{row['start']}-{row['end']} | {row['phenotypes_present']} | {best} | {row['best_by_phenotype']} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation Gate",
            "",
            "These neighborhoods are appropriate for prioritizing regions for LD clumping, gene mapping and cross-disease comparison. They should be reported as coordinate-window neighborhoods unless LD evidence confirms independent/shared signals.",
            "",
            "## Output",
            "",
            f"- Neighborhood table: `{OUT_TSV.relative_to(PROJECT_ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    loci = read_loci()
    rows: list[dict[str, str]] = []
    for threshold in sorted({row["threshold"] for row in loci}):
        rows.extend(merge_for_threshold(loci, threshold))
    rows.sort(key=lambda row: (row["threshold"], chrom_sort_value(row["chrom"]), int(row["start"])))
    write_outputs(rows)
    print(f"Wrote {OUT_TSV}")
    print(f"Wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
