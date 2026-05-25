#!/usr/bin/env python3
"""Compare coordinate-window GWAS loci across phenotypes."""

from __future__ import annotations

import csv
from itertools import combinations
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GWAS_DIR = PROJECT_ROOT / "results" / "gwas"
LOCI = GWAS_DIR / "coordinate_window_loci.tsv"
OUT_TSV = GWAS_DIR / "coordinate_window_locus_overlaps.tsv"
OUT_MD = GWAS_DIR / "coordinate_window_locus_overlaps_summary.md"

PHENOTYPES = [
    "adenomyosis_EUR",
    "endometriosis_EUR_wo_23andMe",
    "endometriosis_wo_adenomyosis_EUR",
]


def read_loci() -> list[dict[str, str]]:
    with LOCI.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def overlap_bp(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def compare_loci(loci: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for threshold in sorted({row["threshold"] for row in loci}):
        threshold_loci = [row for row in loci if row["threshold"] == threshold]
        by_phenotype = {
            phenotype: [row for row in threshold_loci if row["phenotype"] == phenotype]
            for phenotype in PHENOTYPES
        }
        for phenotype_a, phenotype_b in combinations(PHENOTYPES, 2):
            for locus_a in by_phenotype[phenotype_a]:
                for locus_b in by_phenotype[phenotype_b]:
                    if locus_a["chrom"] != locus_b["chrom"]:
                        continue
                    a_start = int(locus_a["window_start"])
                    a_end = int(locus_a["window_end"])
                    b_start = int(locus_b["window_start"])
                    b_end = int(locus_b["window_end"])
                    ov = overlap_bp(a_start, a_end, b_start, b_end)
                    lead_distance = abs(int(locus_a["lead_pos"]) - int(locus_b["lead_pos"]))
                    if ov == 0 and lead_distance > 2_000_000:
                        continue
                    a_size = a_end - a_start + 1
                    b_size = b_end - b_start + 1
                    rows.append(
                        {
                            "threshold": threshold,
                            "phenotype_a": phenotype_a,
                            "locus_a": locus_a["locus_id"],
                            "lead_snp_a": locus_a["lead_snp"],
                            "lead_pos_a": locus_a["lead_pos"],
                            "lead_p_a": locus_a["lead_p"],
                            "phenotype_b": phenotype_b,
                            "locus_b": locus_b["locus_id"],
                            "lead_snp_b": locus_b["lead_snp"],
                            "lead_pos_b": locus_b["lead_pos"],
                            "lead_p_b": locus_b["lead_p"],
                            "chrom": locus_a["chrom"],
                            "window_overlap_bp": str(ov),
                            "reciprocal_overlap_min": f"{min(ov / a_size, ov / b_size):.3f}",
                            "lead_distance_bp": str(lead_distance),
                            "relationship": "overlapping_windows" if ov > 0 else "nearby_leads_no_window_overlap",
                        }
                    )
    rows.sort(key=lambda row: (row["threshold"], row["chrom"], int(row["lead_pos_a"]), int(row["lead_pos_b"])))
    return rows


def write_outputs(rows: list[dict[str, str]], loci: list[dict[str, str]]) -> None:
    fieldnames = [
        "threshold",
        "phenotype_a",
        "locus_a",
        "lead_snp_a",
        "lead_pos_a",
        "lead_p_a",
        "phenotype_b",
        "locus_b",
        "lead_snp_b",
        "lead_pos_b",
        "lead_p_b",
        "chrom",
        "window_overlap_bp",
        "reciprocal_overlap_min",
        "lead_distance_bp",
        "relationship",
    ]
    with OUT_TSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# Coordinate-Window Locus Overlap",
        "",
        "## Scope",
        "",
        "This compares lead-SNP +/- 1 Mb coordinate windows across phenotypes. Overlap here is evidence for shared genomic neighborhoods, not proof of the same LD signal or causal gene.",
        "",
        "## Pairwise Overlaps",
        "",
        "| Threshold | Phenotype pair | Overlapping windows | Nearby non-overlapping lead pairs |",
        "|---|---|---:|---:|",
    ]
    for threshold in sorted({row["threshold"] for row in loci}):
        for phenotype_a, phenotype_b in combinations(PHENOTYPES, 2):
            pair_rows = [
                row
                for row in rows
                if row["threshold"] == threshold
                and row["phenotype_a"] == phenotype_a
                and row["phenotype_b"] == phenotype_b
            ]
            overlapping = sum(row["relationship"] == "overlapping_windows" for row in pair_rows)
            nearby = sum(row["relationship"] == "nearby_leads_no_window_overlap" for row in pair_rows)
            lines.append(f"| {threshold} | {phenotype_a} vs {phenotype_b} | {overlapping:,} | {nearby:,} |")

    lines.extend(
        [
            "",
            "## Strongest Genome-Wide Shared Neighborhoods",
            "",
            "| Pair | Chr | Lead A | P A | Lead B | P B | Lead distance bp |",
            "|---|---|---|---:|---|---:|---:|",
        ]
    )
    top = sorted(
        [
            row
            for row in rows
            if row["threshold"] == "genome_wide" and row["relationship"] == "overlapping_windows"
        ],
        key=lambda row: max(float(row["lead_p_a"]), float(row["lead_p_b"])),
    )[:15]
    for row in top:
        pair = f"{row['phenotype_a']} vs {row['phenotype_b']}"
        lines.append(
            f"| {pair} | {row['chrom']} | {row['lead_snp_a']} | {row['lead_p_a']} | {row['lead_snp_b']} | {row['lead_p_b']} | {row['lead_distance_bp']} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation Gate",
            "",
            "Use overlapping windows to prioritize cross-phenotype regions for LD-aware clumping and gene mapping. Do not interpret these overlaps as colocalization.",
            "",
            "## Output",
            "",
            f"- Overlap table: `{OUT_TSV.relative_to(PROJECT_ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    loci = read_loci()
    rows = compare_loci(loci)
    write_outputs(rows, loci)
    print(f"Wrote {OUT_TSV}")
    print(f"Wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
