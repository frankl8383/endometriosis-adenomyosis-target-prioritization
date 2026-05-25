#!/usr/bin/env python3
"""Build coordinate-window GWAS loci before LD-aware clumping is available.

This is a conservative planning layer, not a substitute for LD clumping. For
each phenotype and threshold, variants are sorted by P value and greedily
assigned to the best unassigned lead SNP within +/- 1 Mb on the same chromosome.
"""

from __future__ import annotations

import csv
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GWAS_DIR = PROJECT_ROOT / "results" / "gwas"
SNP_EVIDENCE = GWAS_DIR / "snps_for_coordinate_mapping.tsv"
COORDS = GWAS_DIR / "ensembl_grch38_variant_coordinates.tsv"
OUT_TSV = GWAS_DIR / "coordinate_window_loci.tsv"
OUT_MEMBERS_TSV = GWAS_DIR / "coordinate_window_locus_members.tsv"
OUT_MD = GWAS_DIR / "coordinate_window_loci_summary.md"

PHENOTYPES = [
    "adenomyosis_EUR",
    "endometriosis_EUR_wo_23andMe",
    "endometriosis_wo_adenomyosis_EUR",
]

THRESHOLDS = [
    ("genome_wide", 5e-8),
    ("suggestive", 1e-5),
]

WINDOW_BP = 1_000_000


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def chrom_sort_value(chrom: str) -> int:
    order = {str(i): i for i in range(1, 23)} | {"X": 23, "Y": 24, "MT": 25}
    return order.get(chrom, 99)


def float_or_none(value: str) -> float | None:
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def make_variant_records(evidence: list[dict[str, str]], coords_by_snp: dict[str, dict[str, str]]) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for row in evidence:
        coord = coords_by_snp.get(row["SNP"])
        if coord is None or not coord["mapping_status"].startswith("mapped_"):
            continue
        if not coord["start"].isdigit():
            continue
        for phenotype in PHENOTYPES:
            p_value = float_or_none(row.get(f"{phenotype}_P", ""))
            if p_value is None:
                continue
            records.append(
                {
                    "phenotype": phenotype,
                    "SNP": row["SNP"],
                    "chrom": coord["chrom"],
                    "pos": int(coord["start"]),
                    "P": p_value,
                    "Z": row.get(f"{phenotype}_Z", ""),
                    "N": row.get(f"{phenotype}_N", ""),
                    "Allele1": row.get(f"{phenotype}_Allele1", ""),
                    "Allele2": row.get(f"{phenotype}_Allele2", ""),
                    "Freq1": row.get(f"{phenotype}_Freq1", ""),
                    "mapping_status": coord["mapping_status"],
                }
            )
    return records


def clump_variants(records: list[dict[str, object]], phenotype: str, threshold_name: str, threshold: float) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    eligible = [
        record
        for record in records
        if record["phenotype"] == phenotype and float(record["P"]) < threshold
    ]
    eligible.sort(key=lambda record: (float(record["P"]), chrom_sort_value(str(record["chrom"])), int(record["pos"])))

    assigned: set[str] = set()
    loci: list[dict[str, str]] = []
    member_rows: list[dict[str, str]] = []
    locus_index = 0

    for lead in eligible:
        lead_snp = str(lead["SNP"])
        if lead_snp in assigned:
            continue
        locus_index += 1
        chrom = str(lead["chrom"])
        lead_pos = int(lead["pos"])
        members = [
            record
            for record in eligible
            if str(record["SNP"]) not in assigned
            and str(record["chrom"]) == chrom
            and abs(int(record["pos"]) - lead_pos) <= WINDOW_BP
        ]
        for member in members:
            assigned.add(str(member["SNP"]))

        member_positions = [int(member["pos"]) for member in members]
        member_ps = [float(member["P"]) for member in members]
        locus_id = f"{phenotype}.{threshold_name}.cw{locus_index:03d}"
        n_multiple_mapping = sum(
            1 for member in members if str(member["mapping_status"]) == "mapped_multiple_primary_grch38_chromosome"
        )
        loci.append(
            {
                "locus_id": locus_id,
                "phenotype": phenotype,
                "threshold": threshold_name,
                "threshold_p": f"{threshold:.1e}",
                "chrom": chrom,
                "lead_snp": lead_snp,
                "lead_pos": str(lead_pos),
                "lead_p": f"{float(lead['P']):.3e}",
                "lead_z": str(lead["Z"]),
                "lead_allele1": str(lead["Allele1"]),
                "lead_allele2": str(lead["Allele2"]),
                "lead_freq1": str(lead["Freq1"]),
                "window_start": str(max(1, lead_pos - WINDOW_BP)),
                "window_end": str(lead_pos + WINDOW_BP),
                "member_span_start": str(min(member_positions)),
                "member_span_end": str(max(member_positions)),
                "n_snps": str(len(members)),
                "n_multiple_mapping_snps": str(n_multiple_mapping),
                "min_member_p": f"{min(member_ps):.3e}",
                "max_member_p": f"{max(member_ps):.3e}",
                "members": ",".join(str(member["SNP"]) for member in members[:50]),
                "members_truncated": "yes" if len(members) > 50 else "no",
            }
        )
        for member in members:
            member_rows.append(
                {
                    "locus_id": locus_id,
                    "phenotype": phenotype,
                    "threshold": threshold_name,
                    "SNP": str(member["SNP"]),
                    "chrom": str(member["chrom"]),
                    "pos": str(member["pos"]),
                    "P": f"{float(member['P']):.3e}",
                    "Z": str(member["Z"]),
                    "is_lead": "yes" if str(member["SNP"]) == lead_snp else "no",
                }
            )
    loci.sort(key=lambda row: (chrom_sort_value(row["chrom"]), int(row["lead_pos"]), row["phenotype"], row["threshold"]))
    return loci, member_rows


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(loci: list[dict[str, str]], member_rows: list[dict[str, str]]) -> None:
    lines = [
        "# Coordinate-Window GWAS Loci",
        "",
        "## Scope",
        "",
        f"Mapped SNPs were greedily grouped by phenotype and threshold using lead SNP +/- {WINDOW_BP:,} bp windows. This is a coordinate-window screen, not an LD-independent locus definition.",
        "",
        "## Locus Counts",
        "",
        "| Phenotype | Threshold | Loci | Member SNPs |",
        "|---|---|---:|---:|",
    ]
    for phenotype in PHENOTYPES:
        for threshold_name, _threshold in THRESHOLDS:
            n_loci = sum(1 for row in loci if row["phenotype"] == phenotype and row["threshold"] == threshold_name)
            n_members = sum(1 for row in member_rows if row["phenotype"] == phenotype and row["threshold"] == threshold_name)
            lines.append(f"| {phenotype} | {threshold_name} | {n_loci:,} | {n_members:,} |")

    lines.extend(
        [
            "",
            "## Top Genome-Wide Coordinate Windows",
            "",
            "| Phenotype | Locus | Lead SNP | Chr:pos | P | SNPs |",
            "|---|---|---|---|---:|---:|",
        ]
    )
    top_gw = sorted(
        [row for row in loci if row["threshold"] == "genome_wide"],
        key=lambda row: float(row["lead_p"]),
    )[:15]
    for row in top_gw:
        lines.append(
            f"| {row['phenotype']} | {row['locus_id']} | {row['lead_snp']} | {row['chrom']}:{row['lead_pos']} | {row['lead_p']} | {row['n_snps']} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation Gate",
            "",
            "Use this table for planning, figure prototyping and positional gene-window preparation. Manuscript-level independent loci still require LD-aware clumping with a matched reference panel, or an explicit sensitivity analysis showing that coordinate windows are robust.",
            "",
            "## Outputs",
            "",
            f"- Locus table: `{OUT_TSV.relative_to(PROJECT_ROOT)}`",
            f"- Locus-member table: `{OUT_MEMBERS_TSV.relative_to(PROJECT_ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    evidence = read_table(SNP_EVIDENCE)
    coords = read_table(COORDS)
    coords_by_snp = {row["SNP"]: row for row in coords}
    variant_records = make_variant_records(evidence, coords_by_snp)

    all_loci: list[dict[str, str]] = []
    all_members: list[dict[str, str]] = []
    for phenotype in PHENOTYPES:
        for threshold_name, threshold in THRESHOLDS:
            loci, members = clump_variants(variant_records, phenotype, threshold_name, threshold)
            all_loci.extend(loci)
            all_members.extend(members)

    locus_fields = [
        "locus_id",
        "phenotype",
        "threshold",
        "threshold_p",
        "chrom",
        "lead_snp",
        "lead_pos",
        "lead_p",
        "lead_z",
        "lead_allele1",
        "lead_allele2",
        "lead_freq1",
        "window_start",
        "window_end",
        "member_span_start",
        "member_span_end",
        "n_snps",
        "n_multiple_mapping_snps",
        "min_member_p",
        "max_member_p",
        "members",
        "members_truncated",
    ]
    member_fields = ["locus_id", "phenotype", "threshold", "SNP", "chrom", "pos", "P", "Z", "is_lead"]
    write_tsv(OUT_TSV, all_loci, locus_fields)
    write_tsv(OUT_MEMBERS_TSV, all_members, member_fields)
    write_summary(all_loci, all_members)
    print(f"Wrote {OUT_TSV}")
    print(f"Wrote {OUT_MEMBERS_TSV}")
    print(f"Wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
