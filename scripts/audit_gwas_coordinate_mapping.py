#!/usr/bin/env python3
"""Independent audit of GWAS rsID coordinate mapping."""

from __future__ import annotations

import csv
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GWAS_DIR = PROJECT_ROOT / "results" / "gwas"
SNP_EVIDENCE = GWAS_DIR / "snps_for_coordinate_mapping.tsv"
COORDS = GWAS_DIR / "ensembl_grch38_variant_coordinates.tsv"
OUT_TSV = GWAS_DIR / "coordinate_mapping_audit.tsv"
OUT_MD = GWAS_DIR / "coordinate_mapping_audit.md"

PHENOTYPES = [
    "adenomyosis_EUR",
    "endometriosis_EUR_wo_23andMe",
    "endometriosis_wo_adenomyosis_EUR",
]

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_phenotypes(value: str) -> set[str]:
    return {item for item in value.split(",") if item}


def allele_tokens(allele_string: str) -> set[str]:
    return {token.upper() for token in allele_string.replace("|", "/").split("/") if token}


def complement(allele: str) -> str:
    return allele.upper().translate(COMPLEMENT)


def allele_status(row: dict[str, str], phenotype: str, coord: dict[str, str]) -> str:
    a1 = row.get(f"{phenotype}_Allele1", "").upper()
    a2 = row.get(f"{phenotype}_Allele2", "").upper()
    if not a1 or not a2:
        return "missing_gwas_alleles"
    if not coord.get("allele_string"):
        return "missing_ensembl_alleles"
    tokens = allele_tokens(coord["allele_string"])
    if a1 in tokens and a2 in tokens:
        return "exact_match"
    if complement(a1) in tokens and complement(a2) in tokens:
        return "strand_flip_match"
    return "incompatible"


def min_p(row: dict[str, str]) -> float:
    try:
        return float(row["min_p"])
    except ValueError:
        return 1.0


def main() -> int:
    evidence_rows = read_table(SNP_EVIDENCE)
    coord_rows = read_table(COORDS)
    coords_by_snp = {row["SNP"]: row for row in coord_rows}

    audit_rows: list[dict[str, str]] = []
    status_counts: dict[str, int] = {}
    genome_wide_total = 0
    genome_wide_mapped = 0
    suggestive_total = 0
    suggestive_mapped = 0
    allele_counts: dict[str, int] = {}

    for row in evidence_rows:
        snp = row["SNP"]
        coord = coords_by_snp.get(snp)
        mapping_status = coord["mapping_status"] if coord else "missing_from_coordinate_table"
        status_counts[mapping_status] = status_counts.get(mapping_status, 0) + 1
        is_mapped = mapping_status.startswith("mapped_")

        genome_wide_phenotypes = parse_phenotypes(row["genome_wide_phenotypes"])
        suggestive_phenotypes = parse_phenotypes(row["suggestive_phenotypes"])
        if genome_wide_phenotypes:
            genome_wide_total += 1
            if is_mapped:
                genome_wide_mapped += 1
        if suggestive_phenotypes:
            suggestive_total += 1
            if is_mapped:
                suggestive_mapped += 1

        per_pheno_statuses: list[str] = []
        if coord and is_mapped:
            for phenotype in PHENOTYPES:
                if phenotype not in parse_phenotypes(row["phenotypes_present"]):
                    continue
                status = allele_status(row, phenotype, coord)
                allele_counts[status] = allele_counts.get(status, 0) + 1
                per_pheno_statuses.append(f"{phenotype}:{status}")

        audit_rows.append(
            {
                "SNP": snp,
                "min_p": row["min_p"],
                "best_phenotype": row["best_phenotype"],
                "mapping_status": mapping_status,
                "chrom": "" if not coord else coord["chrom"],
                "start": "" if not coord else coord["start"],
                "end": "" if not coord else coord["end"],
                "allele_string": "" if not coord else coord["allele_string"],
                "genome_wide_phenotypes": row["genome_wide_phenotypes"],
                "suggestive_phenotypes": row["suggestive_phenotypes"],
                "allele_audit": ";".join(per_pheno_statuses),
            }
        )

    with OUT_TSV.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "SNP",
            "min_p",
            "best_phenotype",
            "mapping_status",
            "chrom",
            "start",
            "end",
            "allele_string",
            "genome_wide_phenotypes",
            "suggestive_phenotypes",
            "allele_audit",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(audit_rows, key=min_p))

    unmapped = [row for row in audit_rows if not row["mapping_status"].startswith("mapped_")]
    multiple = [row for row in audit_rows if row["mapping_status"] == "mapped_multiple_primary_grch38_chromosome"]
    top_unmapped = sorted(unmapped, key=min_p)[:10]
    top_multiple = sorted(multiple, key=min_p)[:10]

    lines = [
        "# Coordinate Mapping Audit",
        "",
        "## Pass/Fail Summary",
        "",
        f"- SNP evidence rows: {len(evidence_rows):,}",
        f"- Coordinate table rows: {len(coord_rows):,}",
        f"- Genome-wide significant rsIDs mapped: {genome_wide_mapped:,}/{genome_wide_total:,}",
        f"- Suggestive rsIDs mapped: {suggestive_mapped:,}/{suggestive_total:,}",
        "",
        "## Mapping Status Counts",
        "",
    ]
    for status, count in sorted(status_counts.items()):
        lines.append(f"- {status}: {count:,}")
    lines.extend(["", "## Allele Compatibility Counts", ""])
    for status, count in sorted(allele_counts.items()):
        lines.append(f"- {status}: {count:,}")
    lines.extend(["", "## Top Unmapped / Not Returned SNPs", ""])
    if top_unmapped:
        lines.extend(["| SNP | Min P | Best phenotype | Status |", "|---|---:|---|---|"])
        for row in top_unmapped:
            lines.append(f"| {row['SNP']} | {row['min_p']} | {row['best_phenotype']} | {row['mapping_status']} |")
    else:
        lines.append("- None.")
    lines.extend(["", "## Top Multiple-Mapping SNPs", ""])
    if top_multiple:
        lines.extend(["| SNP | Min P | Best phenotype | Primary locations |", "|---|---:|---|---|"])
        for row in top_multiple:
            lines.append(f"| {row['SNP']} | {row['min_p']} | {row['best_phenotype']} | {row['chrom']}:{row['start']}-{row['end']} |")
    else:
        lines.append("- None.")
    lines.extend(
        [
            "",
            "## Verdict",
            "",
            "PASS for coordinate-enabled GWAS evidence construction if all genome-wide significant rsIDs are mapped and allele incompatibilities are reviewed before directional cross-phenotype interpretation.",
            "",
            "## Output",
            "",
            f"- Audit table: `{OUT_TSV.relative_to(PROJECT_ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {OUT_TSV}")
    print(f"Wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
