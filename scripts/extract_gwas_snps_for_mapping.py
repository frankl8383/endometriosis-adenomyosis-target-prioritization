#!/usr/bin/env python3
"""Extract GWAS SNPs that need coordinate mapping.

The Koller et al. GWAS files contain rsIDs but no chromosome/base-pair columns.
This script builds a compact rsID universe for downstream coordinate lookup:

- every SNP with P < 1e-5 in any primary GWAS file
- every SNP in the per-phenotype top-5000 tables

It keeps per-phenotype alleles, Z, P, N and direction because these fields are
needed later for harmonization and shared/specific evidence classification.
"""

from __future__ import annotations

import csv
import gzip
import math
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW = PROJECT_ROOT / "data" / "raw_downloads"
OUT_DIR = PROJECT_ROOT / "results" / "gwas"
OUT_TSV = OUT_DIR / "snps_for_coordinate_mapping.tsv"
OUT_TXT = OUT_DIR / "snps_for_coordinate_mapping.txt"

GWAS_FILES = [
    (
        "adenomyosis_EUR",
        RAW / "Zenodo_18983492__adenomyosis_EUR.txt.gz",
    ),
    (
        "endometriosis_EUR_wo_23andMe",
        RAW / "Zenodo_18983492__endometriosis_EUR_wo_23andMe.txt.gz",
    ),
    (
        "endometriosis_wo_adenomyosis_EUR",
        RAW / "Zenodo_18983492__endometriosis_wo.adenomyosis_EUR.txt.gz",
    ),
]

COLUMN_ALIASES = {
    "SNP": ["SNP"],
    "Allele1": ["Allele1"],
    "Allele2": ["Allele2"],
    "Freq1": ["Freq1"],
    "N": ["N"],
    "Z": ["Z"],
    "P": ["P", "P.value"],
    "Direction": ["Direction", "DIRECTION"],
    "HetPVal": ["HetPVal"],
}

TOP_N = 5000
SUGGESTIVE = 1e-5
GENOME_WIDE = 5e-8
PALINDROMIC = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}


def resolve_columns(header: list[str]) -> dict[str, str]:
    header_set = set(header)
    mapping: dict[str, str] = {}
    missing: list[str] = []
    for logical_name, aliases in COLUMN_ALIASES.items():
        match = next((alias for alias in aliases if alias in header_set), None)
        if match is None:
            missing.append(logical_name)
        else:
            mapping[logical_name] = match
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")
    return mapping


def parse_float(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"Non-finite numeric value: {value}")
    return parsed


def is_palindromic(a1: str, a2: str) -> bool:
    return (a1.upper(), a2.upper()) in PALINDROMIC


def ensure_record(records: dict[str, dict[str, object]], snp: str) -> dict[str, object]:
    if snp not in records:
        records[snp] = {
            "SNP": snp,
            "include_reasons": set(),
            "phenotypes_present": set(),
            "genome_wide_phenotypes": set(),
            "suggestive_phenotypes": set(),
            "top5000_phenotypes": set(),
            "min_p": 1.0,
            "best_phenotype": "",
        }
    return records[snp]


def add_row(
    records: dict[str, dict[str, object]],
    phenotype: str,
    row: dict[str, str],
    reason: str,
) -> None:
    snp = row["SNP"].strip()
    if not snp:
        return

    record = ensure_record(records, snp)
    include_reasons: set[str] = record["include_reasons"]  # type: ignore[assignment]
    phenotypes_present: set[str] = record["phenotypes_present"]  # type: ignore[assignment]
    genome_wide: set[str] = record["genome_wide_phenotypes"]  # type: ignore[assignment]
    suggestive: set[str] = record["suggestive_phenotypes"]  # type: ignore[assignment]
    top5000: set[str] = record["top5000_phenotypes"]  # type: ignore[assignment]

    p_value = parse_float(row["P"])
    include_reasons.add(reason)
    phenotypes_present.add(phenotype)
    if p_value < GENOME_WIDE:
        genome_wide.add(phenotype)
    if p_value < SUGGESTIVE:
        suggestive.add(phenotype)
    if reason == "top5000":
        top5000.add(phenotype)

    if p_value < float(record["min_p"]):
        record["min_p"] = p_value
        record["best_phenotype"] = phenotype

    prefix = phenotype
    record[f"{prefix}_P"] = row["P"]
    record[f"{prefix}_Z"] = row["Z"]
    record[f"{prefix}_N"] = row["N"]
    record[f"{prefix}_Allele1"] = row["Allele1"]
    record[f"{prefix}_Allele2"] = row["Allele2"]
    record[f"{prefix}_Freq1"] = row["Freq1"]
    record[f"{prefix}_Direction"] = row["Direction"]
    record[f"{prefix}_HetPVal"] = row["HetPVal"]
    record[f"{prefix}_palindromic"] = "yes" if is_palindromic(row["Allele1"], row["Allele2"]) else "no"


def load_top5000(records: dict[str, dict[str, object]]) -> None:
    for phenotype, _path in GWAS_FILES:
        top_path = OUT_DIR / f"{phenotype}.top_{TOP_N}_snps.tsv"
        with top_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                add_row(records, phenotype, row, "top5000")


def scan_suggestive(records: dict[str, dict[str, object]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for phenotype, path in GWAS_FILES:
        n_suggestive = 0
        with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            column_map = resolve_columns(reader.fieldnames or [])
            for raw_row in reader:
                p_value = parse_float(raw_row[column_map["P"]])
                if p_value >= SUGGESTIVE:
                    continue
                n_suggestive += 1
                row = {
                    "SNP": raw_row[column_map["SNP"]],
                    "Allele1": raw_row[column_map["Allele1"]],
                    "Allele2": raw_row[column_map["Allele2"]],
                    "Freq1": raw_row[column_map["Freq1"]],
                    "N": raw_row[column_map["N"]],
                    "Z": raw_row[column_map["Z"]],
                    "P": raw_row[column_map["P"]],
                    "Direction": raw_row[column_map["Direction"]],
                    "HetPVal": raw_row[column_map["HetPVal"]],
                }
                add_row(records, phenotype, row, "suggestive")
        counts[phenotype] = n_suggestive
    return counts


def join_sorted(value: object) -> str:
    if isinstance(value, set):
        return ",".join(sorted(value))
    return str(value)


def write_outputs(records: dict[str, dict[str, object]], suggestive_counts: dict[str, int]) -> None:
    phenotype_columns: list[str] = []
    for phenotype, _path in GWAS_FILES:
        phenotype_columns.extend(
            [
                f"{phenotype}_P",
                f"{phenotype}_Z",
                f"{phenotype}_N",
                f"{phenotype}_Allele1",
                f"{phenotype}_Allele2",
                f"{phenotype}_Freq1",
                f"{phenotype}_Direction",
                f"{phenotype}_HetPVal",
                f"{phenotype}_palindromic",
            ]
        )

    fieldnames = [
        "SNP",
        "min_p",
        "best_phenotype",
        "include_reasons",
        "phenotypes_present",
        "genome_wide_phenotypes",
        "suggestive_phenotypes",
        "top5000_phenotypes",
        *phenotype_columns,
    ]

    rows = sorted(records.values(), key=lambda record: (float(record["min_p"]), str(record["SNP"])))
    with OUT_TSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for record in rows:
            output_row = {field: join_sorted(record.get(field, "")) for field in fieldnames}
            output_row["min_p"] = f"{float(record['min_p']):.3e}"
            writer.writerow(output_row)

    with OUT_TXT.open("w", encoding="utf-8") as handle:
        for record in rows:
            handle.write(f"{record['SNP']}\n")

    out_md = OUT_DIR / "snps_for_coordinate_mapping_summary.md"
    n_genome_wide = sum(bool(record["genome_wide_phenotypes"]) for record in records.values())
    n_shared_suggestive = sum(
        len(record["suggestive_phenotypes"]) >= 2 for record in records.values()  # type: ignore[arg-type]
    )
    lines = [
        "# GWAS SNP Universe for Coordinate Mapping",
        "",
        "## Scope",
        "",
        "This file contains the compact rsID set for coordinate lookup before locus definition. Inclusion requires P<1e-5 in at least one GWAS or membership in a phenotype-specific top-5000 table.",
        "",
        "## Counts",
        "",
        f"- Unique rsIDs selected: {len(records):,}",
        f"- rsIDs with genome-wide significance in at least one phenotype: {n_genome_wide:,}",
        f"- rsIDs suggestive in at least two phenotypes: {n_shared_suggestive:,}",
        "",
        "## Per-Phenotype Suggestive Counts Recovered",
        "",
    ]
    for phenotype, count in suggestive_counts.items():
        lines.append(f"- {phenotype}: {count:,}")
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- Table: `{OUT_TSV.relative_to(PROJECT_ROOT)}`",
            f"- rsID list: `{OUT_TXT.relative_to(PROJECT_ROOT)}`",
            "",
            "## Interpretation Gate",
            "",
            "This is not a locus table. Independent loci require genomic coordinates plus LD-aware clumping or a documented coordinate-window sensitivity analysis.",
            "",
        ]
    )
    out_md.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    records: dict[str, dict[str, object]] = {}
    load_top5000(records)
    suggestive_counts = scan_suggestive(records)
    write_outputs(records, suggestive_counts)
    print(f"Wrote {OUT_TSV}")
    print(f"Wrote {OUT_TXT}")
    print(f"Wrote {OUT_DIR / 'snps_for_coordinate_mapping_summary.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
