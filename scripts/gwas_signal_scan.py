#!/usr/bin/env python3
"""Initial canonical GWAS signal scan.

This script deliberately stops short of locus definition because the downloaded
GWAS files provide rsIDs but not chromosome/base-pair coordinates. It generates
threshold counts, lambda GC and top-SNP tables that can be used after SNP
coordinate mapping and LD clumping are added.
"""

from __future__ import annotations

import csv
import gzip
import heapq
import math
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW = PROJECT_ROOT / "data" / "raw_downloads"
OUT_DIR = PROJECT_ROOT / "results" / "gwas"
SUMMARY_TSV = OUT_DIR / "gwas_signal_summary.tsv"

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
GENOME_WIDE = 5e-8
SUGGESTIVE = 1e-5
LAMBDA_MEDIAN_CHISQ_1DF = 0.454936423119572


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


def median(values: list[float]) -> float:
    values.sort()
    n = len(values)
    mid = n // 2
    if n % 2:
        return values[mid]
    return (values[mid - 1] + values[mid]) / 2


def scan_file(phenotype: str, path: Path) -> dict[str, str]:
    top_heap: list[tuple[float, int, dict[str, str]]] = []
    chisq_values: list[float] = []
    n_rows = 0
    n_gws = 0
    n_suggestive = 0
    min_p = 1.0
    min_p_snp = ""

    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        column_map = resolve_columns(reader.fieldnames or [])

        for row in reader:
            n_rows += 1
            p_value = parse_float(row[column_map["P"]])
            z_value = parse_float(row[column_map["Z"]])
            chisq_values.append(z_value * z_value)

            if p_value < min_p:
                min_p = p_value
                min_p_snp = row[column_map["SNP"]]
            if p_value < GENOME_WIDE:
                n_gws += 1
            if p_value < SUGGESTIVE:
                n_suggestive += 1

            record = {
                "phenotype": phenotype,
                "SNP": row[column_map["SNP"]],
                "Allele1": row[column_map["Allele1"]],
                "Allele2": row[column_map["Allele2"]],
                "Freq1": row[column_map["Freq1"]],
                "N": row[column_map["N"]],
                "Z": row[column_map["Z"]],
                "P": row[column_map["P"]],
                "Direction": row[column_map["Direction"]],
                "HetPVal": row[column_map["HetPVal"]],
            }
            rank_key = -p_value
            if len(top_heap) < TOP_N:
                heapq.heappush(top_heap, (rank_key, n_rows, record))
            elif rank_key > top_heap[0][0]:
                heapq.heapreplace(top_heap, (rank_key, n_rows, record))

    lambda_gc = median(chisq_values) / LAMBDA_MEDIAN_CHISQ_1DF
    top_records = [item[2] for item in sorted(top_heap, key=lambda item: float(item[2]["P"]))]
    top_path = OUT_DIR / f"{phenotype}.top_{TOP_N}_snps.tsv"
    with top_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(top_records[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(top_records)

    threshold_path = OUT_DIR / f"{phenotype}.threshold_counts.tsv"
    with threshold_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["phenotype", "threshold", "count"], delimiter="\t")
        writer.writeheader()
        writer.writerow({"phenotype": phenotype, "threshold": "P<5e-8", "count": n_gws})
        writer.writerow({"phenotype": phenotype, "threshold": "P<1e-5", "count": n_suggestive})

    return {
        "phenotype": phenotype,
        "rows": str(n_rows),
        "genome_wide_snps_p_lt_5e_8": str(n_gws),
        "suggestive_snps_p_lt_1e_5": str(n_suggestive),
        "lambda_gc_from_z": f"{lambda_gc:.4f}",
        "min_p": f"{min_p:.3e}",
        "min_p_snp": min_p_snp,
        "top_snps_file": str(top_path.relative_to(PROJECT_ROOT)),
        "threshold_counts_file": str(threshold_path.relative_to(PROJECT_ROOT)),
    }


def write_markdown(summary_rows: list[dict[str, str]]) -> None:
    out_md = OUT_DIR / "gwas_signal_summary.md"
    lines = [
        "# GWAS Canonical Signal Scan",
        "",
        "## Scope",
        "",
        "This scan summarizes genome-wide and suggestive SNP-level evidence without defining loci. Locus discovery requires SNP genomic coordinates and an LD reference panel.",
        "",
        "## Summary",
        "",
        "| Phenotype | Rows | P<5e-8 SNPs | P<1e-5 SNPs | Lambda GC | Min P | Top SNP |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {row['phenotype']} | {row['rows']} | {row['genome_wide_snps_p_lt_5e_8']} | {row['suggestive_snps_p_lt_1e_5']} | {row['lambda_gc_from_z']} | {row['min_p']} | {row['min_p_snp']} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation Gate",
            "",
            "These outputs are suitable for QC and planning. They are not yet independent loci, mapped genes or therapeutic targets.",
            "",
        ]
    )
    out_md.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_rows = [scan_file(phenotype, path) for phenotype, path in GWAS_FILES]
    with SUMMARY_TSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    write_markdown(summary_rows)
    print(f"Wrote {SUMMARY_TSV}")
    print(f"Wrote {OUT_DIR / 'gwas_signal_summary.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
