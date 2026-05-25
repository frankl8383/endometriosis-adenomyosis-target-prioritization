#!/usr/bin/env python3
"""Structural QC for phase 1 GWAS summary statistics.

This is a gate before genetics-first candidate nomination. It checks that the
downloaded files are complete gzip streams and that the expected Koller et al.
summary-statistic columns are present and numerically valid.
"""

from __future__ import annotations

import csv
import gzip
import math
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MANIFEST = PROJECT_ROOT / "data" / "download_manifest.tsv"
RAW = PROJECT_ROOT / "data" / "raw_downloads"
OUT_TSV = PROJECT_ROOT / "results" / "phase1_gwas_audit.tsv"
OUT_MD = PROJECT_ROOT / "results" / "phase1_gwas_audit.md"

REQUIRED_COLUMNS = ["SNP", "Allele1", "Allele2", "Freq1", "N", "Z", "P", "DIRECTION", "HetPVal"]
COLUMN_ALIASES = {
    "SNP": ["SNP"],
    "Allele1": ["Allele1"],
    "Allele2": ["Allele2"],
    "Freq1": ["Freq1"],
    "N": ["N"],
    "Z": ["Z"],
    "P": ["P", "P.value"],
    "DIRECTION": ["DIRECTION", "Direction"],
    "HetPVal": ["HetPVal"],
}
PALINDROMIC = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}


def safe_filename(resource: str, file_or_record: str) -> str:
    resource_safe = resource.replace(" ", "_").replace("/", "_")
    return f"{resource_safe}__{file_or_record}"


def parse_float(value: str) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def resolve_columns(header: list[str]) -> tuple[dict[str, str], list[str]]:
    mapping: dict[str, str] = {}
    missing: list[str] = []
    header_set = set(header)
    for logical_name, aliases in COLUMN_ALIASES.items():
        match = next((alias for alias in aliases if alias in header_set), None)
        if match is None:
            missing.append(logical_name)
        else:
            mapping[logical_name] = match
    return mapping, missing


def audit_file(label: str, path: Path) -> dict[str, str]:
    partial = path.with_suffix(path.suffix + ".partial")
    if not path.exists():
        return {
            "file": label,
            "path": str(path),
            "file_state": "partial_file_only" if partial.exists() else "missing",
            "gzip_state": "not_checked",
            "required_columns_present": "no",
            "missing_columns": ",".join(REQUIRED_COLUMNS),
            "rows": "0",
            "unique_snps": "0",
            "duplicate_snps": "0",
            "invalid_p": "0",
            "invalid_z": "0",
            "invalid_n": "0",
            "invalid_freq": "0",
            "missing_snp": "0",
            "missing_allele": "0",
            "palindromic_snps": "0",
            "min_p": "",
            "max_p": "",
            "min_n": "",
            "max_n": "",
            "status": "not_ready",
        }

    seen_snps: set[str] = set()
    duplicate_snps = 0
    invalid_p = invalid_z = invalid_n = invalid_freq = 0
    missing_snp = missing_allele = palindromic = 0
    min_p: float | None = None
    max_p: float | None = None
    min_n: float | None = None
    max_n: float | None = None
    rows = 0

    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = reader.fieldnames or []
            column_map, missing_columns = resolve_columns(header)
            if missing_columns:
                return {
                    "file": label,
                    "path": str(path),
                    "file_state": "complete_file_present",
                    "gzip_state": "ok",
                    "required_columns_present": "no",
                    "missing_columns": ",".join(missing_columns),
                    "rows": "0",
                    "unique_snps": "0",
                    "duplicate_snps": "0",
                    "invalid_p": "0",
                    "invalid_z": "0",
                    "invalid_n": "0",
                    "invalid_freq": "0",
                    "missing_snp": "0",
                    "missing_allele": "0",
                    "palindromic_snps": "0",
                    "min_p": "",
                    "max_p": "",
                    "min_n": "",
                    "max_n": "",
                    "status": "not_ready",
                }

            for row in reader:
                rows += 1
                snp = (row.get(column_map["SNP"]) or "").strip()
                if not snp:
                    missing_snp += 1
                elif snp in seen_snps:
                    duplicate_snps += 1
                else:
                    seen_snps.add(snp)

                allele1 = (row.get(column_map["Allele1"]) or "").upper().strip()
                allele2 = (row.get(column_map["Allele2"]) or "").upper().strip()
                if not allele1 or not allele2:
                    missing_allele += 1
                elif (allele1, allele2) in PALINDROMIC:
                    palindromic += 1

                p_value = parse_float(row.get(column_map["P"], ""))
                if p_value is None or p_value < 0 or p_value > 1:
                    invalid_p += 1
                else:
                    min_p = p_value if min_p is None else min(min_p, p_value)
                    max_p = p_value if max_p is None else max(max_p, p_value)

                z_value = parse_float(row.get(column_map["Z"], ""))
                if z_value is None:
                    invalid_z += 1

                n_value = parse_float(row.get(column_map["N"], ""))
                if n_value is None or n_value <= 0:
                    invalid_n += 1
                else:
                    min_n = n_value if min_n is None else min(min_n, n_value)
                    max_n = n_value if max_n is None else max(max_n, n_value)

                freq = parse_float(row.get(column_map["Freq1"], ""))
                if freq is None or freq < 0 or freq > 1:
                    invalid_freq += 1

    except Exception as exc:  # noqa: BLE001
        return {
            "file": label,
            "path": str(path),
            "file_state": "complete_file_present",
            "gzip_state": f"fail:{type(exc).__name__}:{exc}",
            "required_columns_present": "unknown",
            "missing_columns": "",
            "rows": str(rows),
            "unique_snps": str(len(seen_snps)),
            "duplicate_snps": str(duplicate_snps),
            "invalid_p": str(invalid_p),
            "invalid_z": str(invalid_z),
            "invalid_n": str(invalid_n),
            "invalid_freq": str(invalid_freq),
            "missing_snp": str(missing_snp),
            "missing_allele": str(missing_allele),
            "palindromic_snps": str(palindromic),
            "min_p": "" if min_p is None else f"{min_p:.3e}",
            "max_p": "" if max_p is None else f"{max_p:.3e}",
            "min_n": "" if min_n is None else f"{min_n:.0f}",
            "max_n": "" if max_n is None else f"{max_n:.0f}",
            "status": "not_ready",
        }

    critical_issues = invalid_p + invalid_z + invalid_n + invalid_freq + missing_snp + missing_allele
    if rows > 1_000_000 and critical_issues == 0:
        status = "ready_duplicate_handling_required" if duplicate_snps else "ready"
    else:
        status = "review"

    return {
        "file": label,
        "path": str(path),
        "file_state": "complete_file_present",
        "gzip_state": "ok",
        "required_columns_present": "yes",
        "missing_columns": "",
        "rows": str(rows),
        "unique_snps": str(len(seen_snps)),
        "duplicate_snps": str(duplicate_snps),
        "invalid_p": str(invalid_p),
        "invalid_z": str(invalid_z),
        "invalid_n": str(invalid_n),
        "invalid_freq": str(invalid_freq),
        "missing_snp": str(missing_snp),
        "missing_allele": str(missing_allele),
        "palindromic_snps": str(palindromic),
        "min_p": "" if min_p is None else f"{min_p:.3e}",
        "max_p": "" if max_p is None else f"{max_p:.3e}",
        "min_n": "" if min_n is None else f"{min_n:.0f}",
        "max_n": "" if max_n is None else f"{max_n:.0f}",
        "status": status,
    }


def manifest_gwas_rows() -> list[dict[str, str]]:
    with MANIFEST.open("r", encoding="utf-8") as handle:
        return [row for row in csv.DictReader(handle, delimiter="\t") if row["download_phase"] == "phase1_gwas"]


def write_markdown(rows: list[dict[str, str]]) -> None:
    ready_states = {"ready", "ready_duplicate_handling_required"}
    ready = sum(row["status"] in ready_states for row in rows)
    lines = [
        "# Phase 1 GWAS Structural QC",
        "",
        "## Summary",
        "",
        f"- GWAS files ready: {ready}/{len(rows)}",
        "",
        "## Files",
        "",
        "| File | State | Gzip | Rows | Unique SNPs | Duplicates | Invalid P/Z/N/Freq | Palindromic | Status |",
        "|---|---|---|---:|---:|---:|---|---:|---|",
    ]
    for row in rows:
        invalid = f"{row['invalid_p']}/{row['invalid_z']}/{row['invalid_n']}/{row['invalid_freq']}"
        lines.append(
            f"| `{row['file']}` | {row['file_state']} | {row['gzip_state']} | {row['rows']} | {row['unique_snps']} | {row['duplicate_snps']} | {invalid} | {row['palindromic_snps']} | {row['status']} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation Gate",
            "",
            "All three core GWAS files must be `ready` or `ready_duplicate_handling_required` before lead-SNP/locus definition, MAGMA, or cross-disease genetic comparison. Palindromic SNPs are not an exclusion at this stage, but they must be handled during allele harmonization.",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    rows: list[dict[str, str]] = []
    for row in manifest_gwas_rows():
        file_name = safe_filename(row["resource"], row["file_or_record"])
        rows.append(audit_file(row["file_or_record"], RAW / file_name))

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    write_markdown(rows)
    print(f"Wrote {OUT_TSV}")
    print(f"Wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
