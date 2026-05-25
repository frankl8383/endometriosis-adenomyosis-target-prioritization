#!/usr/bin/env python3
"""Pairwise LD sensitivity for cross-phenotype coordinate-window lead SNPs.

This script queries Ensembl REST pairwise LD for lead-SNP pairs from
`coordinate_window_locus_overlaps.tsv`. It is a sensitivity layer for shared
genomic neighborhoods, not a replacement for PLINK clumping with a local LD
reference panel.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GWAS_DIR = PROJECT_ROOT / "results" / "gwas"
OVERLAPS = GWAS_DIR / "coordinate_window_locus_overlaps.tsv"
OUT_TSV = GWAS_DIR / "ensembl_pairwise_ld_sensitivity.tsv"
OUT_MD = GWAS_DIR / "ensembl_pairwise_ld_sensitivity_summary.md"
INTERIM_DIR = PROJECT_ROOT / "data" / "interim" / "gwas"
RAW_JSONL = INTERIM_DIR / "ensembl_pairwise_ld_raw.jsonl"

ENSEMBL_SERVER = "https://rest.ensembl.org"
DEFAULT_POPULATION = "1000GENOMES:phase_3:EUR"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threshold", default="all", choices=["all", "genome_wide", "suggestive"])
    parser.add_argument("--population", default=DEFAULT_POPULATION)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--sleep", type=float, default=0.2)
    parser.add_argument("--limit", type=int, default=0)
    return parser.parse_args()


def read_overlaps(threshold: str, limit: int) -> list[dict[str, str]]:
    with OVERLAPS.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows = [row for row in rows if row["relationship"] == "overlapping_windows"]
    if threshold != "all":
        rows = [row for row in rows if row["threshold"] == threshold]
    if limit:
        rows = rows[:limit]
    return rows


def pair_key(snp1: str, snp2: str, population: str) -> str:
    ordered = sorted([snp1, snp2])
    return f"{ordered[0]}|{ordered[1]}|{population}"


def load_cache(path: Path) -> dict[str, dict[str, Any]]:
    if not path.exists():
        return {}
    cache: dict[str, dict[str, Any]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            obj = json.loads(line)
            key = obj.get("key")
            if isinstance(key, str):
                cache[key] = obj
    return cache


def append_cache(path: Path, key: str, payload: Any, status: str, error: str = "") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps({"key": key, "payload": payload, "status": status, "error": error}, sort_keys=True))
        handle.write("\n")


def query_pairwise_ld(snp1: str, snp2: str, population: str, max_attempts: int = 5) -> list[dict[str, Any]]:
    encoded_population = urllib.parse.quote(population)
    url = f"{ENSEMBL_SERVER}/ld/human/pairwise/{snp1}/{snp2}?population_name={encoded_population}"
    request = urllib.request.Request(url, headers={"Accept": "application/json", "Content-Type": "application/json"})
    for attempt in range(1, max_attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                payload = json.loads(response.read().decode("utf-8"))
                return payload if isinstance(payload, list) else []
        except urllib.error.HTTPError as exc:
            body = exc.read().decode("utf-8", errors="replace")
            if exc.code in {429, 500, 502, 503, 504} and attempt < max_attempts:
                time.sleep(min(30, 2**attempt))
                continue
            raise RuntimeError(f"HTTP {exc.code}: {body}") from exc
        except urllib.error.URLError as exc:
            if attempt < max_attempts:
                time.sleep(min(30, 2**attempt))
                continue
            raise RuntimeError(str(exc)) from exc
    raise RuntimeError("Exhausted Ensembl LD attempts")


def parse_ld_payload(payload: Any, snp1: str, snp2: str) -> tuple[str, str, str]:
    if snp1 == snp2:
        return "1", "1", "same_snp"
    if not isinstance(payload, list) or not payload:
        return "", "", "not_available"
    first = payload[0]
    r2 = str(first.get("r2", ""))
    d_prime = str(first.get("d_prime", ""))
    return r2, d_prime, "available" if r2 else "not_available"


def classify_r2(r2_text: str, status: str) -> str:
    if status == "same_snp":
        return "same_snp"
    if status != "available" or not r2_text:
        return status
    r2 = float(r2_text)
    if r2 >= 0.8:
        return "high_ld_r2_ge_0.8"
    if r2 >= 0.2:
        return "moderate_ld_r2_0.2_to_0.8"
    return "low_ld_r2_lt_0.2"


def build_rows(overlap_rows: list[dict[str, str]], population: str, force: bool, sleep: float) -> list[dict[str, str]]:
    cache = {} if force else load_cache(RAW_JSONL)
    if force and RAW_JSONL.exists():
        RAW_JSONL.unlink()
    out_rows: list[dict[str, str]] = []

    for index, row in enumerate(overlap_rows, start=1):
        snp1 = row["lead_snp_a"]
        snp2 = row["lead_snp_b"]
        key = pair_key(snp1, snp2, population)

        if snp1 == snp2:
            payload: Any = [{"variation1": snp1, "variation2": snp2, "population_name": population, "r2": "1", "d_prime": "1"}]
            raw_status = "same_snp"
            error = ""
        elif key in cache:
            payload = cache[key].get("payload")
            raw_status = str(cache[key].get("status", "cached"))
            error = str(cache[key].get("error", ""))
        else:
            print(f"Querying Ensembl LD {index}/{len(overlap_rows)} {snp1} {snp2}", flush=True)
            try:
                payload = query_pairwise_ld(snp1, snp2, population)
                raw_status = "queried"
                error = ""
            except RuntimeError as exc:
                payload = []
                raw_status = "error"
                error = str(exc)
            append_cache(RAW_JSONL, key, payload, raw_status, error)
            cache[key] = {"payload": payload, "status": raw_status, "error": error}
            time.sleep(sleep)

        r2, d_prime, ld_status = parse_ld_payload(payload, snp1, snp2)
        if raw_status == "error":
            ld_status = "api_error"
        out_rows.append(
            {
                "threshold": row["threshold"],
                "population": population,
                "phenotype_a": row["phenotype_a"],
                "locus_a": row["locus_a"],
                "lead_snp_a": snp1,
                "lead_p_a": row["lead_p_a"],
                "phenotype_b": row["phenotype_b"],
                "locus_b": row["locus_b"],
                "lead_snp_b": snp2,
                "lead_p_b": row["lead_p_b"],
                "chrom": row["chrom"],
                "lead_distance_bp": row["lead_distance_bp"],
                "window_overlap_bp": row["window_overlap_bp"],
                "r2": r2,
                "d_prime": d_prime,
                "ld_status": ld_status,
                "ld_class": classify_r2(r2, ld_status),
                "query_status": raw_status,
                "query_error": error,
            }
        )
    out_rows.sort(key=lambda item: (item["threshold"], item["chrom"], int(item["lead_distance_bp"]), item["lead_snp_a"], item["lead_snp_b"]))
    return out_rows


def write_outputs(rows: list[dict[str, str]], population: str) -> None:
    fieldnames = [
        "threshold",
        "population",
        "phenotype_a",
        "locus_a",
        "lead_snp_a",
        "lead_p_a",
        "phenotype_b",
        "locus_b",
        "lead_snp_b",
        "lead_p_b",
        "chrom",
        "lead_distance_bp",
        "window_overlap_bp",
        "r2",
        "d_prime",
        "ld_status",
        "ld_class",
        "query_status",
        "query_error",
    ]
    with OUT_TSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# Ensembl Pairwise LD Sensitivity",
        "",
        "## Scope",
        "",
        f"Lead-SNP pairs from overlapping coordinate-window loci were queried with Ensembl REST pairwise LD using `{population}`. This is a sensitivity analysis for cross-phenotype lead-pair relatedness, not full LD clumping.",
        "",
        "## LD Class Counts",
        "",
        "| Threshold | LD class | Count |",
        "|---|---|---:|",
    ]
    for threshold in sorted({row["threshold"] for row in rows}):
        subset = [row for row in rows if row["threshold"] == threshold]
        for ld_class in sorted({row["ld_class"] for row in subset}):
            count = sum(row["ld_class"] == ld_class for row in subset)
            lines.append(f"| {threshold} | {ld_class} | {count:,} |")

    lines.extend(
        [
            "",
            "## Genome-Wide Lead-Pair LD",
            "",
            "| Pair | Chr | Lead A | Lead B | Distance bp | r2 | D' | LD class |",
            "|---|---|---|---|---:|---:|---:|---|",
        ]
    )
    for row in sorted([row for row in rows if row["threshold"] == "genome_wide"], key=lambda item: float(item["r2"] or -1), reverse=True):
        pair = f"{row['phenotype_a']} vs {row['phenotype_b']}"
        lines.append(
            f"| {pair} | {row['chrom']} | {row['lead_snp_a']} | {row['lead_snp_b']} | {row['lead_distance_bp']} | {row['r2'] or 'NA'} | {row['d_prime'] or 'NA'} | {row['ld_class']} |"
        )

    lines.extend(
        [
            "",
            "## Interpretation Gate",
            "",
            "High or moderate pairwise LD strengthens the case that overlapping coordinate windows may reflect related genetic signals. Low or unavailable LD means the shared neighborhood must be treated as potentially distinct signals until local PLINK-based clumping/fine-mapping is done.",
            "",
            "## Outputs",
            "",
            f"- LD table: `{OUT_TSV.relative_to(PROJECT_ROOT)}`",
            f"- Raw cache: `{RAW_JSONL.relative_to(PROJECT_ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    overlap_rows = read_overlaps(args.threshold, args.limit)
    rows = build_rows(overlap_rows, args.population, args.force, args.sleep)
    write_outputs(rows, args.population)
    print(f"Wrote {OUT_TSV}")
    print(f"Wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
