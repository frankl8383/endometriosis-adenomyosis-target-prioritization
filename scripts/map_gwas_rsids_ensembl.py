#!/usr/bin/env python3
"""Map selected GWAS rsIDs to GRCh38 coordinates using Ensembl REST.

Input:
    results/gwas/snps_for_coordinate_mapping.txt

Outputs:
    results/gwas/ensembl_grch38_variant_coordinates.tsv
    results/gwas/ensembl_grch38_variant_coordinates_summary.md
    data/interim/gwas/ensembl_variation_raw.jsonl

The Ensembl POST variation endpoint accepts up to 200 IDs per request. The raw
JSONL cache lets the mapping table be rebuilt without re-querying the API.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_RSIDS = PROJECT_ROOT / "results" / "gwas" / "snps_for_coordinate_mapping.txt"
OUT_DIR = PROJECT_ROOT / "results" / "gwas"
INTERIM_DIR = PROJECT_ROOT / "data" / "interim" / "gwas"
RAW_JSONL = INTERIM_DIR / "ensembl_variation_raw.jsonl"
OUT_TSV = OUT_DIR / "ensembl_grch38_variant_coordinates.tsv"
OUT_MD = OUT_DIR / "ensembl_grch38_variant_coordinates_summary.md"

ENSEMBL_URL = "https://rest.ensembl.org/variation/homo_sapiens"
MAX_BATCH_SIZE = 200
PRIMARY_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT_RSIDS, help="Text file with one rsID per line.")
    parser.add_argument("--raw-jsonl", type=Path, default=RAW_JSONL, help="Raw Ensembl response cache.")
    parser.add_argument("--out-tsv", type=Path, default=OUT_TSV, help="Coordinate table output.")
    parser.add_argument("--batch-size", type=int, default=MAX_BATCH_SIZE, help="Batch size, maximum 200.")
    parser.add_argument("--sleep", type=float, default=0.2, help="Seconds to sleep between API requests.")
    parser.add_argument("--limit", type=int, default=0, help="Limit rsIDs for testing; 0 means all.")
    parser.add_argument("--force", action="store_true", help="Ignore existing raw cache and re-query all IDs.")
    return parser.parse_args()


def read_rsids(path: Path, limit: int = 0) -> list[str]:
    rsids = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    seen: set[str] = set()
    unique = []
    for rsid in rsids:
        if rsid in seen:
            continue
        seen.add(rsid)
        unique.append(rsid)
        if limit and len(unique) >= limit:
            break
    return unique


def batched(items: list[str], batch_size: int) -> list[list[str]]:
    return [items[i : i + batch_size] for i in range(0, len(items), batch_size)]


def load_raw_cache(path: Path) -> dict[str, dict[str, Any]]:
    if not path.exists():
        return {}
    cache: dict[str, dict[str, Any]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            obj = json.loads(line)
            rsid = obj.get("id")
            if isinstance(rsid, str):
                cache[rsid] = obj
    return cache


def query_ensembl(ids: list[str], max_attempts: int = 5) -> dict[str, Any]:
    body = json.dumps({"ids": ids}).encode("utf-8")
    request = urllib.request.Request(
        ENSEMBL_URL,
        data=body,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    for attempt in range(1, max_attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            if exc.code in {429, 500, 502, 503, 504} and attempt < max_attempts:
                time.sleep(min(30, 2**attempt))
                continue
            raise
        except urllib.error.URLError:
            if attempt < max_attempts:
                time.sleep(min(30, 2**attempt))
                continue
            raise
    raise RuntimeError("Exhausted Ensembl query attempts")


def append_cache(path: Path, results: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        for rsid, payload in results.items():
            handle.write(json.dumps({"id": rsid, "payload": payload}, sort_keys=True))
            handle.write("\n")


def grch38_chromosome_mappings(payload: dict[str, Any]) -> list[dict[str, Any]]:
    mappings = payload.get("mappings") or []
    if not isinstance(mappings, list):
        return []
    selected = []
    for mapping in mappings:
        if not isinstance(mapping, dict):
            continue
        if mapping.get("assembly_name") != "GRCh38":
            continue
        if mapping.get("coord_system") != "chromosome":
            continue
        if str(mapping.get("seq_region_name")) not in PRIMARY_CHROMS:
            continue
        selected.append(mapping)
    return selected


def choose_mapping(mappings: list[dict[str, Any]]) -> tuple[str, dict[str, Any] | None]:
    if not mappings:
        return "unmapped_primary_grch38_chromosome", None
    unique_locations = {
        (
            str(mapping.get("seq_region_name")),
            str(mapping.get("start")),
            str(mapping.get("end")),
            str(mapping.get("strand")),
        )
        for mapping in mappings
    }
    if len(unique_locations) == 1:
        return "mapped_unique_primary_grch38_chromosome", mappings[0]
    return "mapped_multiple_primary_grch38_chromosome", mappings[0]


def normalize_payload(obj: dict[str, Any]) -> dict[str, Any] | None:
    payload = obj.get("payload")
    if isinstance(payload, dict):
        return payload
    return None


def payload_field(payload: dict[str, Any] | None, field: str) -> str:
    if not payload:
        return ""
    value = payload.get(field)
    if value is None:
        return ""
    if isinstance(value, list):
        return ",".join(str(item) for item in value)
    return str(value)


def build_coordinate_rows(rsids: list[str], cache: dict[str, dict[str, Any]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for rsid in rsids:
        obj = cache.get(rsid)
        payload = normalize_payload(obj) if obj else None
        if payload is None:
            rows.append(
                {
                    "SNP": rsid,
                    "mapping_status": "not_returned_by_ensembl",
                    "chrom": "",
                    "start": "",
                    "end": "",
                    "location": "",
                    "strand": "",
                    "allele_string": "",
                    "assembly_name": "",
                    "coord_system": "",
                    "var_class": "",
                    "minor_allele": "",
                    "MAF": "",
                    "most_severe_consequence": "",
                    "ambiguity": "",
                    "synonyms": "",
                    "n_total_mappings": "0",
                    "n_primary_grch38_chromosome_mappings": "0",
                    "primary_grch38_locations": "",
                }
            )
            continue

        primary_mappings = grch38_chromosome_mappings(payload)
        status, chosen = choose_mapping(primary_mappings)
        all_mappings = payload.get("mappings") if isinstance(payload.get("mappings"), list) else []
        primary_locations = sorted(
            {
                f"{mapping.get('seq_region_name')}:{mapping.get('start')}-{mapping.get('end')}:{mapping.get('strand')}"
                for mapping in primary_mappings
            }
        )
        rows.append(
            {
                "SNP": rsid,
                "mapping_status": status,
                "chrom": "" if chosen is None else str(chosen.get("seq_region_name", "")),
                "start": "" if chosen is None else str(chosen.get("start", "")),
                "end": "" if chosen is None else str(chosen.get("end", "")),
                "location": "" if chosen is None else str(chosen.get("location", "")),
                "strand": "" if chosen is None else str(chosen.get("strand", "")),
                "allele_string": "" if chosen is None else str(chosen.get("allele_string", "")),
                "assembly_name": "" if chosen is None else str(chosen.get("assembly_name", "")),
                "coord_system": "" if chosen is None else str(chosen.get("coord_system", "")),
                "var_class": payload_field(payload, "var_class"),
                "minor_allele": payload_field(payload, "minor_allele"),
                "MAF": payload_field(payload, "MAF"),
                "most_severe_consequence": payload_field(payload, "most_severe_consequence"),
                "ambiguity": payload_field(payload, "ambiguity"),
                "synonyms": payload_field(payload, "synonyms"),
                "n_total_mappings": str(len(all_mappings)),
                "n_primary_grch38_chromosome_mappings": str(len(primary_mappings)),
                "primary_grch38_locations": ",".join(primary_locations),
            }
        )
    return rows


def chrom_sort_key(row: dict[str, str]) -> tuple[int, int, str]:
    chrom = row["chrom"]
    order = {str(i): i for i in range(1, 23)} | {"X": 23, "Y": 24, "MT": 25}
    pos = int(row["start"]) if row["start"].isdigit() else 0
    return (order.get(chrom, 99), pos, row["SNP"])


def write_coordinate_table(rows: list[dict[str, str]], out_tsv: Path) -> None:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "SNP",
        "mapping_status",
        "chrom",
        "start",
        "end",
        "location",
        "strand",
        "allele_string",
        "assembly_name",
        "coord_system",
        "var_class",
        "minor_allele",
        "MAF",
        "most_severe_consequence",
        "ambiguity",
        "synonyms",
        "n_total_mappings",
        "n_primary_grch38_chromosome_mappings",
        "primary_grch38_locations",
    ]
    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(rows, key=chrom_sort_key))


def write_summary(rows: list[dict[str, str]], out_md: Path, raw_jsonl: Path) -> None:
    status_counts: dict[str, int] = {}
    for row in rows:
        status_counts[row["mapping_status"]] = status_counts.get(row["mapping_status"], 0) + 1
    mapped = sum(
        count for status, count in status_counts.items() if status.startswith("mapped_")
    )
    lines = [
        "# Ensembl GRCh38 Coordinate Mapping",
        "",
        "## Source",
        "",
        "rsIDs were queried against Ensembl REST `POST /variation/homo_sapiens`, retaining primary GRCh38 chromosome mappings only.",
        "",
        "## Mapping Status",
        "",
        f"- Input rsIDs: {len(rows):,}",
        f"- Mapped to primary GRCh38 chromosomes: {mapped:,}",
    ]
    for status, count in sorted(status_counts.items()):
        lines.append(f"- {status}: {count:,}")
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- Coordinate table: `{OUT_TSV.relative_to(PROJECT_ROOT)}`",
            f"- Raw Ensembl JSONL cache: `{raw_jsonl.relative_to(PROJECT_ROOT)}`",
            "",
            "## Interpretation Gate",
            "",
            "These coordinates enable coordinate-window locus construction and preparation for LD-aware clumping. They do not by themselves prove independence of association signals.",
            "",
        ]
    )
    out_md.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    if args.batch_size > MAX_BATCH_SIZE:
        raise ValueError("Ensembl variation POST batch size must be <= 200")

    rsids = read_rsids(args.input, args.limit)
    cache = {} if args.force else load_raw_cache(args.raw_jsonl)
    cached_ids = set(cache)
    missing = [rsid for rsid in rsids if rsid not in cached_ids]

    if missing:
        args.raw_jsonl.parent.mkdir(parents=True, exist_ok=True)
        if args.force and args.raw_jsonl.exists():
            args.raw_jsonl.unlink()
            cache = {}
        for index, batch in enumerate(batched(missing, args.batch_size), start=1):
            print(f"Querying Ensembl batch {index}: {len(batch)} IDs", flush=True)
            results = query_ensembl(batch)
            append_cache(args.raw_jsonl, results)
            for rsid, payload in results.items():
                cache[rsid] = {"id": rsid, "payload": payload}
            time.sleep(args.sleep)

    rows = build_coordinate_rows(rsids, cache)
    write_coordinate_table(rows, args.out_tsv)
    write_summary(rows, OUT_MD, args.raw_jsonl)
    print(f"Wrote {args.out_tsv}")
    print(f"Wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
