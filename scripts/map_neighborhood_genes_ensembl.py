#!/usr/bin/env python3
"""Map coordinate-window neighborhoods to overlapping protein-coding genes."""

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
NEIGHBORHOODS = GWAS_DIR / "coordinate_window_neighborhoods.tsv"
OUT_TSV = GWAS_DIR / "coordinate_window_neighborhood_genes.tsv"
OUT_MD = GWAS_DIR / "coordinate_window_neighborhood_genes_summary.md"
INTERIM_DIR = PROJECT_ROOT / "data" / "interim" / "gwas"
RAW_JSONL = INTERIM_DIR / "ensembl_overlap_gene_raw.jsonl"
ENSEMBL_SERVER = "https://rest.ensembl.org"
MAX_REGION_BP = 4_900_000


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threshold", default="all", choices=["all", "genome_wide", "suggestive"])
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--sleep", type=float, default=0.2)
    return parser.parse_args()


def read_neighborhoods(threshold: str, limit: int) -> list[dict[str, str]]:
    with NEIGHBORHOODS.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if threshold != "all":
        rows = [row for row in rows if row["threshold"] == threshold]
    if limit:
        rows = rows[:limit]
    return rows


def split_region(chrom: str, start: int, end: int) -> list[tuple[str, int, int]]:
    chunks = []
    cursor = start
    while cursor <= end:
        chunk_end = min(end, cursor + MAX_REGION_BP - 1)
        chunks.append((chrom, cursor, chunk_end))
        cursor = chunk_end + 1
    return chunks


def cache_key(chrom: str, start: int, end: int) -> str:
    return f"{chrom}:{start}-{end}"


def load_cache(path: Path) -> dict[str, list[dict[str, Any]]]:
    if not path.exists():
        return {}
    cache: dict[str, list[dict[str, Any]]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            obj = json.loads(line)
            key = obj.get("region")
            payload = obj.get("payload")
            if isinstance(key, str) and isinstance(payload, list):
                cache[key] = payload
    return cache


def append_cache(path: Path, key: str, payload: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps({"region": key, "payload": payload}, sort_keys=True))
        handle.write("\n")


def query_gene_overlap(chrom: str, start: int, end: int, max_attempts: int = 5) -> list[dict[str, Any]]:
    region = urllib.parse.quote(f"{chrom}:{start}-{end}", safe=":-")
    url = f"{ENSEMBL_SERVER}/overlap/region/human/{region}?feature=gene;biotype=protein_coding"
    request = urllib.request.Request(url, headers={"Accept": "application/json", "Content-Type": "application/json"})
    for attempt in range(1, max_attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                payload = json.loads(response.read().decode("utf-8"))
                return payload if isinstance(payload, list) else []
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
    raise RuntimeError("Exhausted Ensembl gene-overlap attempts")


def gene_distance_to_lead(gene: dict[str, Any], lead_pos: int) -> int:
    start = int(gene.get("start") or 0)
    end = int(gene.get("end") or 0)
    if start <= lead_pos <= end:
        return 0
    return min(abs(lead_pos - start), abs(lead_pos - end))


def best_lead_pos(row: dict[str, str]) -> int:
    prefix = f"{row['best_phenotype']}:{row['best_lead_snp']}:"
    entries = row["best_by_phenotype"].split(";")
    for entry in entries:
        if entry.startswith(prefix) and ":pos" in entry:
            return int(entry.rsplit(":pos", 1)[1])
    for entry in entries:
        if ":pos" in entry:
            return int(entry.rsplit(":pos", 1)[1])
    raise ValueError(f"Could not parse best lead position for {row['neighborhood_id']}")


def map_genes(rows: list[dict[str, str]], raw_jsonl: Path, force: bool, sleep: float) -> list[dict[str, str]]:
    cache = {} if force else load_cache(raw_jsonl)
    if force and raw_jsonl.exists():
        raw_jsonl.unlink()
    output_rows: list[dict[str, str]] = []
    seen_pairs: set[tuple[str, str]] = set()

    for index, row in enumerate(rows, start=1):
        chrom = row["chrom"]
        start = int(row["start"])
        end = int(row["end"])
        lead_pos = best_lead_pos(row)
        genes: dict[str, dict[str, Any]] = {}
        for chunk_chrom, chunk_start, chunk_end in split_region(chrom, start, end):
            key = cache_key(chunk_chrom, chunk_start, chunk_end)
            if key not in cache:
                print(f"Querying Ensembl genes {index}/{len(rows)} {key}", flush=True)
                cache[key] = query_gene_overlap(chunk_chrom, chunk_start, chunk_end)
                append_cache(raw_jsonl, key, cache[key])
                time.sleep(sleep)
            for gene in cache[key]:
                gene_id = str(gene.get("id") or gene.get("gene_id") or "")
                if gene_id:
                    genes[gene_id] = gene

        for gene_id, gene in genes.items():
            pair = (row["neighborhood_id"], gene_id)
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)
            output_rows.append(
                {
                    "neighborhood_id": row["neighborhood_id"],
                    "threshold": row["threshold"],
                    "phenotypes_present": row["phenotypes_present"],
                    "n_phenotypes": row["n_phenotypes"],
                    "neighborhood_chrom": row["chrom"],
                    "neighborhood_start": row["start"],
                    "neighborhood_end": row["end"],
                    "best_lead_snp": row["best_lead_snp"],
                    "best_lead_p": row["best_lead_p"],
                    "best_phenotype": row["best_phenotype"],
                    "gene_id": gene_id,
                    "gene_symbol": str(gene.get("external_name") or ""),
                    "gene_biotype": str(gene.get("biotype") or ""),
                    "gene_start": str(gene.get("start") or ""),
                    "gene_end": str(gene.get("end") or ""),
                    "gene_strand": str(gene.get("strand") or ""),
                    "distance_to_best_lead_bp": str(gene_distance_to_lead(gene, lead_pos)),
                    "gene_description": str(gene.get("description") or ""),
                }
            )
    output_rows.sort(
        key=lambda out: (
            out["threshold"],
            out["neighborhood_chrom"],
            int(out["neighborhood_start"]),
            int(out["distance_to_best_lead_bp"]),
            out["gene_symbol"],
        )
    )
    return output_rows


def write_outputs(rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "neighborhood_id",
        "threshold",
        "phenotypes_present",
        "n_phenotypes",
        "neighborhood_chrom",
        "neighborhood_start",
        "neighborhood_end",
        "best_lead_snp",
        "best_lead_p",
        "best_phenotype",
        "gene_id",
        "gene_symbol",
        "gene_biotype",
        "gene_start",
        "gene_end",
        "gene_strand",
        "distance_to_best_lead_bp",
        "gene_description",
    ]
    with OUT_TSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    genome_rows = [row for row in rows if row["threshold"] == "genome_wide"]
    shared_genome_rows = [row for row in genome_rows if int(row["n_phenotypes"]) >= 2]
    lines = [
        "# Coordinate-Window Neighborhood Gene Mapping",
        "",
        "## Scope",
        "",
        "Protein-coding genes overlapping coordinate-window neighborhoods were queried from Ensembl REST. This is positional evidence only; it does not identify causal genes.",
        "",
        "## Counts",
        "",
        f"- Total neighborhood-gene records: {len(rows):,}",
        f"- Genome-wide neighborhood-gene records: {len(genome_rows):,}",
        f"- Shared genome-wide neighborhood-gene records: {len(shared_genome_rows):,}",
        "",
        "## Shared Genome-Wide Neighborhood Genes Nearest Best Lead",
        "",
        "| Neighborhood | Phenotypes | Best lead | Nearest genes |",
        "|---|---|---|---|",
    ]
    by_neighborhood: dict[str, list[dict[str, str]]] = {}
    for row in shared_genome_rows:
        by_neighborhood.setdefault(row["neighborhood_id"], []).append(row)
    for neighborhood_id, neighborhood_rows in sorted(by_neighborhood.items()):
        nearest = sorted(neighborhood_rows, key=lambda row: (int(row["distance_to_best_lead_bp"]), row["gene_symbol"]))[:8]
        first = nearest[0]
        genes = "; ".join(
            f"{row['gene_symbol']}({row['distance_to_best_lead_bp']}bp)" for row in nearest if row["gene_symbol"]
        )
        lead = f"{first['best_lead_snp']}:{first['best_lead_p']}:{first['best_phenotype']}"
        lines.append(f"| {neighborhood_id} | {first['phenotypes_present']} | {lead} | {genes} |")
    lines.extend(
        [
            "",
            "## Interpretation Gate",
            "",
            "Use these genes as positional candidates for MAGMA/TWAS/eQTL/druggability integration. Do not rank therapeutic targets from positional overlap alone.",
            "",
            "## Output",
            "",
            f"- Gene table: `{OUT_TSV.relative_to(PROJECT_ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    neighborhoods = read_neighborhoods(args.threshold, args.limit)
    rows = map_genes(neighborhoods, RAW_JSONL, args.force, args.sleep)
    write_outputs(rows)
    print(f"Wrote {OUT_TSV}")
    print(f"Wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
