#!/usr/bin/env python3
"""Download GEO family SOFT files with sample-level metadata."""

from __future__ import annotations

import gzip
import urllib.request
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_ROOT / "data" / "raw_metadata"
USER_AGENT = "endo-adeno-target-prioritization/0.1"

GEO_SOFT_URLS = {
    "GSE179640": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179640/soft/GSE179640_family.soft.gz",
    "GSE234354": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234354/soft/GSE234354_family.soft.gz",
    "GSE313775": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313775/soft/GSE313775_family.soft.gz",
    "GSE51981": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51981/soft/GSE51981_family.soft.gz",
    "GSE141549": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141549/soft/GSE141549_family.soft.gz",
}


def download(url: str, destination: Path) -> None:
    if destination.exists() and destination.stat().st_size > 0:
        print(f"SKIP {destination.name}")
        return
    partial = destination.with_suffix(destination.suffix + ".partial")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=90) as response, partial.open("wb") as handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)
            print(f"  {destination.name}: {partial.stat().st_size} bytes", flush=True)
    partial.replace(destination)
    print(f"DOWNLOADED {destination.name} {destination.stat().st_size} bytes")


def decompress(source: Path, destination: Path) -> None:
    if destination.exists() and destination.stat().st_size > 0:
        print(f"SKIP {destination.name}")
        return
    with gzip.open(source, "rt", encoding="utf-8", errors="replace") as gz_handle:
        destination.write_text(gz_handle.read(), encoding="utf-8")
    print(f"WROTE {destination.name} {destination.stat().st_size} bytes")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for accession, url in GEO_SOFT_URLS.items():
        gz_path = OUT_DIR / f"{accession}_family.soft.gz"
        txt_path = OUT_DIR / f"{accession}_family.soft.txt"
        download(url, gz_path)
        decompress(gz_path, txt_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
