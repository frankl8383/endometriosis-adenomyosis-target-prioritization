#!/usr/bin/env python3
"""Download smaller single-cell support files before large atlas objects."""

from __future__ import annotations

import urllib.request
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_ROOT / "data" / "raw_downloads"
USER_AGENT = "endo-adeno-target-prioritization/0.1"

FILES = {
    "GSE203191__GSE203191_Shih_endo_meta_frame.tsv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE203nnn/GSE203191/suppl/GSE203191_Shih_endo_meta_frame.tsv.gz",
    "GSE203191__GSE203191_metadata.xlsx": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE203nnn/GSE203191/suppl/GSE203191_metadata.xlsx",
    "GSE179640__GSE179640_organoid_cellhashing_mapping.xls.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179640/suppl/GSE179640_organoid_cellhashing_mapping.xls.gz",
}


def download(name: str, url: str) -> None:
    destination = OUT_DIR / name
    if destination.exists() and destination.stat().st_size > 0:
        print(f"SKIP {name} {destination.stat().st_size} bytes")
        return
    partial = destination.with_suffix(destination.suffix + ".partial")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=90) as response, partial.open("wb") as handle:
        total = int(response.headers.get("Content-Length", "0") or "0")
        written = 0
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)
            written += len(chunk)
            print(f"{name}: {written}/{total or '?'} bytes", flush=True)
            if total and written >= total:
                break
    partial.replace(destination)
    print(f"DOWNLOADED {name} {destination.stat().st_size} bytes")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        download(name, url)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
