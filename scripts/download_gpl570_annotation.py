#!/usr/bin/env python3
"""Download the GEO GPL570 probe annotation file."""

from __future__ import annotations

import urllib.request
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_ROOT / "data" / "raw_metadata"
URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"
DESTINATION = OUT_DIR / "GPL570.annot.gz"
USER_AGENT = "endo-adeno-target-prioritization/0.1"


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if DESTINATION.exists() and DESTINATION.stat().st_size > 0:
        print(f"SKIP {DESTINATION} {DESTINATION.stat().st_size} bytes")
        return 0

    partial = DESTINATION.with_suffix(DESTINATION.suffix + ".partial")
    req = urllib.request.Request(URL, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=90) as response, partial.open("wb") as handle:
        total = int(response.headers.get("Content-Length", "0") or "0")
        written = 0
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)
            written += len(chunk)
            print(f"GPL570.annot.gz: {written}/{total or '?'} bytes", flush=True)
            if total and written >= total:
                break
    partial.replace(DESTINATION)
    print(f"DOWNLOADED {DESTINATION} {DESTINATION.stat().st_size} bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
