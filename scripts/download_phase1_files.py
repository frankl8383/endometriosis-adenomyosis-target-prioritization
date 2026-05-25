#!/usr/bin/env python3
"""Download the phase 1 public-data set.

The manifest keeps the download list reproducible. This script skips files that
already exist with nonzero size. By default it downloads the core phase 1 rows
(`phase1_small` and `phase1_gwas`). It writes incomplete downloads to
`.partial` files and resumes those files when the server supports byte ranges.
"""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MANIFEST = PROJECT_ROOT / "data" / "download_manifest.tsv"
RAW_DOWNLOADS = PROJECT_ROOT / "data" / "raw_downloads"
LOG_FILE = PROJECT_ROOT / "logs" / "phase1_downloads.tsv"

USER_AGENT = "phase1-public-data-downloader/0.1"
PHASE_GROUPS = {
    "phase1": {"phase1_small", "phase1_gwas"},
    "phase1_core": {"phase1_small", "phase1_gwas"},
    "phase1_all": {"phase1_small", "phase1_gwas", "phase1_optional"},
}


def safe_filename(resource: str, file_or_record: str) -> str:
    resource_safe = resource.replace(" ", "_").replace("/", "_")
    return f"{resource_safe}__{file_or_record}"


def expected_size_bytes(row: dict[str, str]) -> int | None:
    try:
        size_mb = float(row["size_mb"])
    except (KeyError, TypeError, ValueError):
        return None
    if size_mb <= 0:
        return None
    return int(size_mb * 1024 * 1024)


def file_is_plausibly_complete(destination: Path, expected_bytes: int | None) -> bool:
    if not destination.exists() or destination.stat().st_size == 0:
        return False
    if expected_bytes is None:
        return True
    if expected_bytes < 1024 * 1024:
        return True
    observed = destination.stat().st_size
    lower = int(expected_bytes * 0.97)
    upper = int(expected_bytes * 1.03)
    return lower <= observed <= upper


def download_urllib(url: str, partial_destination: Path) -> int:
    headers = {"User-Agent": USER_AGENT}
    resume_from = partial_destination.stat().st_size if partial_destination.exists() else 0
    mode = "ab" if resume_from > 0 else "wb"
    if resume_from > 0:
        headers["Range"] = f"bytes={resume_from}-"

    req = urllib.request.Request(url, headers=headers)
    try:
        response = urllib.request.urlopen(req, timeout=60)
    except urllib.error.HTTPError as exc:
        if exc.code == 416 and resume_from > 0:
            return resume_from
        raise

    with response, partial_destination.open(mode) as out:
        status = getattr(response, "status", 200)
        if resume_from > 0 and status == 200:
            # Server ignored Range. Restart to avoid duplicated content.
            out.close()
            partial_destination.unlink()
            return download_urllib(url, partial_destination)

        total = 0
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
            total += len(chunk)
    return resume_from + total


def download_curl(url: str, partial_destination: Path) -> int:
    curl = shutil.which("curl")
    if curl is None:
        raise RuntimeError("curl was requested but was not found on PATH")
    partial_destination.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        curl,
        "--location",
        "--fail",
        "--continue-at",
        "-",
        "--retry",
        "8",
        "--retry-delay",
        "5",
        "--connect-timeout",
        "30",
        "--speed-limit",
        "1024",
        "--speed-time",
        "90",
        "--user-agent",
        USER_AGENT,
        "--output",
        str(partial_destination),
        url,
    ]
    subprocess.run(cmd, check=True)
    return partial_destination.stat().st_size


def download_aria2(url: str, partial_destination: Path, relaxed: bool = False) -> int:
    aria2 = shutil.which("aria2c")
    if aria2 is None:
        raise RuntimeError("aria2c was requested but was not found on PATH")
    partial_destination.parent.mkdir(parents=True, exist_ok=True)
    max_tries = "0" if relaxed else "8"
    retry_wait = "15" if relaxed else "5"
    cmd = [
        aria2,
        "--continue=true",
        "--max-connection-per-server=8",
        "--split=8",
        "--min-split-size=1M",
        f"--max-tries={max_tries}",
        f"--retry-wait={retry_wait}",
        "--connect-timeout=30",
        "--allow-overwrite=true",
        "--auto-file-renaming=false",
        "--summary-interval=30",
        "--user-agent",
        USER_AGENT,
        "--dir",
        str(partial_destination.parent),
        "--out",
        partial_destination.name,
        url,
    ]
    if not relaxed:
        cmd.insert(8, "--lowest-speed-limit=1K")
    subprocess.run(cmd, check=True)
    return partial_destination.stat().st_size


def download(url: str, partial_destination: Path, backend: str) -> int:
    if backend == "urllib":
        return download_urllib(url, partial_destination)
    if backend == "curl":
        return download_curl(url, partial_destination)
    if backend == "aria2":
        return download_aria2(url, partial_destination)
    if backend == "aria2_relaxed":
        return download_aria2(url, partial_destination, relaxed=True)
    raise ValueError(f"Unsupported backend: {backend}")


def select_rows(rows: list[dict[str, str]], requested_phase: str) -> list[dict[str, str]]:
    phases = PHASE_GROUPS.get(requested_phase, {requested_phase})
    return [row for row in rows if row["download_phase"] in phases]


def main() -> int:
    phase = sys.argv[1] if len(sys.argv) > 1 else "phase1"
    backend = sys.argv[2] if len(sys.argv) > 2 else "urllib"
    RAW_DOWNLOADS.mkdir(parents=True, exist_ok=True)
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

    with MANIFEST.open("r", encoding="utf-8") as handle:
        rows = select_rows(list(csv.DictReader(handle, delimiter="\t")), phase)

    if not rows:
        print(f"No manifest rows for phase: {phase}", file=sys.stderr)
        return 1

    log_exists = LOG_FILE.exists()
    with LOG_FILE.open("a", encoding="utf-8", newline="") as log_handle:
        fieldnames = ["timestamp", "phase", "resource", "file", "destination", "status", "bytes"]
        writer = csv.DictWriter(log_handle, fieldnames=fieldnames, delimiter="\t")
        if not log_exists:
            writer.writeheader()

        for row in rows:
            destination = RAW_DOWNLOADS / safe_filename(row["resource"], row["file_or_record"])
            partial_destination = destination.with_suffix(destination.suffix + ".partial")
            expected_bytes = expected_size_bytes(row)
            timestamp = time.strftime("%Y-%m-%d %H:%M:%S")

            if file_is_plausibly_complete(destination, expected_bytes):
                status = "skipped_exists"
                byte_count = destination.stat().st_size
                print(f"SKIP {destination.name} ({byte_count} bytes)")
            else:
                print(f"DOWNLOAD {destination.name}")
                try:
                    byte_count = download(row["url"], partial_destination, backend)
                    partial_destination.replace(destination)
                    status = f"downloaded_needs_audit:{backend}"
                except Exception as exc:  # noqa: BLE001 - log and continue
                    status = f"failed:{type(exc).__name__}:{exc}"
                    byte_count = partial_destination.stat().st_size if partial_destination.exists() else 0
                    print(f"FAILED {destination.name}: {exc}", file=sys.stderr)

            writer.writerow(
                {
                    "timestamp": timestamp,
                    "phase": phase,
                    "resource": row["resource"],
                    "file": row["file_or_record"],
                    "destination": str(destination),
                    "status": status,
                    "bytes": byte_count,
                }
            )
            log_handle.flush()

    print(f"Log written to {LOG_FILE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
