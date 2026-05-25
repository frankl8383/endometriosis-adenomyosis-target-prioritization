#!/usr/bin/env python3
"""Audit phase 1 downloads before biological analysis.

The audit is intentionally conservative. A file is not considered analysis-ready
unless it exists, has a plausible size relative to the manifest, and passes a
lightweight format check when the extension makes that possible.
"""

from __future__ import annotations

import csv
import gzip
import zipfile
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MANIFEST = PROJECT_ROOT / "data" / "download_manifest.tsv"
RAW_DOWNLOADS = PROJECT_ROOT / "data" / "raw_downloads"
OUT_TSV = PROJECT_ROOT / "results" / "phase1_download_audit.tsv"
OUT_MD = PROJECT_ROOT / "results" / "phase1_download_audit.md"

CORE_PHASES = {"phase1_small", "phase1_gwas"}


def safe_filename(resource: str, file_or_record: str) -> str:
    resource_safe = resource.replace(" ", "_").replace("/", "_")
    return f"{resource_safe}__{file_or_record}"


def expected_size_bytes(size_mb: str) -> int | None:
    try:
        size = float(size_mb)
    except ValueError:
        return None
    if size <= 0:
        return None
    return int(size * 1024 * 1024)


def size_status(observed: int, expected: int | None) -> str:
    if observed <= 0:
        return "empty"
    if expected is None:
        return "unknown_expected_size"
    if expected < 1024 * 1024:
        return "small_file_present"
    ratio = observed / expected
    if 0.97 <= ratio <= 1.03:
        return "plausible"
    if ratio < 0.97:
        return "too_small"
    return "too_large"


def inspect_gzip(path: Path) -> tuple[str, str]:
    try:
        preview: list[str] = []
        with gzip.open(path, "rb") as handle:
            raw_preview = handle.read(4096)
            for line in raw_preview.decode("utf-8", errors="replace").splitlines()[:3]:
                preview.append(line.strip())
            while handle.read(1024 * 1024):
                pass
        return "ok", " | ".join(preview)[:500]
    except Exception as exc:  # noqa: BLE001 - audit should not stop on one file
        return "fail", f"{type(exc).__name__}: {exc}"


def inspect_plain_text(path: Path) -> tuple[str, str]:
    try:
        with path.open("rt", encoding="utf-8", errors="replace") as handle:
            preview = [handle.readline().strip() for _ in range(3)]
        return "ok", " | ".join([line for line in preview if line])[:500]
    except Exception as exc:  # noqa: BLE001
        return "fail", f"{type(exc).__name__}: {exc}"


def inspect_xlsx(path: Path) -> tuple[str, str]:
    try:
        if not zipfile.is_zipfile(path):
            return "fail", "not_a_zip_container"
        with zipfile.ZipFile(path) as workbook:
            bad_member = workbook.testzip()
            if bad_member:
                return "fail", f"corrupt_zip_member:{bad_member}"
            names = workbook.namelist()
        sheet_count = sum(name.startswith("xl/worksheets/") for name in names)
        return "ok", f"zip_container_ok; worksheet_files={sheet_count}"
    except Exception as exc:  # noqa: BLE001
        return "fail", f"{type(exc).__name__}: {exc}"


def inspect_format(path: Path) -> tuple[str, str]:
    name = path.name.lower()
    if name.endswith(".gz"):
        return inspect_gzip(path)
    if name.endswith(".xlsx"):
        return inspect_xlsx(path)
    if name.endswith(".txt") or name.endswith(".tsv") or name.endswith(".csv"):
        return inspect_plain_text(path)
    return "not_checked", "no_lightweight_format_check_for_extension"


def readiness(file_state: str, size_state: str, format_state: str) -> str:
    if file_state != "complete_file_present":
        return "not_ready"
    if size_state not in {"plausible", "unknown_expected_size", "small_file_present"}:
        return "not_ready"
    if format_state == "fail":
        return "not_ready"
    return "ready"


def audit_row(row: dict[str, str]) -> dict[str, str]:
    destination = RAW_DOWNLOADS / safe_filename(row["resource"], row["file_or_record"])
    partial = destination.with_suffix(destination.suffix + ".partial")
    expected = expected_size_bytes(row["size_mb"])

    if destination.exists():
        file_state = "complete_file_present"
        observed = destination.stat().st_size
        format_state, preview = inspect_format(destination)
    elif partial.exists():
        file_state = "partial_file_only"
        observed = partial.stat().st_size
        format_state, preview = "not_checked", "partial download; do not analyze"
    else:
        file_state = "missing"
        observed = 0
        format_state, preview = "not_checked", "missing"

    size_state = size_status(observed, expected)
    return {
        "priority": row["priority"],
        "resource": row["resource"],
        "file": row["file_or_record"],
        "download_phase": row["download_phase"],
        "purpose": row["purpose"],
        "expected_bytes": str(expected or ""),
        "observed_bytes": str(observed),
        "file_state": file_state,
        "size_status": size_state,
        "format_status": format_state,
        "analysis_readiness": readiness(file_state, size_state, format_state),
        "preview_or_error": preview,
    }


def write_markdown(rows: list[dict[str, str]]) -> None:
    ready = sum(row["analysis_readiness"] == "ready" for row in rows)
    core = [row for row in rows if row["download_phase"] in CORE_PHASES]
    core_ready = sum(row["analysis_readiness"] == "ready" for row in core)

    lines = [
        "# Phase 1 Download Audit",
        "",
        "## Summary",
        "",
        f"- Ready files: {ready}/{len(rows)}",
        f"- Core phase 1 ready files: {core_ready}/{len(core)}",
        "",
        "## Interpretation",
        "",
        "Only rows marked `ready` can enter analysis. `partial_file_only`, size-mismatched files, and failed format checks must be downloaded again or resumed before modeling.",
        "",
        "## Core Files",
        "",
        "| Resource | File | Phase | State | Size | Format | Ready |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in rows:
        if row["download_phase"] not in CORE_PHASES:
            continue
        lines.append(
            "| {resource} | `{file}` | {download_phase} | {file_state} | {size_status} | {format_status} | {analysis_readiness} |".format(
                **row
            )
        )

    lines.extend(
        [
            "",
            "## Required Review Gate",
            "",
            "Phase 1 cannot pass until the three core GWAS files and the four core bulk/immune validation files are `ready`. The adenomyosis h5ad remains a conditional large-file gate and should be audited separately before single-cell/spatial modeling.",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    with MANIFEST.open("r", encoding="utf-8") as handle:
        rows = [audit_row(row) for row in csv.DictReader(handle, delimiter="\t")]

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
