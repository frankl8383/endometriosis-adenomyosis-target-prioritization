#!/usr/bin/env python3
"""Phase 0 metadata/accessibility audit for public sources.

This script intentionally downloads only small metadata pages. It does not
download expression matrices, FASTQ files, or GWAS summary-statistic archives.
"""

from __future__ import annotations

import csv
import html
import re
import sys
import textwrap
import time
import urllib.error
import urllib.request
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INVENTORY = PROJECT_ROOT / "data" / "dataset_inventory.tsv"
RAW_DIR = PROJECT_ROOT / "data" / "raw_metadata"
OUT_TSV = PROJECT_ROOT / "results" / "phase0_source_audit.tsv"
OUT_MD = PROJECT_ROOT / "results" / "phase0_source_audit.md"

USER_AGENT = "phase0-public-source-audit/0.1 (metadata only)"


def fetch_text(url: str, timeout: int = 30) -> tuple[int | None, str, str]:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as handle:
            status = getattr(handle, "status", 200)
            final_url = handle.geturl()
            body = handle.read().decode("utf-8", errors="replace")
            return status, final_url, body
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        return exc.code, url, body
    except Exception as exc:  # noqa: BLE001 - preserve audit failure in output
        return None, url, f"ERROR: {type(exc).__name__}: {exc}"


def first_value(lines: list[str], key: str) -> str:
    prefix = f"!Series_{key} = "
    for line in lines:
        if line.startswith(prefix):
            return line[len(prefix) :].strip()
    return ""


def all_values(lines: list[str], key: str) -> list[str]:
    prefix = f"!Series_{key} = "
    values: list[str] = []
    for line in lines:
        if line.startswith(prefix):
            values.append(line[len(prefix) :].strip())
    return values


def audit_geo(accession: str) -> dict[str, str]:
    url = (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        f"?acc={accession}&targ=self&form=text&view=quick"
    )
    status, final_url, text = fetch_text(url)
    (RAW_DIR / f"{accession}.soft.txt").write_text(text, encoding="utf-8")

    lines = text.splitlines()
    sample_ids = all_values(lines, "sample_id")
    supp_files = all_values(lines, "supplementary_file")
    pubmed_ids = all_values(lines, "pubmed_id")
    series_type = all_values(lines, "type")

    title = first_value(lines, "title")
    summary = first_value(lines, "summary")
    design = first_value(lines, "overall_design")

    processed_hint_terms = (
        "matrix",
        "count",
        "counts",
        "h5",
        "h5ad",
        "rds",
        "csv",
        "txt",
        "tsv",
        "xls",
        "xlsx",
    )
    raw_only_terms = ("RAW.tar", ".sra", ".fastq", ".fq")
    lower_supp = " ".join(supp_files).lower()
    has_processed_hint = any(term in lower_supp for term in processed_hint_terms)
    has_raw_hint = any(term.lower() in lower_supp for term in raw_only_terms)

    if status == 200 and title:
        status_note = "ok"
    elif status == 200:
        status_note = "fetched_but_no_series_title"
    else:
        status_note = f"http_{status}" if status is not None else "fetch_failed"

    matrix_status = "unknown"
    if has_processed_hint:
        matrix_status = "processed_or_table_supplement_hint"
    elif has_raw_hint:
        matrix_status = "raw_supplement_hint_only"
    elif "array" in " ".join(series_type).lower():
        matrix_status = "geo_series_matrix_likely"

    return {
        "audit_kind": "geo",
        "http_status": str(status) if status is not None else "",
        "final_url": final_url,
        "access_status": status_note,
        "title": title,
        "series_status": first_value(lines, "status"),
        "last_update": first_value(lines, "last_update_date"),
        "data_type_detected": "; ".join(series_type),
        "sample_count_detected": str(len(sample_ids)),
        "pubmed_ids": ";".join(pubmed_ids),
        "supplementary_file_count": str(len(supp_files)),
        "matrix_status": matrix_status,
        "supplementary_files": "; ".join(supp_files[:8]),
        "summary_short": shorten(summary or design, 260),
    }


def strip_html(text: str) -> str:
    text = re.sub(r"<script.*?</script>", " ", text, flags=re.I | re.S)
    text = re.sub(r"<style.*?</style>", " ", text, flags=re.I | re.S)
    text = re.sub(r"<[^>]+>", " ", text)
    text = html.unescape(text)
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def find_title(text: str) -> str:
    match = re.search(r"<title[^>]*>(.*?)</title>", text, flags=re.I | re.S)
    if match:
        return shorten(strip_html(match.group(1)), 180)
    h1 = re.search(r"<h1[^>]*>(.*?)</h1>", text, flags=re.I | re.S)
    if h1:
        return shorten(strip_html(h1.group(1)), 180)
    return ""


def extract_data_availability(plain: str) -> str:
    idx = plain.lower().find("data availability")
    if idx == -1:
        return ""
    return shorten(plain[idx : idx + 900], 500)


def audit_web(resource: str, url: str) -> dict[str, str]:
    status, final_url, text = fetch_text(url)
    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", resource)[:80]
    (RAW_DIR / f"{safe_name}.html.txt").write_text(text, encoding="utf-8")
    plain = strip_html(text)
    zenodo = ";".join(sorted(set(re.findall(r"https?://zenodo\.org/records/\d+", text))))
    data_availability = extract_data_availability(plain)

    return {
        "audit_kind": "web",
        "http_status": str(status) if status is not None else "",
        "final_url": final_url,
        "access_status": "ok" if status and 200 <= status < 400 else "check_manually",
        "title": find_title(text),
        "series_status": "",
        "last_update": "",
        "data_type_detected": "",
        "sample_count_detected": "",
        "pubmed_ids": "",
        "supplementary_file_count": "",
        "matrix_status": "web_source_check",
        "supplementary_files": zenodo,
        "summary_short": data_availability or shorten(plain, 260),
    }


def shorten(text: str, width: int) -> str:
    text = re.sub(r"\s+", " ", text or "").strip()
    if len(text) <= width:
        return text
    return textwrap.shorten(text, width=width, placeholder="...")


def classify_decision(row: dict[str, str], audit: dict[str, str]) -> tuple[str, str]:
    role = row["role"]
    resource = row["accession_or_resource"]
    status = audit["access_status"]
    sample_count = int(audit["sample_count_detected"] or 0)
    matrix_status = audit["matrix_status"]
    summary = audit["summary_short"].lower()
    supp_files = audit["supplementary_files"].lower()

    if status != "ok":
        return "manual_check", "Metadata fetch failed or source requires manual inspection."

    if role == "genetics":
        if "zenodo" in supp_files or "summary" in summary or "data availability" in summary:
            return "primary_if_downloadable", "Genetics source is visible; verify downloadable summary-stat files next."
        return "manual_check", "Need to find exact summary-statistic download link."

    if resource.startswith("GSE"):
        if sample_count >= 20 and matrix_status != "raw_supplement_hint_only":
            return "likely_pass", "GEO metadata looks usable; inspect series matrix or supplementary file format next."
        if sample_count >= 20 and matrix_status == "raw_supplement_hint_only":
            return "conditional", "Sample count is strong, but current metadata suggests raw archive; inspect tar contents before committing."
        if sample_count > 0:
            return "conditional", "Metadata exists but sample count or matrix hint is limited."
        return "manual_check", "GEO accession fetched but samples were not detected."

    if "data availability" in summary or "gsa-human" in summary or "zenodo" in summary:
        return "conditional", "Article/source has a data availability signal; exact matrix accessibility needs manual confirmation."

    return "manual_check", "Needs manual inspection for actual data access."


def main() -> int:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)

    with INVENTORY.open("r", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    out_rows: list[dict[str, str]] = []
    for row in rows:
        resource = row["accession_or_resource"]
        url = row["source_url"]
        print(f"Auditing {resource}...", file=sys.stderr)
        if resource.startswith("GSE"):
            audit = audit_geo(resource)
        else:
            audit = audit_web(resource, url)

        decision, rationale = classify_decision(row, audit)
        out_rows.append(
            {
                **row,
                **audit,
                "phase0_decision": decision,
                "decision_rationale": rationale,
            }
        )
        time.sleep(0.34)

    fieldnames = list(out_rows[0].keys())
    with OUT_TSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    write_markdown(out_rows)
    print(f"Wrote {OUT_TSV}", file=sys.stderr)
    print(f"Wrote {OUT_MD}", file=sys.stderr)
    return 0


def write_markdown(rows: list[dict[str, str]]) -> None:
    lines = [
        "# Phase 0 Source Audit",
        "",
        "This audit checks lightweight public metadata only. It does not download large matrices, FASTQ files, or GWAS archives.",
        "",
        "## Summary Table",
        "",
        "| Resource | Role | Decision | Samples | Matrix/access hint | Rationale |",
        "|---|---|---|---:|---|---|",
    ]
    for row in rows:
        lines.append(
            "| {resource} | {role} | {decision} | {samples} | {matrix} | {rationale} |".format(
                resource=row["accession_or_resource"],
                role=row["role"],
                decision=row["phase0_decision"],
                samples=row["sample_count_detected"] or "",
                matrix=row["matrix_status"],
                rationale=row["decision_rationale"],
            )
        )

    lines.extend(["", "## Notes", ""])
    for row in rows:
        lines.extend(
            [
                f"### {row['accession_or_resource']}",
                "",
                f"- Title: {row['title'] or 'not detected'}",
                f"- URL: {row['final_url']}",
                f"- Summary: {row['summary_short'] or 'not detected'}",
                f"- Supplementary/link hints: {row['supplementary_files'] or 'not detected'}",
                "",
            ]
        )

    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())

