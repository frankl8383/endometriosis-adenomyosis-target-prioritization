#!/usr/bin/env python3
"""Audit GSE179640 tar contents, sample metadata, and 10x h5 readability."""

from __future__ import annotations

import csv
import json
import re
import shutil
import sys
import tarfile
import tempfile
from collections import Counter
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW = PROJECT_ROOT / "data" / "raw_downloads"
META = PROJECT_ROOT / "data" / "raw_metadata"
RESULTS = PROJECT_ROOT / "results" / "singlecell"
PY_PKG_DIR = PROJECT_ROOT / "tools" / "py_packages"

if sys.version_info[:2] != (3, 12):
    raise SystemExit(
        "This script needs the project Python 3.12 runtime because local h5py/anndata "
        "wheels were built for CPython 3.12. Use: "
        "python3 "
        "scripts/audit_gse179640_tar_and_metadata.py"
    )

sys.path.insert(0, str(PY_PKG_DIR))

import h5py  # noqa: E402


def decode_array(values) -> list[str]:
    out: list[str] = []
    for value in values:
        if isinstance(value, bytes):
            out.append(value.decode("utf-8", errors="replace"))
        else:
            out.append(str(value))
    return out


def parse_soft_samples(path: Path) -> list[dict[str, object]]:
    samples: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    characteristics: dict[str, str] = {}
    supplementary: list[str] = []

    def finish() -> None:
        if current is None:
            return
        current.update(characteristics)
        current["supplementary_files"] = "|".join(supplementary)
        current["supplementary_basenames"] = "|".join(Path(item).name for item in supplementary)
        samples.append(current.copy())

    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if line.startswith("^SAMPLE = "):
            finish()
            accession = line.split("=", 1)[1].strip()
            current = {"geo_accession": accession}
            characteristics = {}
            supplementary = []
            continue
        if current is None:
            continue
        if line.startswith("!Sample_title = "):
            current["title"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Sample_geo_accession = "):
            current["geo_accession"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Sample_source_name_ch1 = "):
            current["source_name"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Sample_description = "):
            current["description"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Sample_characteristics_ch1 = "):
            payload = line.split("=", 1)[1].strip()
            if ":" in payload:
                key, value = payload.split(":", 1)
                characteristics[key.strip().replace(" ", "_")] = value.strip()
        elif line.startswith("!Sample_supplementary_file"):
            supplementary.append(line.split("=", 1)[1].strip())
        elif line.startswith("!Sample_library_strategy = "):
            current["library_strategy"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Sample_library_source = "):
            current["library_source"] = line.split("=", 1)[1].strip()
    finish()
    return samples


def load_candidates(path: Path) -> tuple[list[dict[str, str]], set[str], set[str]]:
    records: list[dict[str, str]] = []
    symbols: set[str] = set()
    ensembl_ids: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            symbol = row.get("gene_symbol", "")
            ensembl_id = row.get("gene_id", "").split(".")[0]
            records.append({"gene_symbol": symbol, "gene_id": ensembl_id})
            if symbol:
                symbols.add(symbol)
            if ensembl_id:
                ensembl_ids.add(ensembl_id)
    return records, symbols, ensembl_ids


def classify_entry(name: str) -> dict[str, str]:
    accession = ""
    sample_token = ""
    match = re.match(r"^(GSM\d+)_([^/]+)$", name)
    if match:
        accession = match.group(1)
        sample_token = match.group(2)

    if name.endswith("_filtered_feature_bc_matrix.h5"):
        kind = "filtered_feature_bc_matrix_h5"
    elif name.endswith("_bulk.featurecounts.txt.gz"):
        kind = "bulk_featurecounts"
    elif name.endswith("_barcodes.tsv.gz") or name.endswith("_features.tsv.gz") or name.endswith("_matrix.mtx.gz"):
        kind = "cell_hashing_mtx_component"
    else:
        kind = "other"

    code = ""
    tissue_code = ""
    token_match = re.search(r"_(C\d+|E\d+|EOR\d+)_?([A-Za-z0-9]*)", name)
    if token_match:
        code = token_match.group(1)
        tissue_code = token_match.group(2)
    return {
        "tar_entry": name,
        "geo_accession": accession,
        "sample_token": sample_token,
        "entry_kind": kind,
        "inferred_subject_code": code,
        "inferred_tissue_code": tissue_code,
    }


def inspect_10x_h5_from_tar(
    tar: tarfile.TarFile,
    member: tarfile.TarInfo,
    candidate_records: list[dict[str, str]],
    candidate_symbols: set[str],
    candidate_ensembl_ids: set[str],
) -> dict[str, object]:
    extracted = tar.extractfile(member)
    if extracted is None:
        raise RuntimeError(f"Could not extract {member.name}")
    with extracted, tempfile.NamedTemporaryFile(suffix=".h5") as tmp:
        shutil.copyfileobj(extracted, tmp)
        tmp.flush()
        with h5py.File(tmp.name, "r") as h5:
            if "matrix" not in h5:
                raise RuntimeError("missing /matrix group")
            matrix = h5["matrix"]
            shape = [int(x) for x in matrix["shape"][()]]
            genes = int(shape[0])
            cells = int(shape[1])
            nnz = int(matrix["data"].shape[0])
            barcodes = int(matrix["barcodes"].shape[0])
            features = matrix["features"]
            feature_ids = decode_array(features["id"][:])
            feature_names = decode_array(features["name"][:])
            feature_types = decode_array(features["feature_type"][:]) if "feature_type" in features else []
            genomes = decode_array(features["genome"][:]) if "genome" in features else []

    feature_id_set = {item.split(".")[0] for item in feature_ids}
    feature_name_set = set(feature_names)
    present_symbol = sorted(candidate_symbols & feature_name_set)
    present_ensembl = sorted(candidate_ensembl_ids & feature_id_set)
    present_records = [
        row
        for row in candidate_records
        if row["gene_symbol"] in feature_name_set or row["gene_id"] in feature_id_set
    ]
    return {
        "tar_entry": member.name,
        "h5_status": "readable",
        "n_features": genes,
        "n_barcodes": cells,
        "n_nonzero": nnz,
        "barcodes_dataset_n": barcodes,
        "feature_type_counts": json.dumps(Counter(feature_types), sort_keys=True),
        "genome_counts": json.dumps(Counter(genomes), sort_keys=True),
        "candidate_symbol_present_n": len(present_symbol),
        "candidate_ensembl_present_n": len(present_ensembl),
        "candidate_record_present_n": len(present_records),
        "candidate_symbols_present": "|".join(present_symbol),
        "candidate_ensembl_present": "|".join(present_ensembl),
        "candidate_records_present": "|".join(
            f"{row['gene_symbol']}:{row['gene_id']}" for row in present_records
        ),
    }


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    tar_path = RAW / "GSE179640__GSE179640_RAW.tar"
    soft_path = META / "GSE179640_family.soft.txt"
    candidate_path = PROJECT_ROOT / "results" / "gwas" / "gwas_candidate_gene_universe.tsv"
    if not tar_path.exists():
        raise SystemExit(f"Missing {tar_path}")
    if not soft_path.exists():
        raise SystemExit(f"Missing {soft_path}")

    samples = parse_soft_samples(soft_path)
    sample_by_accession = {str(row["geo_accession"]): row for row in samples}
    candidate_records, candidate_symbols, candidate_ensembl_ids = load_candidates(candidate_path)

    manifest_rows: list[dict[str, object]] = []
    h5_rows: list[dict[str, object]] = []
    with tarfile.open(tar_path, "r") as tar:
        members = [member for member in tar.getmembers() if member.isfile()]
        for member in members:
            row: dict[str, object] = classify_entry(member.name)
            row["size_bytes"] = member.size
            row["tar_mtime"] = member.mtime
            meta = sample_by_accession.get(str(row["geo_accession"]), {})
            for field in [
                "title",
                "source_name",
                "condition",
                "sample_location",
                "tissue",
                "method",
                "library_type",
                "description",
                "library_strategy",
                "library_source",
            ]:
                row[field] = meta.get(field, "")
            row["metadata_joined"] = bool(meta)
            manifest_rows.append(row)

            if row["entry_kind"] == "filtered_feature_bc_matrix_h5":
                try:
                    h5_row = inspect_10x_h5_from_tar(
                        tar,
                        member,
                        candidate_records,
                        candidate_symbols,
                        candidate_ensembl_ids,
                    )
                except Exception as exc:  # noqa: BLE001
                    h5_row = {"tar_entry": member.name, "h5_status": f"failed:{type(exc).__name__}:{exc}"}
                h5_row.update(row)
                h5_rows.append(h5_row)

    sc_meta = [
        row
        for row in samples
        if row.get("method") == "scRNA-seq" and row.get("library_type") == "Gene Expression"
    ]
    h5_accessions = {str(row["geo_accession"]) for row in h5_rows}
    sc_meta_accessions = {str(row["geo_accession"]) for row in sc_meta}
    missing_from_tar = sorted(sc_meta_accessions - h5_accessions)
    h5_not_in_sc_meta = sorted(h5_accessions - sc_meta_accessions)

    write_tsv(RESULTS / "GSE179640_family_sample_metadata.tsv", samples)
    write_tsv(RESULTS / "GSE179640_tar_manifest.tsv", manifest_rows)
    write_tsv(RESULTS / "GSE179640_10x_h5_audit.tsv", h5_rows)

    entry_counts = Counter(str(row["entry_kind"]) for row in manifest_rows)
    sample_location_counts = Counter(str(row.get("sample_location", "")) for row in h5_rows)
    tissue_counts = Counter(str(row.get("tissue", "")) for row in h5_rows)
    condition_counts = Counter(str(row.get("condition", "")) for row in h5_rows)
    readable_h5 = [row for row in h5_rows if row.get("h5_status") == "readable"]
    total_cells = sum(int(row.get("n_barcodes", 0) or 0) for row in readable_h5)
    total_nnz = sum(int(row.get("n_nonzero", 0) or 0) for row in readable_h5)
    candidate_min = min((int(row.get("candidate_record_present_n", 0) or 0) for row in readable_h5), default=0)
    candidate_max = max((int(row.get("candidate_record_present_n", 0) or 0) for row in readable_h5), default=0)

    verdict = "PASS_WITH_CONDITIONS"
    blockers: list[str] = []
    if missing_from_tar:
        blockers.append(f"Sample metadata accessions missing from tar h5: {', '.join(missing_from_tar)}")
    if h5_not_in_sc_meta:
        blockers.append(f"Tar h5 accessions not classified as scRNA-seq Gene Expression in metadata: {', '.join(h5_not_in_sc_meta)}")
    if any(row.get("h5_status") != "readable" for row in h5_rows):
        blockers.append("At least one 10x h5 file failed readability checks.")
    if blockers:
        verdict = "REVIEW_REQUIRED"

    lines = [
        "# GSE179640 tar and metadata audit",
        "",
        f"- Tar file: `{tar_path}`",
        f"- Family SOFT: `{soft_path}`",
        f"- Parsed GEO samples: {len(samples)}",
        f"- Tar file entries: {len(manifest_rows)}",
        f"- 10x h5 gene-expression entries inspected: {len(h5_rows)}",
        f"- Readable 10x h5 entries: {len(readable_h5)}",
        f"- Total barcodes across readable h5 entries: {total_cells}",
        f"- Total non-zero counts across readable h5 entries: {total_nnz}",
        f"- Candidate gene-record coverage per h5, symbol or Ensembl match: min {candidate_min}, max {candidate_max} of {len(candidate_records)} GWAS candidate records",
        f"- Entry type counts: `{json.dumps(entry_counts, sort_keys=True)}`",
        f"- h5 condition counts: `{json.dumps(condition_counts, sort_keys=True)}`",
        f"- h5 sample location counts: `{json.dumps(sample_location_counts, sort_keys=True)}`",
        f"- h5 tissue counts: `{json.dumps(tissue_counts, sort_keys=True)}`",
        f"- Metadata scRNA-seq Gene Expression samples missing from tar h5: `{','.join(missing_from_tar) if missing_from_tar else 'none'}`",
        f"- Tar h5 accessions not matching scRNA-seq Gene Expression metadata: `{','.join(h5_not_in_sc_meta) if h5_not_in_sc_meta else 'none'}`",
        "",
        "## Self-review",
        "",
        f"Verdict: **{verdict}**",
        "",
        "Pass criteria checked:",
        "",
        "- The completed tar is readable by `tarfile` and contains expected GEO supplementary files.",
        "- Each `filtered_feature_bc_matrix.h5` is opened with `h5py` and has a 10x `/matrix` group.",
        "- GEO sample metadata are parsed from official family SOFT and joined by GSM accession rather than inferred only from filenames.",
        "- Candidate gene coverage is measured against both HGNC symbols and Ensembl IDs.",
        "",
        "Conditions before biological interpretation:",
        "",
        "- These are raw/filtered 10x matrices without author cell-type annotations; cell-state claims require downstream QC, clustering/annotation, or independently sourced author labels.",
        "- Organoid and cell-hashing entries should be excluded from primary lesion cell-state localization unless analyzed in a dedicated organoid subanalysis.",
        "- Ectopic peritoneum, ectopic adjacent peritoneum, ectopic ovary, eutopic endometrium, and control endometrium must remain separate until a donor-aware model is specified.",
        "",
    ]
    if blockers:
        lines.extend(["Review issues:", ""])
        lines.extend(f"- {item}" for item in blockers)
        lines.append("")
    (RESULTS / "GSE179640_tar_metadata_audit.md").write_text("\n".join(lines), encoding="utf-8")
    print(RESULTS / "GSE179640_tar_metadata_audit.md")
    print(f"VERDICT {verdict}")
    return 0 if verdict == "PASS_WITH_CONDITIONS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
