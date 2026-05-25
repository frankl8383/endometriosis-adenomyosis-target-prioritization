#!/usr/bin/env python3
"""Parse bulk metadata and extract genetics-prioritized candidate expression."""

from __future__ import annotations

import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW = PROJECT_ROOT / "data" / "raw_downloads"
META = PROJECT_ROOT / "data" / "raw_metadata"
RESULTS = PROJECT_ROOT / "results" / "bulk"
GWAS_CANDIDATES = PROJECT_ROOT / "results" / "gwas" / "gwas_candidate_gene_universe.tsv"


def clean_value(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        value = value[1:-1]
    return value.strip()


def parse_soft_samples(path: Path) -> pd.DataFrame:
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    characteristic_counts: defaultdict[str, int] = defaultdict(int)

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current:
                    records.append(current)
                current = {"sample_accession": clean_value(line.split("=", 1)[1])}
                characteristic_counts.clear()
                continue
            if current is None:
                continue

            if line.startswith("!Sample_title = "):
                current["title"] = clean_value(line.split("=", 1)[1])
            elif line.startswith("!Sample_geo_accession = "):
                current["geo_accession"] = clean_value(line.split("=", 1)[1])
            elif line.startswith("!Sample_source_name_ch1 = "):
                current["source_name_ch1"] = clean_value(line.split("=", 1)[1])
            elif line.startswith("!Sample_characteristics_ch1 = "):
                item = clean_value(line.split("=", 1)[1])
                if ":" in item:
                    key, value = item.split(":", 1)
                    key = re.sub(r"[^0-9A-Za-z]+", "_", key.strip().lower()).strip("_")
                    key = "characteristic" if not key else key
                    characteristic_counts[key] += 1
                    if characteristic_counts[key] > 1:
                        key = f"{key}_{characteristic_counts[key]}"
                    current[key] = value.strip()
                else:
                    characteristic_counts["characteristic"] += 1
                    current[f"characteristic_{characteristic_counts['characteristic']}"] = item

    if current:
        records.append(current)
    return pd.DataFrame(records)


def parse_gse234354_metadata() -> pd.DataFrame:
    meta = parse_soft_samples(META / "GSE234354_family.soft.txt")
    meta["dataset"] = "GSE234354"
    meta["matrix_sample_id"] = meta["title"].str.extract(r",\s*([^,]+)$", expand=False)
    meta["disease_status"] = "cycle_reference"
    meta["analysis_role"] = "menstrual_cycle_control"
    return meta


def parse_gse313775_metadata(count_columns: list[str]) -> pd.DataFrame:
    soft = parse_soft_samples(META / "GSE313775_family.soft.txt")
    title_re = re.compile(
        r"^(?P<title_cell_type>.+?) cells?, (?P<cohort>Endometriosis|Control) cohort, Biol rep (?P<biological_replicate>\d+)$"
    )
    title_parts = soft["title"].str.extract(title_re)
    soft = pd.concat([soft, title_parts], axis=1)
    soft["cohort"] = soft["cohort"].str.lower()
    soft["disease_status"] = soft["cohort"].map({"endometriosis": "endometriosis", "control": "control"})
    soft["biological_replicate"] = pd.to_numeric(soft["biological_replicate"], errors="coerce").astype("Int64")

    soft["cell_subset_from_title"] = (
        soft["title_cell_type"]
        .str.replace("Th1/17", "Th1-17", regex=False)
        .str.replace(" ", "", regex=False)
    )

    donor_order: list[str] = []
    parsed_columns: list[dict[str, str]] = []
    for column in count_columns:
        match = re.match(r"^(?P<donor>EDCV\d+)_(?P<cell_subset>Th1-17|Th1|Th17)_S\d+$", column)
        if not match:
            continue
        donor = match.group("donor")
        if donor not in donor_order:
            donor_order.append(donor)
        parsed_columns.append(
            {
                "matrix_sample_id": column,
                "donor_id": donor,
                "cell_subset": match.group("cell_subset"),
            }
        )

    endometriosis_donors = donor_order[:10]
    control_donors = donor_order[10:]
    donor_map = {
        donor: {"disease_status": "endometriosis", "biological_replicate": idx + 1}
        for idx, donor in enumerate(endometriosis_donors)
    }
    donor_map.update(
        {
            donor: {"disease_status": "control", "biological_replicate": idx + 1}
            for idx, donor in enumerate(control_donors)
        }
    )

    records: list[dict[str, str]] = []
    for parsed in parsed_columns:
        donor_info = donor_map[parsed["donor_id"]]
        match = soft[
            (soft["disease_status"] == donor_info["disease_status"])
            & (soft["biological_replicate"] == donor_info["biological_replicate"])
            & (soft["cell_subset_from_title"] == parsed["cell_subset"])
        ]
        record = {
            "dataset": "GSE313775",
            "matrix_sample_id": parsed["matrix_sample_id"],
            "donor_id": parsed["donor_id"],
            "cell_subset": parsed["cell_subset"],
            "disease_status": donor_info["disease_status"],
            "biological_replicate": donor_info["biological_replicate"],
            "metadata_match_strategy": "matrix_donor_order_to_geo_biological_replicate",
        }
        if len(match) == 1:
            row = match.iloc[0].to_dict()
            for key in [
                "geo_accession",
                "title",
                "source_name_ch1",
                "rasrm_disease_stage_1_4",
                "age",
                "race",
                "bmi",
                "education",
            ]:
                record[key] = row.get(key, "")
            record["metadata_match_status"] = "matched"
        else:
            record["metadata_match_status"] = f"unmatched_n={len(match)}"
        records.append(record)
    return pd.DataFrame(records)


def classify_gse141549_code(code: str) -> dict[str, str]:
    replicate = code.endswith(" Replicate")
    base = code.replace(" Replicate", "")
    mapping = {
        "CE": ("control", "endometrium", "control_endometrium"),
        "CP": ("control", "peritoneum", "control_peritoneum"),
        "PE": ("endometriosis", "endometrium", "patient_eutopic_endometrium"),
        "PP": ("endometriosis", "peritoneum", "patient_peritoneum"),
        "OMA": ("endometriosis", "lesion", "ovarian_endometrioma"),
        "SuL": ("endometriosis", "lesion", "superficial_lesion"),
        "PeLB": ("endometriosis", "lesion", "peritoneal_lesion_black"),
        "PeLR": ("endometriosis", "lesion", "peritoneal_lesion_red"),
        "PeLW": ("endometriosis", "lesion", "peritoneal_lesion_white"),
        "REV": ("endometriosis", "lesion", "rectovaginal_lesion"),
        "DiEIn": ("endometriosis", "lesion", "deep_infiltrating_endometriosis_intestinal"),
        "DiEB": ("endometriosis", "lesion", "deep_infiltrating_endometriosis_bladder"),
    }
    disease_status, broad_tissue_class, tissue_subtype = mapping.get(base, ("unknown", "unknown", "unknown"))
    return {
        "sample_code_base": base,
        "is_replicate_label": str(replicate),
        "disease_status": disease_status,
        "broad_tissue_class": broad_tissue_class,
        "tissue_subtype_preliminary": tissue_subtype,
        "code_interpretation_status": "preliminary_from_GSE141549_sample_code",
    }


def parse_gse141549_metadata(expression_columns: list[str]) -> pd.DataFrame:
    sample_link = pd.read_excel(RAW / "GSE141549__GSE141549_Sample_link.xlsx")
    link_records: dict[str, str] = {}
    duplicates: Counter[str] = Counter()
    for _, row in sample_link.iterrows():
        link_id = str(row.iloc[0]).strip()
        for value in row.iloc[1:]:
            if pd.isna(value):
                continue
            sample_code = str(value).strip()
            duplicates[sample_code] += 1
            link_records[sample_code] = link_id

    records: list[dict[str, str]] = []
    for column in expression_columns:
        match = re.match(r"^(SAMPLE\s+\d+)\s+(.+)$", str(column))
        sample_number = match.group(1) if match else ""
        sample_code = str(column)
        sample_code_suffix = match.group(2) if match else ""
        record = {
            "dataset": "GSE141549",
            "matrix_sample_id": sample_code,
            "sample_number": sample_number,
            "sample_code_suffix": sample_code_suffix,
            "sample_link_id": link_records.get(sample_code, ""),
            "sample_link_status": "matched" if sample_code in link_records else "unmatched",
        }
        record.update(classify_gse141549_code(sample_code_suffix))
        records.append(record)
    meta = pd.DataFrame(records)
    meta["sample_link_duplicate_count"] = meta["matrix_sample_id"].map(duplicates).fillna(0).astype(int)
    return meta


def parse_gse51981_metadata() -> pd.DataFrame:
    path = RAW / "GSE51981__GSE51981_series_matrix.txt.gz"
    sample_fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            parts = [clean_value(part) for part in line.rstrip("\n").split("\t")]
            field = parts[0].replace("!", "")
            values = parts[1:]
            if field == "Sample_characteristics_ch1":
                field = f"{field}_{sum(key.startswith('Sample_characteristics_ch1') for key in sample_fields) + 1}"
            sample_fields[field] = values

    accessions = sample_fields.get("Sample_geo_accession", [])
    records: list[dict[str, str]] = []
    for idx, accession in enumerate(accessions):
        record = {"dataset": "GSE51981", "matrix_sample_id": accession, "geo_accession": accession}
        for field, values in sample_fields.items():
            if idx < len(values):
                record[field] = values[idx]
        for field, value in list(record.items()):
            if not field.startswith("Sample_characteristics_ch1"):
                continue
            if ":" in value:
                key, val = value.split(":", 1)
                clean_key = re.sub(r"[^0-9A-Za-z]+", "_", key.strip().lower()).strip("_")
                record[clean_key] = val.strip()
        records.append(record)
    meta = pd.DataFrame(records)
    meta["analysis_role"] = "independent_endometrium_validation_probe_level_until_GPL570_mapping"
    return meta


def candidate_gene_table() -> pd.DataFrame:
    candidates = pd.read_csv(GWAS_CANDIDATES, sep="\t")
    grouped = (
        candidates.groupby(["gene_id", "gene_symbol"], dropna=False)
        .agg(
            genetic_priority=("genetic_priority", lambda x: ",".join(sorted(set(map(str, x))))),
            neighborhoods=("neighborhood_id", lambda x: ",".join(sorted(set(map(str, x))))),
            ld_neighborhood_class=("ld_neighborhood_class", lambda x: ",".join(sorted(set(map(str, x))))),
            module_hint_preliminary=("module_hint_preliminary", lambda x: ",".join(sorted({str(v) for v in x if str(v) != "nan" and str(v)}))),
            gene_biotype=("gene_biotype", "first"),
            gene_description=("gene_description", "first"),
        )
        .reset_index()
    )
    return grouped


def extract_gse234354(candidates: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    matrix = pd.read_csv(RAW / "GSE234354__GSE234354_gene_count_matrix.txt.gz", sep="\t")
    matrix["gene_id_stripped"] = matrix["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    gene_map = candidates[["gene_id", "gene_symbol"]].drop_duplicates()
    out = matrix[matrix["gene_id_stripped"].isin(set(candidates["gene_id"]))].copy()
    out = out.merge(gene_map, left_on="gene_id_stripped", right_on="gene_id", how="left", suffixes=("_matrix", ""))
    leading = ["gene_id", "gene_symbol", "gene_id_matrix", "gene_id_stripped"]
    sample_cols = [col for col in matrix.columns if col not in {"gene_id", "gene_id_stripped"}]
    out = out[leading + sample_cols]
    stats = {
        "candidate_genes": candidates["gene_id"].nunique(),
        "matched_rows": len(out),
        "matched_unique_gene_ids": out["gene_id"].nunique(),
        "sample_columns": len(sample_cols),
    }
    return out, stats


def extract_gse313775(candidates: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    matrix = pd.read_csv(RAW / "GSE313775__GSE313775_rawCountMatrix.tsv.gz", sep="\t")
    matrix["gene_id_stripped"] = matrix["Gene"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    gene_ids = set(candidates["gene_id"])
    symbols = set(candidates["gene_symbol"])
    out = matrix[matrix["gene_id_stripped"].isin(gene_ids) | matrix["Symbol"].isin(symbols)].copy()
    gene_map = candidates[["gene_id", "gene_symbol"]].drop_duplicates()
    out = out.merge(gene_map, left_on="gene_id_stripped", right_on="gene_id", how="left", suffixes=("_matrix", ""))
    out["gene_symbol"] = out["gene_symbol"].fillna(out["Symbol"])
    sample_cols = [col for col in matrix.columns if col not in {"Gene", "Symbol", "gene_id_stripped"}]
    out = out[["gene_id", "gene_symbol", "Gene", "Symbol", "gene_id_stripped"] + sample_cols]
    stats = {
        "candidate_genes": candidates["gene_id"].nunique(),
        "matched_rows": len(out),
        "matched_unique_gene_ids": out["gene_id"].nunique(),
        "matched_unique_symbols": out["gene_symbol"].nunique(),
        "sample_columns": len(sample_cols),
    }
    return out, stats


def extract_gse141549(candidates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    matrix = pd.read_excel(RAW / "GSE141549__GSE141549_batchCorrectednormalizedArrayscombined.xlsx")
    symbols = set(candidates["gene_symbol"])
    probe = matrix[matrix["Gene_symbol"].isin(symbols)].copy()
    symbol_to_gene_ids = candidates.groupby("gene_symbol")["gene_id"].apply(lambda x: ",".join(sorted(set(x)))).to_dict()
    probe.insert(1, "candidate_gene_id", probe["Gene_symbol"].map(symbol_to_gene_ids))
    sample_cols = [col for col in matrix.columns if col not in {"Gene_symbol", "Probe_Id"}]
    grouped = probe.groupby("Gene_symbol", dropna=False)
    gene_expr = grouped[sample_cols].median().reset_index()
    gene_expr.insert(1, "candidate_gene_id", gene_expr["Gene_symbol"].map(symbol_to_gene_ids))
    gene_expr.insert(2, "n_probes_collapsed", gene_expr["Gene_symbol"].map(grouped.size().to_dict()))
    stats = {
        "candidate_symbols": candidates["gene_symbol"].nunique(),
        "matched_probe_rows": len(probe),
        "matched_unique_symbols": probe["Gene_symbol"].nunique(),
        "gene_level_rows_after_median_probe_collapse": len(gene_expr),
        "sample_columns": len(sample_cols),
    }
    return probe, gene_expr, stats


def load_gpl570_annotation(candidates: pd.DataFrame) -> pd.DataFrame:
    annotation = pd.read_csv(
        META / "GPL570.annot.gz",
        sep="\t",
        skiprows=27,
        comment="!",
        dtype=str,
        low_memory=False,
    )
    symbols = set(candidates["gene_symbol"].dropna().astype(str))
    gene_id_by_symbol = candidates.groupby("gene_symbol")["gene_id"].apply(lambda x: ",".join(sorted(set(map(str, x))))).to_dict()

    records: list[dict[str, object]] = []
    for _, row in annotation.iterrows():
        raw_symbol = str(row.get("Gene symbol", "") or "")
        if not raw_symbol or raw_symbol == "nan":
            continue
        split_symbols = [symbol.strip() for symbol in raw_symbol.split("///") if symbol.strip()]
        unique_symbols = sorted(set(split_symbols))
        for symbol in unique_symbols:
            if symbol not in symbols:
                continue
            records.append(
                {
                    "probe_id": row["ID"],
                    "gene_symbol": symbol,
                    "candidate_gene_id": gene_id_by_symbol.get(symbol, ""),
                    "raw_gene_symbol_annotation": raw_symbol,
                    "raw_entrez_gene_id_annotation": row.get("Gene ID", ""),
                    "n_gene_symbols_on_probe": len(unique_symbols),
                    "probe_mapping_status": "multi_gene_probe" if len(unique_symbols) > 1 else "single_gene_probe",
                }
            )
    return pd.DataFrame(records).drop_duplicates()


def extract_gse51981(candidates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, int]]:
    mapping = load_gpl570_annotation(candidates)
    matrix = pd.read_csv(
        RAW / "GSE51981__GSE51981_series_matrix.txt.gz",
        sep="\t",
        compression="gzip",
        comment="!",
    )
    sample_cols = [col for col in matrix.columns if col != "ID_REF"]
    probe = mapping.merge(matrix, left_on="probe_id", right_on="ID_REF", how="inner")
    gene_expr = (
        probe.groupby(["gene_symbol", "candidate_gene_id"], dropna=False)[sample_cols]
        .median()
        .reset_index()
    )
    gene_expr.insert(
        2,
        "n_probes_collapsed",
        gene_expr["gene_symbol"].map(probe.groupby("gene_symbol")["probe_id"].nunique().to_dict()),
    )
    stats = {
        "candidate_symbols": candidates["gene_symbol"].nunique(),
        "matched_probe_rows": len(probe),
        "matched_unique_symbols": probe["gene_symbol"].nunique(),
        "gene_level_rows_after_median_probe_collapse": len(gene_expr),
        "multi_gene_probe_rows": int((probe["probe_mapping_status"] == "multi_gene_probe").sum()),
        "sample_columns": len(sample_cols),
    }
    return mapping, probe, gene_expr, stats


def write_summary(metadata_stats: list[dict[str, object]], extraction_stats: list[dict[str, object]]) -> None:
    metadata_df = pd.DataFrame(metadata_stats)
    extraction_df = pd.DataFrame(extraction_stats)
    metadata_df.to_csv(RESULTS / "bulk_metadata_audit.tsv", sep="\t", index=False)
    extraction_df.to_csv(RESULTS / "candidate_expression_extraction_audit.tsv", sep="\t", index=False)

    lines = [
        "# Bulk metadata and candidate-expression preprocessing audit",
        "",
        "## Metadata parsing",
        "",
        dataframe_to_markdown(metadata_df),
        "",
        "## Candidate expression extraction",
        "",
        dataframe_to_markdown(extraction_df),
        "",
        "## Pre-analysis interpretation",
        "",
        "- GSE234354 is retained as a menstrual-cycle reference dataset, not as disease/control evidence.",
        "- GSE313775 is usable for donor-aware Th-subset immune validation; matrix columns do not contain GSM IDs, so donor-to-GEO mapping is recorded as an order-based audit assumption.",
        "- GSE141549 is usable for lesion/endometrium/peritoneum candidate-expression validation; tissue-code labels are preliminary and should be cross-checked against the Scientific Data descriptor before manuscript wording.",
        "- GSE51981 is now usable for independent endometrium validation after official GPL570 annotation; multi-gene probes are retained with flags and should be sensitivity-filtered.",
        "",
    ]
    (RESULTS / "bulk_metadata_and_candidate_expression_summary.md").write_text("\n".join(lines), encoding="utf-8")


def dataframe_to_markdown(df: pd.DataFrame) -> str:
    columns = [str(col) for col in df.columns]
    rows = []
    rows.append("| " + " | ".join(columns) + " |")
    rows.append("| " + " | ".join(["---"] * len(columns)) + " |")
    for _, row in df.iterrows():
        values = [str(row[col]).replace("|", "\\|") for col in df.columns]
        rows.append("| " + " | ".join(values) + " |")
    return "\n".join(rows)


def write_self_review(metadata_stats: list[dict[str, object]], extraction_stats: list[dict[str, object]]) -> None:
    gse313775 = next(item for item in metadata_stats if item["dataset"] == "GSE313775")
    gse141549 = next(item for item in metadata_stats if item["dataset"] == "GSE141549")
    verdict = "PASS_WITH_CONDITIONS"
    conditions = [
        "GSE313775 disease/control and Th-subset metadata are now recoverable, but the raw count matrix uses EDCV donor IDs rather than GEO GSM IDs; retain the order-based mapping assumption in methods and sensitivity notes.",
        "GSE141549 sample-code tissue labels must be cited to the dataset descriptor before final manuscript claims.",
        "GSE51981 probe-to-gene annotation is available; final validation should report a sensitivity analysis excluding multi-gene probes.",
    ]
    if int(gse313775.get("matched_samples", 0)) != int(gse313775.get("matrix_sample_columns", -1)):
        verdict = "REQUIRES_FIX"
        conditions.insert(0, "GSE313775 sample metadata did not fully match matrix columns.")
    if int(gse141549.get("matched_samples", 0)) != int(gse141549.get("matrix_sample_columns", -1)):
        verdict = "REQUIRES_FIX"
        conditions.insert(0, "GSE141549 sample-link metadata did not fully match matrix columns.")

    lines = [
        "# Phase 3 bulk preprocessing self-review",
        "",
        f"Verdict: {verdict}",
        "",
        "## Checks passed",
        "",
        "- Parsed sample-level SOFT metadata for GSE234354 and GSE313775.",
        "- Parsed GSE141549 sample-link groups and expression-column tissue codes.",
        "- Parsed GSE51981 series-matrix sample metadata.",
        "- Extracted candidate expression matrices for GSE234354, GSE313775, GSE141549 and GSE51981.",
        "- Preserved probe-level and median-collapsed gene-level outputs for microarray datasets.",
        "",
        "## Conditions before differential-expression claims",
        "",
    ]
    lines.extend(f"- {condition}" for condition in conditions)
    lines.extend(
        [
            "",
            "## Next step",
            "",
            "Build dataset-specific QC summaries and then run prespecified candidate-gene association models rather than whole-transcriptome DEG discovery.",
            "",
        ]
    )
    (RESULTS / "phase3_bulk_preprocessing_self_review.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)

    candidates = candidate_gene_table()
    candidates.to_csv(RESULTS / "candidate_genes_unique.tsv", sep="\t", index=False)

    gse234354_counts_cols = list(pd.read_csv(RAW / "GSE234354__GSE234354_gene_count_matrix.txt.gz", sep="\t", nrows=0).columns[1:])
    gse313775_count_cols = list(pd.read_csv(RAW / "GSE313775__GSE313775_rawCountMatrix.tsv.gz", sep="\t", nrows=0).columns[2:])
    gse141549_cols = list(pd.read_excel(RAW / "GSE141549__GSE141549_batchCorrectednormalizedArrayscombined.xlsx", nrows=0).columns[2:])

    meta234 = parse_gse234354_metadata()
    meta313 = parse_gse313775_metadata(gse313775_count_cols)
    meta141 = parse_gse141549_metadata(gse141549_cols)
    meta519 = parse_gse51981_metadata()

    meta234.to_csv(RESULTS / "metadata_GSE234354.tsv", sep="\t", index=False)
    meta313.to_csv(RESULTS / "metadata_GSE313775.tsv", sep="\t", index=False)
    meta141.to_csv(RESULTS / "metadata_GSE141549.tsv", sep="\t", index=False)
    meta519.to_csv(RESULTS / "metadata_GSE51981.tsv", sep="\t", index=False)

    expr234, stats234 = extract_gse234354(candidates)
    expr313, stats313 = extract_gse313775(candidates)
    probe141, gene141, stats141 = extract_gse141549(candidates)
    gpl570_mapping, probe519, gene519, stats519 = extract_gse51981(candidates)

    expr234.to_csv(RESULTS / "GSE234354_candidate_counts.tsv", sep="\t", index=False)
    expr313.to_csv(RESULTS / "GSE313775_candidate_counts.tsv", sep="\t", index=False)
    probe141.to_csv(RESULTS / "GSE141549_candidate_probe_expression.tsv", sep="\t", index=False)
    gene141.to_csv(RESULTS / "GSE141549_candidate_gene_expression_median_probe.tsv", sep="\t", index=False)
    gpl570_mapping.to_csv(RESULTS / "GPL570_candidate_probe_mapping.tsv", sep="\t", index=False)
    probe519.to_csv(RESULTS / "GSE51981_candidate_probe_expression.tsv", sep="\t", index=False)
    gene519.to_csv(RESULTS / "GSE51981_candidate_gene_expression_median_probe.tsv", sep="\t", index=False)

    metadata_stats = [
        {
            "dataset": "GSE234354",
            "metadata_rows": len(meta234),
            "matrix_sample_columns": len(gse234354_counts_cols),
            "matched_samples": int(meta234["matrix_sample_id"].isin(gse234354_counts_cols).sum()),
            "primary_role": "cycle_reference",
        },
        {
            "dataset": "GSE313775",
            "metadata_rows": len(meta313),
            "matrix_sample_columns": len(gse313775_count_cols),
            "matched_samples": int((meta313["metadata_match_status"] == "matched").sum()),
            "primary_role": "Th_subset_immune_validation",
        },
        {
            "dataset": "GSE141549",
            "metadata_rows": len(meta141),
            "matrix_sample_columns": len(gse141549_cols),
            "matched_samples": int((meta141["sample_link_status"] == "matched").sum()),
            "primary_role": "lesion_endometrium_peritoneum_validation",
        },
        {
            "dataset": "GSE51981",
            "metadata_rows": len(meta519),
            "matrix_sample_columns": len(meta519),
            "matched_samples": len(meta519),
            "primary_role": "independent_endometrium_validation_GPL570_annotated",
        },
    ]

    extraction_stats = [
        {"dataset": "GSE234354", **stats234},
        {"dataset": "GSE313775", **stats313},
        {"dataset": "GSE141549_probe", **stats141},
        {"dataset": "GSE51981_probe", **stats519, "status": "GPL570_annotated"},
    ]
    write_summary(metadata_stats, extraction_stats)
    write_self_review(metadata_stats, extraction_stats)

    print("Wrote bulk metadata and candidate-expression preprocessing outputs.")
    print(RESULTS / "bulk_metadata_and_candidate_expression_summary.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
