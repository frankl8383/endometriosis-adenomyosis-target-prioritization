#!/usr/bin/env python3
"""Query PubMed literature evidence for stable target hypotheses."""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT_ROOT / "results"
INTEGRATION = RESULTS / "integration"
LIT = RESULTS / "literature"
PY_PKG_DIR = PROJECT_ROOT / "tools" / "py_packages"

if sys.version_info[:2] != (3, 12):
    raise SystemExit(
        "Use the bundled Python 3.12 runtime: "
        "python3 "
        "scripts/query_phase14_pubmed_literature.py"
    )

sys.path.insert(0, str(PY_PKG_DIR))

import pandas as pd  # noqa: E402


INPUT = INTEGRATION / "rank_stability_combined_matrix.tsv"
RAW_JSON = LIT / "phase14_pubmed_query_raw.json"
OUT_TABLE = LIT / "phase14_pubmed_literature_hits.tsv"
OUT_SUMMARY = LIT / "phase14_pubmed_literature_summary.md"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

TARGET_GENES = [
    "KDM1A",
    "KDR",
    "LY96",
    "PDGFRA",
    "ECE1",
    "STK38L",
    "SRD5A3",
    "KIT",
    "WNT4",
    "ITPR2",
    "C1QA",
    "HSPG2",
    "SSPN",
    "ESR1",
]

MECHANISM_QUERIES = {
    "VEGF_KDR_angiogenesis": '("VEGF"[Title/Abstract] OR "VEGFR2"[Title/Abstract] OR "KDR"[Title/Abstract]) AND (endometriosis[Title/Abstract] OR adenomyosis[Title/Abstract]) AND angiogenesis[Title/Abstract]',
    "PDGF_fibrosis": '("PDGF"[Title/Abstract] OR "PDGFRA"[Title/Abstract] OR "PDGFR"[Title/Abstract]) AND (endometriosis[Title/Abstract] OR adenomyosis[Title/Abstract]) AND (fibrosis[Title/Abstract] OR fibroblast[Title/Abstract] OR stromal[Title/Abstract])',
    "TLR4_LY96_inflammation": '("LY96"[Title/Abstract] OR "MD-2"[Title/Abstract] OR "TLR4"[Title/Abstract] OR "Toll-like receptor 4"[Title/Abstract]) AND (endometriosis[Title/Abstract] OR adenomyosis[Title/Abstract])',
    "mast_cell_KIT": '("KIT"[Title/Abstract] OR "c-Kit"[Title/Abstract] OR "mast cell"[Title/Abstract] OR "mast cells"[Title/Abstract]) AND (endometriosis[Title/Abstract] OR adenomyosis[Title/Abstract])',
    "WNT4_GWAS": '"WNT4"[Title/Abstract] AND (endometriosis[Title/Abstract] OR adenomyosis[Title/Abstract] OR uterus[Title/Abstract] OR endometrium[Title/Abstract]) AND (GWAS[Title/Abstract] OR genetic[Title/Abstract] OR risk[Title/Abstract] OR polymorphism[Title/Abstract])',
    "endothelin_ECE1": '("ECE1"[Title/Abstract] OR "endothelin converting enzyme"[Title/Abstract] OR endothelin[Title/Abstract]) AND (endometriosis[Title/Abstract] OR adenomyosis[Title/Abstract])',
    "KDM1A_epigenetic": '("KDM1A"[Title/Abstract] OR "LSD1"[Title/Abstract] OR "lysine-specific demethylase 1"[Title/Abstract]) AND (endometriosis[Title/Abstract] OR adenomyosis[Title/Abstract] OR endometrium[Title/Abstract])',
    "steroid_SRD5A3": '("SRD5A3"[Title/Abstract] OR steroidogenesis[Title/Abstract] OR "androgen metabolism"[Title/Abstract]) AND (endometriosis[Title/Abstract] OR adenomyosis[Title/Abstract] OR endometrium[Title/Abstract])',
    "complement_C1Q": '("C1QA"[Title/Abstract] OR "C1QB"[Title/Abstract] OR "C1QC"[Title/Abstract] OR complement[Title/Abstract]) AND (endometriosis[Title/Abstract] OR adenomyosis[Title/Abstract])',
    "ECM_HSPG2_SSPN": '("HSPG2"[Title/Abstract] OR perlecan[Title/Abstract] OR "SSPN"[Title/Abstract] OR sarcospan[Title/Abstract]) AND (endometriosis[Title/Abstract] OR adenomyosis[Title/Abstract] OR endometrium[Title/Abstract])',
    "KIT_CSF1R_pexidartinib": '(pexidartinib[Title/Abstract] OR PLX3397[Title/Abstract] OR "CSF1R"[Title/Abstract]) AND (KIT[Title/Abstract] OR "c-Kit"[Title/Abstract] OR KITL[Title/Abstract]) AND endometriosis[Title/Abstract]',
    "KDR_genetic_endometriosis": '(KDR[Title/Abstract] OR VEGFR2[Title/Abstract] OR "VEGFR-2"[Title/Abstract]) AND endometriosis[Title/Abstract] AND (variant[Title/Abstract] OR variants[Title/Abstract] OR locus[Title/Abstract] OR loci[Title/Abstract] OR GWAS[Title/Abstract] OR genetic[Title/Abstract])',
    "TLR4_title_endometriosis": '("TLR4"[Title] OR "Toll-like receptor 4"[Title] OR "Toll-like receptor system"[Title]) AND endometriosis[Title/Abstract]',
    "PDGF_endometrial_cell_function": '(PDGF[Title/Abstract] OR PDGFR[Title/Abstract]) AND endometrial[Title/Abstract] AND (proliferation[Title/Abstract] OR motility[Title/Abstract] OR stromal[Title/Abstract])',
}

GENE_ALIASES = {
    "KDM1A": ["KDM1A", "LSD1", "lysine-specific demethylase 1"],
    "KDR": ["KDR", "VEGFR2", "VEGFR-2"],
    "LY96": ["LY96", "MD-2"],
    "PDGFRA": ["PDGFRA", "PDGFR alpha", "PDGFR-alpha"],
    "ECE1": ["ECE1", "endothelin converting enzyme"],
    "STK38L": ["STK38L", "NDR2"],
    "SRD5A3": ["SRD5A3"],
    "KIT": ["KIT", "c-Kit"],
    "WNT4": ["WNT4"],
    "ITPR2": ["ITPR2", "IP3 receptor type 2"],
    "C1QA": ["C1QA", "C1q"],
    "HSPG2": ["HSPG2", "perlecan"],
    "SSPN": ["SSPN", "sarcospan"],
    "ESR1": ["ESR1", "estrogen receptor alpha"],
}


def fetch_json(url: str, timeout: int = 30) -> dict[str, Any]:
    with urllib.request.urlopen(url, timeout=timeout) as response:
        return json.loads(response.read().decode())


def fetch_text(url: str, timeout: int = 30) -> str:
    with urllib.request.urlopen(url, timeout=timeout) as response:
        return response.read().decode()


def esearch(query: str, retmax: int, sleep: float) -> tuple[int, list[str]]:
    params = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": str(retmax),
        "sort": "relevance",
    }
    url = f"{EUTILS}/esearch.fcgi?{urllib.parse.urlencode(params)}"
    data = fetch_json(url)
    time.sleep(sleep)
    result = data.get("esearchresult", {})
    count = int(result.get("count", 0))
    ids = list(result.get("idlist", []))
    return count, ids


def efetch_details(pmids: list[str], sleep: float) -> list[dict[str, Any]]:
    if not pmids:
        return []
    params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"}
    url = f"{EUTILS}/efetch.fcgi?{urllib.parse.urlencode(params)}"
    text = fetch_text(url)
    time.sleep(sleep)
    root = ET.fromstring(text)
    records: list[dict[str, Any]] = []
    for article in root.findall(".//PubmedArticle"):
        pmid = "".join(article.findtext(".//PMID") or "").strip()
        title = "".join(article.findtext(".//ArticleTitle") or "").strip()
        journal = "".join(article.findtext(".//Journal/Title") or "").strip()
        year = (
            article.findtext(".//JournalIssue/PubDate/Year")
            or article.findtext(".//ArticleDate/Year")
            or article.findtext(".//PubDate/MedlineDate")
            or ""
        )
        doi = ""
        for article_id in article.findall(".//ArticleId"):
            if article_id.attrib.get("IdType") == "doi":
                doi = "".join(article_id.itertext()).strip()
                break
        abstract_parts = [" ".join(node.itertext()).strip() for node in article.findall(".//Abstract/AbstractText")]
        abstract = " ".join(part for part in abstract_parts if part)
        authors = []
        for author in article.findall(".//Author")[:5]:
            lastname = author.findtext("LastName") or ""
            initials = author.findtext("Initials") or ""
            collective = author.findtext("CollectiveName") or ""
            if collective:
                authors.append(collective)
            elif lastname:
                authors.append(f"{lastname} {initials}".strip())
        records.append(
            {
                "pmid": pmid,
                "title": title,
                "journal": journal,
                "year": str(year)[:4],
                "doi": doi,
                "authors": "; ".join(authors),
                "abstract": abstract,
                "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
            }
        )
    return records


def build_queries(stable_genes: list[str]) -> dict[str, dict[str, str]]:
    queries: dict[str, dict[str, str]] = {}
    for gene in stable_genes:
        aliases = GENE_ALIASES.get(gene, [gene])
        gene_terms = " OR ".join(f'"{alias}"[Title/Abstract]' for alias in aliases)
        queries[f"gene_disease_{gene}"] = {
            "category": "gene_disease",
            "gene": gene,
            "query": f"({gene_terms}) AND (endometriosis[Title/Abstract] OR adenomyosis[Title/Abstract] OR endometrium[Title/Abstract])",
        }
    for name, query in MECHANISM_QUERIES.items():
        queries[f"mechanism_{name}"] = {"category": "mechanism", "gene": "", "query": query}
    queries["disease_target_prioritization"] = {
        "category": "framework",
        "gene": "",
        "query": "(endometriosis[Title/Abstract] OR adenomyosis[Title/Abstract]) AND (\"drug target\"[Title/Abstract] OR \"therapeutic target\"[Title/Abstract] OR \"target prioritization\"[Title/Abstract] OR \"multi-omics\"[Title/Abstract])",
    }
    return queries


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retmax", type=int, default=8)
    parser.add_argument("--sleep", type=float, default=0.35)
    args = parser.parse_args()

    LIT.mkdir(parents=True, exist_ok=True)
    ranked = pd.read_csv(INPUT, sep="\t", keep_default_na=False)
    stable = ranked[
        ranked["rank_stability_class"].isin(["rank_stable_target_hypothesis", "moderately_stable_target_hypothesis"])
    ]
    stable_genes = [gene for gene in stable["gene_symbol"].tolist() if gene] or TARGET_GENES[:9]
    for gene in TARGET_GENES:
        if gene not in stable_genes:
            stable_genes.append(gene)

    queries = build_queries(stable_genes)
    raw: dict[str, Any] = {"queries": queries, "records": {}}
    rows: list[dict[str, Any]] = []
    for key, spec in queries.items():
        count, pmids = esearch(spec["query"], args.retmax, args.sleep)
        details = efetch_details(pmids, args.sleep)
        raw["records"][key] = {"count": count, "pmids": pmids, "details": details}
        for rank, record in enumerate(details, start=1):
            item = {
                "query_key": key,
                "query_category": spec["category"],
                "query_gene": spec["gene"],
                "query": spec["query"],
                "pubmed_total_count": count,
                "query_rank": rank,
            }
            item.update(record)
            rows.append(item)

    RAW_JSON.write_text(json.dumps(raw, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    table = pd.DataFrame(rows)
    if table.empty:
        table = pd.DataFrame(
            columns=[
                "query_key",
                "query_category",
                "query_gene",
                "query",
                "pubmed_total_count",
                "query_rank",
                "pmid",
                "title",
                "journal",
                "year",
                "doi",
                "authors",
                "abstract",
                "pubmed_url",
            ]
        )
    table.to_csv(OUT_TABLE, sep="\t", index=False)

    counts = table.groupby(["query_key", "query_category", "query_gene"], dropna=False)["pmid"].count().reset_index(name="fetched_records")
    lines = [
        "# Phase 14 PubMed literature query summary",
        "",
        f"- Query count: {len(queries)}",
        f"- Fetched citation rows: {len(table)}",
        f"- Stable genes included: {', '.join(stable_genes)}",
        "",
        "| query | category | gene | PubMed total | fetched |",
        "| --- | --- | --- | ---: | ---: |",
    ]
    for key, spec in queries.items():
        fetched = int(counts.loc[counts["query_key"] == key, "fetched_records"].iloc[0]) if key in set(counts["query_key"]) else 0
        total = int(raw["records"][key]["count"])
        lines.append(f"| {key} | {spec['category']} | {spec['gene']} | {total} | {fetched} |")
    lines.extend(
        [
            "",
            "The table contains PubMed metadata and abstracts for directionality/safety review. It is not a systematic review and should be treated as targeted evidence retrieval.",
        ]
    )
    OUT_SUMMARY.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {OUT_TABLE}")
    print(f"Wrote {OUT_SUMMARY}")
    print(f"Wrote {RAW_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
