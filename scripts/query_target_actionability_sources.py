#!/usr/bin/env python3
"""Query and score target actionability evidence for candidate genes.

The script caches Open Targets, DGIdb and ChEMBL records so downstream
prioritisation can be reproduced without repeatedly hitting live APIs.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT_ROOT / "results"
INTEGRATION = RESULTS / "integration"
OUT = RESULTS / "druggability"
PY_PKG_DIR = PROJECT_ROOT / "tools" / "py_packages"

if sys.version_info[:2] != (3, 12):
    raise SystemExit(
        "Use the bundled Python 3.12 runtime: "
        "python3 "
        "scripts/query_target_actionability_sources.py"
    )

sys.path.insert(0, str(PY_PKG_DIR))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


INPUT = INTEGRATION / "cross_disease_singlecell_localization_matrix.tsv"
OT_RAW = OUT / "raw_open_targets_target_actionability.jsonl"
DGIDB_RAW = OUT / "raw_dgidb_gene_actionability.jsonl"
CHEMBL_RAW = OUT / "raw_chembl_target_actionability.jsonl"
SCORES = OUT / "target_actionability_scores.tsv"
SUMMARY = OUT / "target_actionability_summary.md"
REVIEW = OUT / "phase11_druggability_actionability_self_review.md"
SOURCE_AUDIT = OUT / "actionability_source_api_audit.md"

OPEN_TARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"
DGIDB_URL = "https://dgidb.org/api/graphql"
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"

OT_QUERY = """
query target($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    approvedName
    biotype
    isEssential
    targetClass { label level }
    tractability { label modality value }
    drugAndClinicalCandidates {
      count
      rows {
        id
        maxClinicalStage
        drug { id name drugType maximumClinicalStage }
      }
    }
    safetyLiabilities {
      event
      eventId
      datasource
      effects { direction dosing }
      biosamples { tissueLabel cellLabel cellFormat }
      studies { name type description }
      literature
      url
    }
    depMapEssentiality { tissueId tissueName }
  }
}
"""

DGIDB_QUERY = """
query q($genes: [String!]!) {
  genes(names: $genes, first: 10) {
    totalCount
    nodes {
      name
      conceptId
      longName
      geneCategories { name }
      geneCategoriesWithSources { name sourceNames }
    }
  }
  interactions(geneNames: $genes, first: 500) {
    totalCount
    nodes {
      id
      interactionScore
      evidenceScore
      gene { name conceptId }
      drug { name conceptId approved immunotherapy antiNeoplastic }
      interactionTypes { type directionality }
      sources { sourceDbName }
      publications { pmid }
    }
  }
}
"""

STAGE_RANK = {
    "UNKNOWN": 0,
    "NA": 0,
    "NOT_PROVIDED": 0,
    "EARLY_PHASE1": 0.5,
    "PHASE_1": 1,
    "PHASE_1/PHASE_2": 1.5,
    "PHASE_2": 2,
    "PHASE_2/PHASE_3": 2.5,
    "PHASE_3": 3,
    "PHASE_4": 4,
    "APPROVAL": 4,
}

HIGH_SYSTEMIC_RISK_SYMBOLS = {
    "ANGPT1",
    "ANGPT2",
    "ESR1",
    "FLT1",
    "FLT4",
    "KDR",
    "PGR",
    "PRLR",
    "TEK",
    "TGFB1",
    "TGFBR1",
    "TGFBR2",
    "VEGFA",
}

STRUCTURAL_MARKER_TERMS = [
    "basement membrane",
    "collagen",
    "extracellular matrix",
    "heparan sulfate proteoglycan",
    "proteoglycan",
    "sarcomere",
    "structural constituent",
]

ACTIONABLE_CLASS_TERMS = [
    "channel",
    "enzyme",
    "external side",
    "g protein-coupled receptor",
    "kinase",
    "ligand",
    "membrane",
    "phosphatase",
    "plasma membrane",
    "protease",
    "receptor",
    "secreted",
    "surface",
    "transporter",
]


def post_graphql(url: str, query: str, variables: dict[str, Any], timeout: int = 60) -> dict[str, Any]:
    payload = json.dumps({"query": query, "variables": variables}).encode()
    request = urllib.request.Request(
        url,
        data=payload,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            data = json.loads(response.read().decode())
    except urllib.error.HTTPError as exc:
        body = exc.read().decode(errors="replace")
        return {"_error": f"HTTPError {exc.code}", "_body": body[:5000]}
    except Exception as exc:  # noqa: BLE001
        return {"_error": type(exc).__name__, "_body": str(exc)}
    return data


def get_json(url: str, timeout: int = 60) -> dict[str, Any]:
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return json.loads(response.read().decode())
    except urllib.error.HTTPError as exc:
        body = exc.read().decode(errors="replace")
        return {"_error": f"HTTPError {exc.code}", "_body": body[:5000], "_url": url}
    except Exception as exc:  # noqa: BLE001
        return {"_error": type(exc).__name__, "_body": str(exc), "_url": url}


def write_jsonl(path: Path, records: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=True, sort_keys=True) + "\n")


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    records: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                records.append(json.loads(line))
    return records


def stage_value(stage: object) -> float:
    if stage is None:
        return 0.0
    if isinstance(stage, (int, float)) and not math.isnan(float(stage)):
        return float(stage)
    key = str(stage).strip().upper().replace(" ", "_")
    return float(STAGE_RANK.get(key, 0.0))


def flatten_strings(value: Any) -> list[str]:
    strings: list[str] = []
    if value is None:
        return strings
    if isinstance(value, str):
        return [value]
    if isinstance(value, dict):
        for item in value.values():
            strings.extend(flatten_strings(item))
        return strings
    if isinstance(value, list):
        for item in value:
            strings.extend(flatten_strings(item))
        return strings
    return [str(value)]


def true_tractability_labels(target: dict[str, Any]) -> list[str]:
    labels: list[str] = []
    for item in target.get("tractability") or []:
        if item.get("value") is True:
            labels.append(f"{item.get('modality', '')}:{item.get('label', '')}")
    return labels


def class_labels(target: dict[str, Any], dgidb: dict[str, Any]) -> list[str]:
    labels: list[str] = []
    for item in target.get("targetClass") or []:
        labels.append(str(item.get("label", "")))
    for node in (dgidb.get("data") or {}).get("genes", {}).get("nodes", []) or []:
        for category in node.get("geneCategories") or []:
            labels.append(str(category.get("name", "")))
        for category in node.get("geneCategoriesWithSources") or []:
            labels.append(str(category.get("name", "")))
    return sorted({label for label in labels if label})


def select_chembl_human_target(symbol: str, search_result: dict[str, Any]) -> dict[str, Any] | None:
    candidates = []
    symbol_upper = symbol.upper()
    for target in search_result.get("targets") or []:
        if target.get("organism") != "Homo sapiens":
            continue
        synonyms = []
        for component in target.get("target_components") or []:
            for synonym in component.get("target_component_synonyms") or []:
                synonyms.append(str(synonym.get("component_synonym", "")).upper())
        pref = str(target.get("pref_name", "")).upper()
        target_type = str(target.get("target_type", ""))
        exact_synonym = symbol_upper in synonyms
        single_protein = target_type == "SINGLE PROTEIN"
        score = float(target.get("score") or 0)
        rank = (
            100 if exact_synonym and single_protein else 0,
            50 if exact_synonym else 0,
            25 if single_protein else 0,
            score,
            5 if symbol_upper in pref else 0,
        )
        candidates.append((rank, target))
    if not candidates:
        return None
    return sorted(candidates, key=lambda item: item[0], reverse=True)[0][1]


def query_open_targets(candidates: pd.DataFrame, sleep_seconds: float) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for _, row in candidates.iterrows():
        gene_id = str(row["gene_id"]).split(".")[0]
        symbol = str(row["gene_symbol"])
        data = post_graphql(OPEN_TARGETS_URL, OT_QUERY, {"ensemblId": gene_id})
        records.append({"gene_id": gene_id, "gene_symbol": symbol, "response": data})
        time.sleep(sleep_seconds)
    return records


def query_dgidb(candidates: pd.DataFrame, sleep_seconds: float) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for _, row in candidates.iterrows():
        symbol = str(row["gene_symbol"])
        gene_id = str(row["gene_id"]).split(".")[0]
        if not symbol.strip():
            data = {"_skipped": "missing_gene_symbol"}
        else:
            data = post_graphql(DGIDB_URL, DGIDB_QUERY, {"genes": [symbol]})
        records.append({"gene_id": gene_id, "gene_symbol": symbol, "response": data})
        time.sleep(sleep_seconds)
    return records


def query_chembl(candidates: pd.DataFrame, sleep_seconds: float) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for _, row in candidates.iterrows():
        symbol = str(row["gene_symbol"])
        gene_id = str(row["gene_id"]).split(".")[0]
        if not symbol.strip():
            records.append(
                {
                    "gene_id": gene_id,
                    "gene_symbol": symbol,
                    "search": {"_skipped": "missing_gene_symbol"},
                    "selected_target": None,
                    "mechanisms": {},
                }
            )
            continue
        search_url = f"{CHEMBL_BASE}/target/search.json?q={urllib.parse.quote(symbol)}&limit=100"
        search_result = get_json(search_url)
        selected = None if "_error" in search_result else select_chembl_human_target(symbol, search_result)
        mechanisms: dict[str, Any] = {}
        if selected is not None:
            target_id = selected.get("target_chembl_id")
            mech_url = f"{CHEMBL_BASE}/mechanism.json?target_chembl_id={urllib.parse.quote(str(target_id))}&limit=1000"
            mechanisms = get_json(mech_url)
            time.sleep(sleep_seconds)
        records.append(
            {
                "gene_id": gene_id,
                "gene_symbol": symbol,
                "search": search_result,
                "selected_target": selected,
                "mechanisms": mechanisms,
            }
        )
        time.sleep(sleep_seconds)
    return records


def summarize_open_targets(record: dict[str, Any]) -> dict[str, Any]:
    response = record.get("response") or {}
    target = ((response.get("data") or {}).get("target")) or {}
    if not target:
        return {
            "ot_target_found": False,
            "ot_error": response.get("_error", ""),
            "ot_target_class_labels": "",
            "ot_true_tractability_labels": "",
            "ot_drug_candidate_count": 0,
            "ot_max_clinical_stage": 0,
            "ot_safety_liability_count": 0,
            "ot_depmap_tissue_row_count": 0,
            "ot_is_essential": "",
        }

    candidates = (target.get("drugAndClinicalCandidates") or {}).get("rows") or []
    stage_values = [stage_value(item.get("maxClinicalStage")) for item in candidates]
    drug_stage_values = [stage_value((item.get("drug") or {}).get("maximumClinicalStage")) for item in candidates]
    all_stage_values = stage_values + drug_stage_values
    classes = [str(item.get("label", "")) for item in target.get("targetClass") or [] if item.get("label")]
    safety = target.get("safetyLiabilities") or []
    depmap = target.get("depMapEssentiality") or []
    return {
        "ot_target_found": True,
        "ot_error": "",
        "ot_approved_symbol": target.get("approvedSymbol", ""),
        "ot_approved_name": target.get("approvedName", ""),
        "ot_biotype": target.get("biotype", ""),
        "ot_target_class_labels": "|".join(sorted(set(classes))),
        "ot_true_tractability_labels": "|".join(true_tractability_labels(target)),
        "ot_drug_candidate_count": int((target.get("drugAndClinicalCandidates") or {}).get("count") or len(candidates)),
        "ot_max_clinical_stage": max(all_stage_values) if all_stage_values else 0,
        "ot_safety_liability_count": len(safety),
        "ot_safety_events": "|".join(sorted({str(item.get("event", "")) for item in safety if item.get("event")})),
        "ot_depmap_tissue_row_count": len(depmap),
        "ot_depmap_tissue_rows": "|".join(sorted({str(item.get("tissueName", "")) for item in depmap if item.get("tissueName")})),
        "ot_is_essential": target.get("isEssential") if target.get("isEssential") is not None else "",
    }


def summarize_dgidb(record: dict[str, Any]) -> dict[str, Any]:
    response = record.get("response") or {}
    if response.get("_skipped"):
        return {
            "dgidb_gene_found": False,
            "dgidb_query_status": str(response.get("_skipped")),
            "dgidb_error": "",
            "dgidb_gene_category_labels": "",
            "dgidb_gene_categories_with_sources": "",
            "dgidb_interaction_total_count": 0,
            "dgidb_interaction_returned_count": 0,
            "dgidb_approved_interaction_count": 0,
            "dgidb_antineoplastic_interaction_count": 0,
            "dgidb_immunotherapy_interaction_count": 0,
            "dgidb_interaction_type_labels": "",
            "dgidb_source_labels": "",
            "dgidb_max_evidence_score": 0,
            "dgidb_max_interaction_score": 0,
        }
    data = response.get("data") or {}
    genes = data.get("genes") or {}
    interactions = data.get("interactions") or {}
    nodes = interactions.get("nodes") or []
    gene_nodes = genes.get("nodes") or []
    categories: list[str] = []
    categories_with_sources: list[str] = []
    for node in gene_nodes:
        for category in node.get("geneCategories") or []:
            categories.append(str(category.get("name", "")))
        for category in node.get("geneCategoriesWithSources") or []:
            name = str(category.get("name", ""))
            sources = ",".join(category.get("sourceNames") or [])
            categories_with_sources.append(f"{name}:{sources}" if sources else name)
    approved_count = sum(1 for item in nodes if (item.get("drug") or {}).get("approved") is True)
    anti_neoplastic_count = sum(1 for item in nodes if (item.get("drug") or {}).get("antiNeoplastic") is True)
    immunotherapy_count = sum(1 for item in nodes if (item.get("drug") or {}).get("immunotherapy") is True)
    directions: list[str] = []
    sources: list[str] = []
    for item in nodes:
        for interaction_type in item.get("interactionTypes") or []:
            text = str(interaction_type.get("type", ""))
            direction = str(interaction_type.get("directionality", ""))
            directions.append(f"{text}:{direction}" if direction else text)
        for source in item.get("sources") or []:
            sources.append(str(source.get("sourceDbName", "")))
    evidence_scores = [float(item.get("evidenceScore") or 0) for item in nodes]
    interaction_scores = [float(item.get("interactionScore") or 0) for item in nodes]
    return {
        "dgidb_gene_found": bool(gene_nodes),
        "dgidb_query_status": "ok" if not response.get("_error") else "error",
        "dgidb_error": response.get("_error", ""),
        "dgidb_gene_category_labels": "|".join(sorted(set(filter(None, categories)))),
        "dgidb_gene_categories_with_sources": "|".join(sorted(set(filter(None, categories_with_sources)))),
        "dgidb_interaction_total_count": int(interactions.get("totalCount") or 0),
        "dgidb_interaction_returned_count": len(nodes),
        "dgidb_approved_interaction_count": approved_count,
        "dgidb_antineoplastic_interaction_count": anti_neoplastic_count,
        "dgidb_immunotherapy_interaction_count": immunotherapy_count,
        "dgidb_interaction_type_labels": "|".join(sorted(set(filter(None, directions)))),
        "dgidb_source_labels": "|".join(sorted(set(filter(None, sources)))),
        "dgidb_max_evidence_score": max(evidence_scores) if evidence_scores else 0,
        "dgidb_max_interaction_score": max(interaction_scores) if interaction_scores else 0,
    }


def summarize_chembl(record: dict[str, Any]) -> dict[str, Any]:
    if not str(record.get("gene_symbol", "")).strip():
        return {
            "chembl_target_found": False,
            "chembl_query_status": "missing_gene_symbol",
            "chembl_error": "",
            "chembl_target_id": "",
            "chembl_target_pref_name": "",
            "chembl_target_type": "",
            "chembl_mechanism_count": 0,
            "chembl_returned_mechanism_count": 0,
            "chembl_max_phase": 0,
            "chembl_action_types": "",
            "chembl_mechanism_molecule_count": 0,
        }
    selected = record.get("selected_target")
    mechanisms = record.get("mechanisms") or {}
    mechanism_rows = mechanisms.get("mechanisms") or []
    phases = [float(item.get("max_phase") or 0) for item in mechanism_rows]
    action_types = [str(item.get("action_type", "")) for item in mechanism_rows]
    molecule_ids = [str(item.get("molecule_chembl_id", "")) for item in mechanism_rows if item.get("molecule_chembl_id")]
    return {
        "chembl_target_found": selected is not None,
        "chembl_query_status": "ok" if not ((record.get("search") or {}).get("_error") or mechanisms.get("_error")) else "error",
        "chembl_error": (record.get("search") or {}).get("_error", "") or mechanisms.get("_error", ""),
        "chembl_target_id": "" if selected is None else selected.get("target_chembl_id", ""),
        "chembl_target_pref_name": "" if selected is None else selected.get("pref_name", ""),
        "chembl_target_type": "" if selected is None else selected.get("target_type", ""),
        "chembl_mechanism_count": int((mechanisms.get("page_meta") or {}).get("total_count") or len(mechanism_rows)),
        "chembl_returned_mechanism_count": len(mechanism_rows),
        "chembl_max_phase": max(phases) if phases else 0,
        "chembl_action_types": "|".join(sorted(set(filter(None, action_types)))),
        "chembl_mechanism_molecule_count": len(set(molecule_ids)),
    }


def has_any_label(labels: list[str] | str, terms: list[str]) -> bool:
    text = " ".join(labels) if isinstance(labels, list) else str(labels)
    text = text.lower()
    return any(term in text for term in terms)


def calculate_scores(row: pd.Series) -> dict[str, Any]:
    all_labels = "|".join(
        [
            str(row.get("ot_target_class_labels", "")),
            str(row.get("ot_true_tractability_labels", "")),
            str(row.get("dgidb_gene_category_labels", "")),
            str(row.get("dgidb_interaction_type_labels", "")),
            str(row.get("chembl_action_types", "")),
        ]
    )
    actionable_class = has_any_label(all_labels, ACTIONABLE_CLASS_TERMS)
    approved_tractability = "Approved Drug" in str(row.get("ot_true_tractability_labels", ""))
    clinical_tractability = any(
        token in str(row.get("ot_true_tractability_labels", ""))
        for token in ["Advanced Clinical", "Phase 1 Clinical", "High-Quality Ligand", "Structure with Ligand", "Druggable Family"]
    )
    approved_drug_count = int(row.get("dgidb_approved_interaction_count", 0) or 0)
    ot_candidate_count = int(row.get("ot_drug_candidate_count", 0) or 0)
    dgidb_interaction_count = int(row.get("dgidb_interaction_total_count", 0) or 0)
    chembl_mechanism_count = int(row.get("chembl_mechanism_count", 0) or 0)
    direct_clinical_phase = max(
        float(row.get("ot_max_clinical_stage", 0) or 0),
        float(row.get("chembl_max_phase", 0) or 0),
    )

    target_class_score = 0
    if approved_tractability or direct_clinical_phase >= 4:
        target_class_score = 5
    elif clinical_tractability or "CLINICALLY ACTIONABLE" in str(row.get("dgidb_gene_category_labels", "")):
        target_class_score = 4
    elif actionable_class:
        target_class_score = 3
    elif str(row.get("gene_biotype", "")) == "protein_coding":
        target_class_score = 1

    drug_precedent_score = 0
    if direct_clinical_phase >= 4:
        drug_precedent_score = 5
    elif direct_clinical_phase >= 3:
        drug_precedent_score = 4
    elif direct_clinical_phase >= 2 or (approved_drug_count > 0 and row.get("chembl_target_found") is True):
        drug_precedent_score = 3
    elif direct_clinical_phase >= 1 or approved_drug_count > 0 or chembl_mechanism_count > 0:
        drug_precedent_score = 2
    elif dgidb_interaction_count > 0 or ot_candidate_count > 0:
        drug_precedent_score = 1

    strong_sources = 0
    if row.get("ot_target_found") is True and (approved_tractability or ot_candidate_count > 0 or actionable_class):
        strong_sources += 1
    if dgidb_interaction_count > 0 or "DRUGGABLE GENOME" in str(row.get("dgidb_gene_category_labels", "")):
        strong_sources += 1
    if row.get("chembl_target_found") is True and (chembl_mechanism_count > 0 or direct_clinical_phase > 0):
        strong_sources += 1
    evidence_convergence_score = min(4, strong_sources + (1 if strong_sources >= 2 else 0))

    safety_count = int(row.get("ot_safety_liability_count", 0) or 0)
    is_essential = str(row.get("ot_is_essential", "")).lower() == "true"
    safety_score = 4
    if is_essential:
        safety_score -= 2
    if safety_count >= 10:
        safety_score -= 2
    elif safety_count >= 3:
        safety_score -= 1
    if row.get("gene_symbol") in HIGH_SYSTEMIC_RISK_SYMBOLS:
        safety_score -= 1
    if safety_count == 0 and row.get("ot_is_essential") == "":
        safety_score = min(safety_score, 2)
    safety_score = max(0, min(4, safety_score))

    localization_axis = str(row.get("dominant_cross_disease_axis", ""))
    local_feasibility_score = 0
    if actionable_class and localization_axis in {"fibrovascular", "epithelial", "immune", "mixed"}:
        local_feasibility_score = 2
    elif actionable_class or direct_clinical_phase > 0:
        local_feasibility_score = 1

    raw_score = (
        target_class_score
        + drug_precedent_score
        + evidence_convergence_score
        + safety_score
        + local_feasibility_score
    )

    penalty_reasons: list[str] = []
    penalty = 0
    marker_text = " ".join(
        [
            str(row.get("gene_description", "")),
            str(row.get("module_hint_preliminary", "")),
            str(row.get("ot_target_class_labels", "")),
        ]
    ).lower()
    structural_marker = any(term in marker_text for term in STRUCTURAL_MARKER_TERMS)
    if structural_marker and direct_clinical_phase < 1 and ot_candidate_count == 0 and chembl_mechanism_count == 0:
        penalty += 5
        penalty_reasons.append("structural_or_ecm_readout_without_direct_drug_precedent")
    if row.get("gene_symbol") in HIGH_SYSTEMIC_RISK_SYMBOLS and direct_clinical_phase >= 1:
        penalty += 2
        penalty_reasons.append("systemic_modulation_safety_caution")
    if is_essential:
        penalty += 2
        penalty_reasons.append("open_targets_common_essentiality_caution")
    if safety_count >= 10:
        penalty += 2
        penalty_reasons.append("multiple_safety_liabilities")
    if row.get("chembl_target_found") is False and row.get("ot_target_found") is False and row.get("dgidb_gene_found") is False:
        penalty += 2
        penalty_reasons.append("no_external_actionability_record")

    net_score = max(0, min(20, raw_score - penalty))
    if net_score >= 16 and penalty <= 2:
        actionability_class = "high_actionability"
    elif net_score >= 16:
        actionability_class = "high_druggability_with_major_safety_caution"
    elif net_score >= 12:
        actionability_class = "moderate_actionability"
    elif net_score >= 8:
        actionability_class = "limited_actionability"
    else:
        actionability_class = "minimal_actionability"

    if chembl_mechanism_count > 0:
        direction = str(row.get("chembl_action_types", ""))
    elif dgidb_interaction_count > 0:
        direction = str(row.get("dgidb_interaction_type_labels", ""))
    else:
        direction = "no_drug_direction_evidence"

    guardrail = "actionability_score_not_a_therapeutic_claim_requires_expression_directionality_and_local_delivery_review"
    if penalty_reasons:
        guardrail += ";penalties=" + ",".join(penalty_reasons)

    return {
        "target_class_score_5": target_class_score,
        "drug_precedent_score_5": drug_precedent_score,
        "evidence_convergence_score_4": evidence_convergence_score,
        "safety_acceptability_score_4": safety_score,
        "local_repurposing_feasibility_score_2": local_feasibility_score,
        "raw_druggability_score_20": raw_score,
        "druggability_penalty_points": penalty,
        "druggability_penalty_reasons": "|".join(penalty_reasons),
        "druggability_score_20": round(net_score, 3),
        "actionability_class": actionability_class,
        "drug_direction_evidence": direction,
        "actionability_guardrail": guardrail,
    }


def build_scores(
    candidates: pd.DataFrame,
    ot_records: list[dict[str, Any]],
    dgidb_records: list[dict[str, Any]],
    chembl_records: list[dict[str, Any]],
) -> pd.DataFrame:
    ot_by_gene = {record["gene_id"]: summarize_open_targets(record) for record in ot_records}
    dgidb_by_gene = {record["gene_id"]: summarize_dgidb(record) for record in dgidb_records}
    chembl_by_gene = {record["gene_id"]: summarize_chembl(record) for record in chembl_records}

    rows: list[dict[str, Any]] = []
    for _, candidate in candidates.iterrows():
        gene_id = str(candidate["gene_id"]).split(".")[0]
        row = candidate.to_dict()
        row.update(ot_by_gene.get(gene_id, summarize_open_targets({})))
        row.update(dgidb_by_gene.get(gene_id, summarize_dgidb({})))
        row.update(chembl_by_gene.get(gene_id, summarize_chembl({})))
        row.update(calculate_scores(pd.Series(row)))
        row["final_target_priority_score_100_pre_rank_stability"] = round(
            float(row.get("pre_druggability_biologic_evidence_score_80", 0) or 0)
            + float(row.get("druggability_score_20", 0) or 0),
            3,
        )
        rows.append(row)
    return pd.DataFrame(rows)


def write_summary(scores: pd.DataFrame, args: argparse.Namespace) -> None:
    class_counts = scores["actionability_class"].value_counts().to_dict()
    top = scores.sort_values("final_target_priority_score_100_pre_rank_stability", ascending=False).head(20)
    source_errors = {
        "open_targets": int((scores["ot_error"].astype(str) != "").sum()),
        "dgidb": int((scores["dgidb_error"].astype(str) != "").sum()),
        "chembl": int((scores["chembl_error"].astype(str) != "").sum()),
    }
    lines = [
        "# Target actionability summary",
        "",
        f"- Input candidates: {len(scores)}",
        f"- Refresh live APIs: {args.refresh}",
        f"- Limited run: {args.limit if args.limit else 'no'}",
        f"- Source errors by summarized target: {source_errors}",
        f"- Actionability classes: {class_counts}",
        "",
        "## Top preliminary targets after adding druggability",
        "",
        "| rank | gene | final_pre_stability_score | druggability | class | main axis | penalties |",
        "| ---: | --- | ---: | ---: | --- | --- | --- |",
    ]
    for idx, (_, row) in enumerate(top.iterrows(), start=1):
        lines.append(
            "| "
            + " | ".join(
                [
                    str(idx),
                    str(row["gene_symbol"]),
                    f"{float(row['final_target_priority_score_100_pre_rank_stability']):.3f}",
                    f"{float(row['druggability_score_20']):.3f}",
                    str(row["actionability_class"]),
                    str(row.get("dominant_cross_disease_axis", "")),
                    str(row.get("druggability_penalty_reasons", "")),
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "## Guardrails",
            "",
            "- The 100-point score here is still pre-rank-stability and pre-manual-directionality review.",
            "- Structural ECM or basement-membrane genes are penalized when they lack drug precedent, even if their biological localization score is high.",
            "- High druggability with systemic angiogenesis, hormone or immune-liability concerns is flagged rather than removed.",
        ]
    )
    SUMMARY.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_review(scores: pd.DataFrame) -> None:
    max_drug = float(scores["druggability_score_20"].max())
    max_final = float(scores["final_target_priority_score_100_pre_rank_stability"].max())
    failures: list[str] = []
    if len(scores) != 102:
        failures.append(f"Expected 102 candidate genes, observed {len(scores)}.")
    if max_drug > 20:
        failures.append(f"Druggability score exceeds 20: {max_drug}.")
    if max_final > 100:
        failures.append(f"Final preliminary score exceeds 100: {max_final}.")
    if scores["actionability_guardrail"].astype(str).str.len().min() == 0:
        failures.append("One or more rows lack an actionability guardrail.")
    if int(scores["druggability_penalty_points"].gt(0).sum()) == 0:
        failures.append("No penalty was applied to any gene; check whether penalties are too permissive.")

    status = "PASS" if not failures else "FAIL"
    lines = [
        "# Phase 11 self-review: druggability/actionability layer",
        "",
        f"Status: {status}",
        "",
        "## Checks",
        "",
        f"- Candidate rows: {len(scores)}",
        f"- Maximum druggability score: {max_drug:.3f} / 20",
        f"- Maximum preliminary total score: {max_final:.3f} / 100",
        f"- Genes with at least one penalty: {int(scores['druggability_penalty_points'].gt(0).sum())}",
        f"- Genes with approved DGIdb interaction: {int(scores['dgidb_approved_interaction_count'].gt(0).sum())}",
        f"- Genes with ChEMBL mechanisms: {int(scores['chembl_mechanism_count'].gt(0).sum())}",
        f"- Genes with Open Targets clinical candidates: {int(scores['ot_drug_candidate_count'].gt(0).sum())}",
        "",
        "## Scientific interpretation limits",
        "",
        "- This layer evaluates tractability, existing clinical or chemical precedent and safety caution; it does not prove therapeutic efficacy.",
        "- Directionality remains provisional unless disease-expression direction and drug mechanism align in the same disease-relevant cell state.",
        "- Open Targets isEssential, DepMap row availability and safety fields are used as screening flags, not as a complete toxicity assessment.",
        "- ChEMBL target matching prioritizes human single-protein targets with exact gene-symbol synonym matches; ambiguous target-family matches are not treated as direct evidence.",
        "",
    ]
    if failures:
        lines.extend(["## Failures", ""])
        lines.extend(f"- {failure}" for failure in failures)
    else:
        lines.extend(
            [
                "## Decision",
                "",
                "The actionability scoring layer passes automated range and guardrail checks. It is ready for rank-stability sensitivity analysis and manual target-direction review.",
            ]
        )
    REVIEW.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_source_audit() -> None:
    lines = [
        "# Actionability source API audit",
        "",
        "| Source | Endpoint | Evidence used | Notes |",
        "| --- | --- | --- | --- |",
        "| Open Targets Platform | https://api.platform.opentargets.org/api/v4/graphql | target class, tractability, clinical candidates, safety liabilities, isEssential and DepMap row availability | Queried by Ensembl gene ID. |",
        "| DGIdb | https://dgidb.org/api/graphql | drug-gene interactions, approved-drug flags, interaction direction, gene druggability categories | Queried by HGNC symbol. |",
        "| ChEMBL Web Services | https://www.ebi.ac.uk/chembl/api/data | human target matching and mechanism-of-action records | Human single-protein exact symbol matches prioritized. |",
        "",
        "This file records live API sources used to construct the cached actionability evidence. Raw API responses are stored as JSONL in this directory.",
    ]
    SOURCE_AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=0, help="Process only the first N candidates for API testing.")
    parser.add_argument("--refresh", action="store_true", help="Refresh all API caches instead of using existing JSONL files.")
    parser.add_argument("--sleep", type=float, default=0.10, help="Delay between API calls in seconds.")
    parser.add_argument("--skip-network", action="store_true", help="Use existing caches only.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    candidates = pd.read_csv(INPUT, sep="\t", keep_default_na=False)
    candidates["gene_id"] = candidates["gene_id"].astype(str).str.split(".").str[0]
    if args.limit:
        candidates = candidates.head(args.limit).copy()

    if args.skip_network:
        ot_records = read_jsonl(OT_RAW)
        dgidb_records = read_jsonl(DGIDB_RAW)
        chembl_records = read_jsonl(CHEMBL_RAW)
    else:
        if args.refresh or not OT_RAW.exists():
            ot_records = query_open_targets(candidates, args.sleep)
            write_jsonl(OT_RAW, ot_records)
        else:
            ot_records = read_jsonl(OT_RAW)
        if args.refresh or not DGIDB_RAW.exists():
            dgidb_records = query_dgidb(candidates, args.sleep)
            write_jsonl(DGIDB_RAW, dgidb_records)
        else:
            dgidb_records = read_jsonl(DGIDB_RAW)
        if args.refresh or not CHEMBL_RAW.exists():
            chembl_records = query_chembl(candidates, args.sleep)
            write_jsonl(CHEMBL_RAW, chembl_records)
        else:
            chembl_records = read_jsonl(CHEMBL_RAW)

    scores = build_scores(candidates, ot_records, dgidb_records, chembl_records)
    scores = scores.sort_values(
        ["final_target_priority_score_100_pre_rank_stability", "druggability_score_20"],
        ascending=False,
    )
    scores.to_csv(SCORES, sep="\t", index=False)
    write_source_audit()
    write_summary(scores, args)
    write_review(scores)

    print(f"Wrote {SCORES}")
    print(f"Wrote {SUMMARY}")
    print(f"Wrote {REVIEW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
