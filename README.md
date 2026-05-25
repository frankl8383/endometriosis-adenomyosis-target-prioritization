# Genetics-guided experimental-hypothesis prioritisation in endometriosis and adenomyosis

This is the minimal reproducible code repository for the manuscript
`Genetics-guided cell-contextual prioritisation of experimental hypotheses in
endometriosis and adenomyosis`.

The repository is intentionally small. It includes analysis scripts, public-data
download manifests, workflow files and environment specifications. It does not
include the manuscript, journal submission package, audit trail, raw public
datasets, local API caches, local Python package caches, rendered Word/PDF/PNG
QA artifacts, or scratch intermediates.

## What this repository can reproduce

- construction and audit of the GWAS coordinate-window candidate universe;
- candidate-level bulk-expression support layers;
- conservative single-cell/cell-label context layers for GSE179640,
  GSE203191 and Zenodo 17078290;
- cross-disease integration, druggability/actionability scoring and rank
  stability checks;
- source-data tables, shortlist tables and manuscript figures when the public
  inputs have been downloaded locally.

## What must be downloaded separately

Large public inputs must be downloaded from the cited source repositories using
`data_manifests/download_manifest.tsv` and the scripts in `scripts/`. Raw data
are excluded from GitHub by design. The largest inputs are the Zenodo 17078290
h5ad and the GSE179640 RAW tar file.

## Reproduction outline

1. Create the Python/R environment from `env/environment.yml`, or install the
   packages listed in `env/requirements-python.txt` and `env/requirements-r.txt`.
2. Download public inputs with `scripts/download_phase1_files.py` and the
   supporting download scripts listed in `workflow/reproduction_manifest.tsv`.
3. Run the commands in `workflow/reproduction_manifest.tsv` in order.
4. Compare generated outputs with the processed matrices supplied separately in
   the journal submission package.

The original analysis deliberately avoids diagnostic-model claims, genome-wide
DEG discovery as the main endpoint, causal fine-mapping claims and spatial
claims for the adenomyosis h5ad resource.

## Upload policy

Do not upload `data/raw_downloads/`, `data/raw_metadata/`, `tools/py_packages/`,
API cache JSONL files, rendered DOCX/PDF/PNG QA folders or other scratch output.
The `.gitignore` in this staging repository is configured accordingly.
