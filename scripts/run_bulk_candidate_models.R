#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

project_root <- normalizePath(getwd())
raw_dir <- file.path(project_root, "data", "raw_downloads")
results_dir <- file.path(project_root, "results", "bulk")
model_dir <- file.path(results_dir, "models")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

clean_na <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "nan", "Unknown", "unknown")] <- NA_character_
  x
}

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

aggregate_counts <- function(counts, gene_ids) {
  gene_ids <- sub("\\.\\d+$", "", as.character(gene_ids))
  keep <- !is.na(gene_ids) & gene_ids != ""
  counts <- counts[keep, , drop = FALSE]
  gene_ids <- gene_ids[keep]
  rowsum(as.matrix(counts), group = gene_ids, reorder = FALSE)
}

merge_candidate_results <- function(candidates, result, by_field = "gene_id") {
  if (nrow(result) == 0) {
    candidates$analysis_status <- "not_tested"
    return(candidates)
  }
  merged <- merge(candidates, result, by = by_field, all.x = TRUE, sort = FALSE)
  merged$analysis_status <- ifelse(is.na(merged$P.Value), "not_tested_or_filtered", "tested")
  merged
}

write_result <- function(df, path) {
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

bind_summary_frames <- function(frames) {
  frames <- lapply(frames, as.data.frame)
  all_names <- unique(unlist(lapply(frames, names)))
  frames <- lapply(frames, function(frame) {
    missing <- setdiff(all_names, names(frame))
    for (name in missing) {
      frame[[name]] <- NA
    }
    frame[, all_names, drop = FALSE]
  })
  do.call(rbind, frames)
}

candidates <- read_tsv(file.path(results_dir, "candidate_genes_unique.tsv"))
candidate_gene_ids <- unique(candidates$gene_id)
candidate_symbols <- unique(candidates$gene_symbol)

summaries <- list()

# GSE234354: menstrual-cycle dependence in endometrial reference RNA-seq.
run_gse234354 <- function() {
  counts <- read.delim(
    gzfile(file.path(raw_dir, "GSE234354__GSE234354_gene_count_matrix.txt.gz")),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  gene_ids <- counts$gene_id
  count_mat <- aggregate_counts(counts[, setdiff(names(counts), "gene_id"), drop = FALSE], gene_ids)
  meta <- read_tsv(file.path(results_dir, "metadata_GSE234354.tsv"))
  samples <- intersect(colnames(count_mat), meta$matrix_sample_id)
  meta <- meta[match(samples, meta$matrix_sample_id), , drop = FALSE]
  count_mat <- count_mat[, samples, drop = FALSE]
  meta$cycle_stage <- clean_na(meta$cycle_stage)
  meta$batch <- clean_na(meta$batch)
  usable <- !is.na(meta$cycle_stage) & !is.na(meta$batch)
  meta <- meta[usable, , drop = FALSE]
  count_mat <- count_mat[, meta$matrix_sample_id, drop = FALSE]
  meta$cycle_stage <- factor(paste0("stage_", meta$cycle_stage))
  meta$batch <- factor(make.names(meta$batch))

  design <- model.matrix(~ cycle_stage + batch, data = meta)
  y <- DGEList(counts = count_mat)
  y <- calcNormFactors(y)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  v <- voom(y, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  cycle_coefs <- grep("^cycle_stage", colnames(design))
  tab <- topTable(fit, coef = cycle_coefs, number = Inf, sort.by = "none")
  tab$gene_id <- rownames(tab)
  tab <- tab[, c("gene_id", setdiff(names(tab), "gene_id"))]
  out <- merge_candidate_results(candidates, tab, "gene_id")
  out$analysis <- "GSE234354_cycle_stage_F_test"
  write_result(out, file.path(model_dir, "GSE234354_cycle_candidate_results.tsv"))
  data.frame(
    analysis = "GSE234354_cycle_stage_F_test",
    samples_used = ncol(count_mat),
    genes_tested_after_filter = nrow(tab),
    candidates_tested = sum(out$analysis_status == "tested"),
    nominal_p_lt_0_05 = sum(out$P.Value < 0.05, na.rm = TRUE),
    fdr_lt_0_10 = sum(out$adj.P.Val < 0.10, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# GSE313775: disease contrasts within sorted Th subsets.
run_gse313775 <- function() {
  counts <- read.delim(
    gzfile(file.path(raw_dir, "GSE313775__GSE313775_rawCountMatrix.tsv.gz")),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  gene_ids <- counts$Gene
  sample_cols <- setdiff(names(counts), c("Gene", "Symbol"))
  symbol_map <- counts[, c("Gene", "Symbol")]
  symbol_map$Gene <- sub("\\.\\d+$", "", symbol_map$Gene)
  symbol_map <- symbol_map[!duplicated(symbol_map$Gene), , drop = FALSE]
  count_mat <- aggregate_counts(counts[, sample_cols, drop = FALSE], gene_ids)
  meta <- read_tsv(file.path(results_dir, "metadata_GSE313775.tsv"))
  samples <- intersect(colnames(count_mat), meta$matrix_sample_id)
  meta <- meta[match(samples, meta$matrix_sample_id), , drop = FALSE]
  count_mat <- count_mat[, samples, drop = FALSE]
  meta$cell_subset_clean <- gsub("-", "_", meta$cell_subset)
  meta$group <- factor(make.names(paste(meta$disease_status, meta$cell_subset_clean, sep = ".")))
  design <- model.matrix(~ 0 + group, data = meta)
  colnames(design) <- sub("^group", "", colnames(design))

  y <- DGEList(counts = count_mat)
  y <- calcNormFactors(y)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  v <- voom(y, design, plot = FALSE)
  corfit <- duplicateCorrelation(v, design, block = meta$donor_id)
  fit <- lmFit(v, design, block = meta$donor_id, correlation = corfit$consensus)
  contrast_matrix <- makeContrasts(
    endometriosis_vs_control_Th1 = endometriosis.Th1 - control.Th1,
    endometriosis_vs_control_Th1_17 = endometriosis.Th1_17 - control.Th1_17,
    endometriosis_vs_control_Th17 = endometriosis.Th17 - control.Th17,
    Th17_interaction_vs_Th1 = (endometriosis.Th17 - control.Th17) - (endometriosis.Th1 - control.Th1),
    Th1_17_interaction_vs_Th1 = (endometriosis.Th1_17 - control.Th1_17) - (endometriosis.Th1 - control.Th1),
    levels = design
  )
  fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))
  all_out <- list()
  summary_rows <- list()
  for (contrast in colnames(contrast_matrix)) {
    tab <- topTable(fit2, coef = contrast, number = Inf, sort.by = "none")
    tab$gene_id <- rownames(tab)
    tab <- merge(tab, symbol_map, by.x = "gene_id", by.y = "Gene", all.x = TRUE, sort = FALSE)
    names(tab)[names(tab) == "Symbol"] <- "matrix_symbol"
    tab <- tab[, c("gene_id", "matrix_symbol", setdiff(names(tab), c("gene_id", "matrix_symbol")))]
    out <- merge_candidate_results(candidates, tab, "gene_id")
    out$analysis <- "GSE313775_Th_subset_disease_contrast"
    out$contrast <- contrast
    all_out[[contrast]] <- out
    summary_rows[[contrast]] <- data.frame(
      analysis = paste0("GSE313775_", contrast),
      samples_used = ncol(count_mat),
      genes_tested_after_filter = nrow(tab),
      candidates_tested = sum(out$analysis_status == "tested"),
      nominal_p_lt_0_05 = sum(out$P.Value < 0.05, na.rm = TRUE),
      fdr_lt_0_10 = sum(out$adj.P.Val < 0.10, na.rm = TRUE),
      duplicate_correlation = corfit$consensus,
      stringsAsFactors = FALSE
    )
  }
  combined <- do.call(rbind, all_out)
  write_result(combined, file.path(model_dir, "GSE313775_th_subset_candidate_results.tsv"))
  do.call(rbind, summary_rows)
}

run_array_limma <- function(expr, meta, sample_col, group, block = NULL, contrasts, by_field, analysis_name, output_name) {
  sample_cols <- intersect(names(expr), meta[[sample_col]])
  meta <- meta[match(sample_cols, meta[[sample_col]]), , drop = FALSE]
  expr_mat <- as.matrix(expr[, sample_cols, drop = FALSE])
  rownames(expr_mat) <- expr[[by_field]]
  meta$model_group <- factor(make.names(group[match(sample_cols, meta[[sample_col]])]))
  design <- model.matrix(~ 0 + model_group, data = meta)
  colnames(design) <- sub("^model_group", "", colnames(design))
  if (!is.null(block)) {
    block_vec <- meta[[block]]
    corfit <- duplicateCorrelation(expr_mat, design, block = block_vec)
    fit <- lmFit(expr_mat, design, block = block_vec, correlation = corfit$consensus)
    duplicate_correlation <- corfit$consensus
  } else {
    fit <- lmFit(expr_mat, design)
    duplicate_correlation <- NA_real_
  }
  contrast_matrix <- makeContrasts(contrasts = contrasts, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))
  all_out <- list()
  summary_rows <- list()
  for (contrast in colnames(contrast_matrix)) {
    tab <- topTable(fit2, coef = contrast, number = Inf, sort.by = "none")
    tab[[by_field]] <- rownames(tab)
    tab <- tab[, c(by_field, setdiff(names(tab), by_field))]
    if (by_field == "gene_symbol") {
      base <- candidates[, c("gene_symbol", "gene_id", "genetic_priority", "neighborhoods", "ld_neighborhood_class", "module_hint_preliminary", "gene_biotype", "gene_description")]
      out <- merge(base, tab, by = "gene_symbol", all.x = TRUE, sort = FALSE)
      out$analysis_status <- ifelse(is.na(out$P.Value), "not_tested_or_filtered", "tested")
    } else {
      out <- merge_candidate_results(candidates, tab, by_field)
    }
    out$analysis <- analysis_name
    out$contrast <- contrast
    all_out[[contrast]] <- out
    summary_rows[[contrast]] <- data.frame(
      analysis = paste0(analysis_name, "_", contrast),
      samples_used = ncol(expr_mat),
      genes_tested_after_filter = nrow(tab),
      candidates_tested = sum(out$analysis_status == "tested"),
      nominal_p_lt_0_05 = sum(out$P.Value < 0.05, na.rm = TRUE),
      fdr_lt_0_10 = sum(out$adj.P.Val < 0.10, na.rm = TRUE),
      duplicate_correlation = duplicate_correlation,
      stringsAsFactors = FALSE
    )
  }
  write_result(do.call(rbind, all_out), file.path(model_dir, output_name))
  do.call(rbind, summary_rows)
}

# GSE141549: lesion/endometrium/peritoneum contrasts.
run_gse141549 <- function() {
  expr <- read_tsv(file.path(results_dir, "GSE141549_candidate_gene_expression_median_probe.tsv"))
  names(expr)[names(expr) == "Gene_symbol"] <- "gene_symbol"
  meta <- read_tsv(file.path(results_dir, "metadata_GSE141549.tsv"))
  meta <- meta[meta$is_replicate_label != "True", , drop = FALSE]
  group <- ifelse(meta$broad_tissue_class == "lesion", "lesion", meta$tissue_subtype_preliminary)
  keep <- group %in% c("lesion", "control_endometrium", "patient_eutopic_endometrium", "control_peritoneum", "patient_peritoneum")
  meta <- meta[keep, , drop = FALSE]
  group <- group[keep]
  contrasts <- c(
    "lesion-control_endometrium",
    "lesion-patient_eutopic_endometrium",
    "patient_eutopic_endometrium-control_endometrium",
    "patient_peritoneum-control_peritoneum",
    "lesion-patient_peritoneum"
  )
  run_array_limma(
    expr = expr,
    meta = meta,
    sample_col = "matrix_sample_id",
    group = group,
    block = "sample_link_id",
    contrasts = contrasts,
    by_field = "gene_symbol",
    analysis_name = "GSE141549_tissue_contrast",
    output_name = "GSE141549_tissue_candidate_results.tsv"
  )
}

# GSE51981: independent endometrium validation, adjusted for cycle phase.
run_gse51981 <- function() {
  meta <- read_tsv(file.path(results_dir, "metadata_GSE51981.tsv"))
  meta$endo_status <- clean_na(meta$endometriosis_no_endometriosis)
  meta$cycle_phase <- clean_na(gsub(" Endometrial tissue", "", meta$tissue, fixed = TRUE))
  keep <- meta$endo_status %in% c("Endometriosis", "Non-Endometriosis") & !is.na(meta$cycle_phase)
  meta <- meta[keep, , drop = FALSE]
  meta$endo_status <- factor(meta$endo_status, levels = c("Non-Endometriosis", "Endometriosis"))
  meta$cycle_phase <- factor(make.names(meta$cycle_phase))

  fit_expr <- function(expr, output_name, analysis_label) {
    names(expr)[names(expr) == "Gene_symbol"] <- "gene_symbol"
    sample_cols <- intersect(names(expr), meta$matrix_sample_id)
    meta_fit <- meta[match(sample_cols, meta$matrix_sample_id), , drop = FALSE]
    expr_mat <- as.matrix(expr[, sample_cols, drop = FALSE])
    rownames(expr_mat) <- expr$gene_symbol
    design <- model.matrix(~ endo_status + cycle_phase, data = meta_fit)
    fit <- eBayes(lmFit(expr_mat, design))
    coef_name <- "endo_statusEndometriosis"
    tab <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
    tab$gene_symbol <- rownames(tab)
    tab <- tab[, c("gene_symbol", setdiff(names(tab), "gene_symbol"))]
    base <- candidates[, c("gene_symbol", "gene_id", "genetic_priority", "neighborhoods", "ld_neighborhood_class", "module_hint_preliminary", "gene_biotype", "gene_description")]
    out <- merge(base, tab, by = "gene_symbol", all.x = TRUE, sort = FALSE)
    out$analysis_status <- ifelse(is.na(out$P.Value), "not_tested_or_filtered", "tested")
    out$analysis <- analysis_label
    out$contrast <- "Endometriosis_vs_Non_Endometriosis_adjusted_cycle"
    write_result(out, file.path(model_dir, output_name))
    data.frame(
      analysis = analysis_label,
      samples_used = ncol(expr_mat),
      genes_tested_after_filter = nrow(tab),
      candidates_tested = sum(out$analysis_status == "tested"),
      nominal_p_lt_0_05 = sum(out$P.Value < 0.05, na.rm = TRUE),
      fdr_lt_0_10 = sum(out$adj.P.Val < 0.10, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  expr_all <- read_tsv(file.path(results_dir, "GSE51981_candidate_gene_expression_median_probe.tsv"))
  all_summary <- fit_expr(
    expr_all,
    "GSE51981_endometrium_candidate_results.tsv",
    "GSE51981_Endometriosis_vs_Non_Endometriosis_adjusted_cycle_all_mapped_probes"
  )

  probe <- read_tsv(file.path(results_dir, "GSE51981_candidate_probe_expression.tsv"))
  sample_cols <- grep("^GSM", names(probe), value = TRUE)
  probe_single <- probe[probe$probe_mapping_status == "single_gene_probe", , drop = FALSE]
  expr_single <- aggregate(
    probe_single[, sample_cols, drop = FALSE],
    by = list(gene_symbol = probe_single$gene_symbol, candidate_gene_id = probe_single$candidate_gene_id),
    FUN = median
  )
  single_summary <- fit_expr(
    expr_single,
    "GSE51981_endometrium_candidate_results_single_gene_probe_sensitivity.tsv",
    "GSE51981_Endometriosis_vs_Non_Endometriosis_adjusted_cycle_single_gene_probes"
  )
  rbind(all_summary, single_summary)
}

summaries[["GSE234354"]] <- run_gse234354()
summaries[["GSE313775"]] <- run_gse313775()
summaries[["GSE141549"]] <- run_gse141549()
summaries[["GSE51981"]] <- run_gse51981()

summary_df <- bind_summary_frames(summaries)
write_result(summary_df, file.path(model_dir, "bulk_candidate_model_summary.tsv"))

summary_lines <- c(
  "# Bulk candidate-gene association model summary",
  "",
  "These models test the genetics-prioritized candidate gene set only. They are not genome-wide DEG screens.",
  "",
  paste(capture.output(print(summary_df, row.names = FALSE)), collapse = "\n"),
  "",
  "Model notes:",
  "",
  "- GSE234354 uses edgeR TMM normalization and limma-voom; the F-test captures menstrual-cycle stage dependence adjusted for batch.",
  "- GSE313775 uses edgeR TMM, limma-voom and duplicateCorrelation by donor for disease contrasts within Th subsets and subset-specific interactions.",
  "- GSE141549 uses limma on median-collapsed candidate gene expression, excludes columns marked as replicate labels for the primary contrast, and blocks by sample-link ID.",
  "- GSE51981 uses limma on median-collapsed GPL570-annotated gene expression, testing endometriosis versus non-endometriosis adjusted for cycle phase; single-gene probe sensitivity is also written.",
  ""
)
writeLines(summary_lines, file.path(model_dir, "bulk_candidate_model_summary.md"))

self_review <- c(
  "# Phase 4 bulk candidate-model self-review",
  "",
  "Verdict: PASS_WITH_CONDITIONS",
  "",
  "Checks passed:",
  "",
  "- RNA-seq count datasets were modeled with edgeR normalization and limma-voom rather than raw count t-tests.",
  "- GSE313775 repeated Th-subset measurements were modeled with donor blocking.",
  "- Microarray datasets were analyzed as normalized continuous expression with limma.",
  "- Model outputs are limited to the GWAS-prioritized candidate universe.",
  "",
  "Conditions before interpretation:",
  "",
  "- GSE141549 tissue-code definitions remain preliminary until explicitly tied to the Scientific Data descriptor.",
  "- GSE313775 donor-to-GSM mapping depends on matrix donor order matching GEO biological replicate order.",
  "- GSE51981 single-gene probe sensitivity has been run and should be reported alongside the all-mapped-probe analysis.",
  "- These candidate models provide expression support, not causal target evidence; they must be integrated with single-cell/spatial and druggability layers.",
  ""
)
writeLines(self_review, file.path(model_dir, "phase4_bulk_candidate_model_self_review.md"))

message("Wrote bulk candidate model outputs to ", model_dir)
