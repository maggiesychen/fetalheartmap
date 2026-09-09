#!/usr/bin/env Rscript
## Step 6b -- attach cell-state covariates and set the analysis parameters.
##
## One script for every analysis variant. What the variant controls:
##   * cis vs trans discovery pairs;
##   * which extra covariates get attached and which of them are factors;
##   * whether trans discovery pairs are restricted to protein-coding responses.
##
## Two things worth knowing before editing this (see the analysis directory's
## 6.sceptre/readmes/20260908-methods.md for the full write-up):
##
## 1. `formula_object` is passed EXPLICITLY. If it were left to sceptre's
##    auto-construction, `auto_construct_formula_object()` would silently drop
##    every continuous covariate with >= 15 distinct values -- i.e. all of
##    S_score, G2M_score, pct_ribo and pct_mito -- with no warning and no
##    indication in the printed summary.
##
## 2. This formula governs the association analyses only (calibration check,
##    power check, discovery). It does NOT reach gRNA assignment:
##    `assign_grnas()` never reads `@formula_object` and builds its own default
##    (transcriptome depth + gRNA depth + batch). To change the assignment
##    model you must pass `--grna_assignment_formula` to the Nextflow pipeline;
##    see sceptre.pipeline.grna_assignment_formula in the config.
##
## The covariate CSV is joined POSITIONALLY: row i of the CSV is sceptre cell
## index i. That is guaranteed by 05_extract_sceptre_covariates.py building it
## from the same barcodes.tsv.gz files in the same order, and is checked by
## 06c_validate_cell_ordering.py. The barcode column is compared here as a
## belt-and-braces check when the sceptre object exposes cell barcodes.

## Locate and source the shared base-R helpers next to this script.
this_script_dir_bootstrap <- function() {
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  getwd()
}
source(file.path(this_script_dir_bootstrap(), "r_utils.R"))

suppressPackageStartupMessages({
  library(sceptre)
  library(ondisc)
})

opt <- parse_cli_args(list(
  sceptre_object            = list(required = TRUE, help = "Base sceptre object .rds"),
  response_odm              = list(required = TRUE, help = "gene.odm"),
  grna_odm                  = list(required = TRUE, help = "grna.odm"),
  covariates                = list(default = "",
    help = "Extra covariate CSV from 05_extract_sceptre_covariates.py"),
  extra_covariates          = list(default = "",
    help = "Comma-separated covariates to attach (may be empty)"),
  categorical_covariates    = list(default = "",
    help = "Comma-separated subset of --extra-covariates to coerce to factor"),
  pair_type                 = list(default = "trans", help = "cis or trans"),
  distance_threshold        = list(default = 1e6, type = "numeric", help = "cis window in bp"),
  protein_coding_only       = list(flag = TRUE,
    help = "Restrict discovery responses to protein-coding genes"),
  gencode_gtf               = list(default = "",
    help = "GENCODE GTF, required with --protein-coding-only"),
  side                      = list(default = "both"),
  resampling_mechanism      = list(default = "permutations"),
  control_group             = list(default = "complement"),
  grna_integration_strategy = list(default = "union"),
  multiple_testing_alpha    = list(default = 0.1, type = "numeric"),
  output                    = list(required = TRUE, help = "Output sceptre object .rds"),
  summary                   = list(default = "",
    help = "Optional path to write the printed object summary to")
))

extra_covariates <- split_csv(opt$extra_covariates)
categorical_covariates <- split_csv(opt$categorical_covariates)
stray <- setdiff(categorical_covariates, extra_covariates)
if (length(stray) > 0) {
  stop(sprintf("--categorical-covariates not in --extra-covariates: %s",
               paste(stray, collapse = ", ")))
}

## ── 1. Load the ondisc-backed object ──────────────────────────────────────
cat("Loading sceptre object...\n")
sceptre_object <- read_ondisc_backed_sceptre_object(
  sceptre_object_fp    = opt$sceptre_object,
  response_odm_file_fp = opt$response_odm,
  grna_odm_file_fp     = opt$grna_odm
)
n_cells <- nrow(sceptre_object@covariate_data_frame)
cat("Cells:", n_cells, "\n")
cat("Built-in covariates:",
    paste(colnames(sceptre_object@covariate_data_frame), collapse = ", "), "\n")

## ── 2. Attach the extra covariates (positional) ───────────────────────────
if (length(extra_covariates) > 0) {
  if (!nzchar(opt$covariates)) {
    stop("--extra-covariates was given but --covariates (the CSV) was not")
  }
  cat("\nAttaching extra covariates from", opt$covariates, "\n")
  extra <- read.csv(opt$covariates, row.names = 1, check.names = FALSE)
  cat("  CSV dimensions:", nrow(extra), "x", ncol(extra), "\n")

  if (nrow(extra) != n_cells) {
    stop(sprintf(
      "ABORT: covariate CSV has %d rows but the sceptre object has %d cells. The positional join would be wrong.",
      nrow(extra), n_cells))
  }
  missing <- setdiff(extra_covariates, colnames(extra))
  if (length(missing) > 0) {
    stop(sprintf("covariate CSV is missing columns: %s",
                 paste(missing, collapse = ", ")))
  }

  cov_df <- sceptre_object@covariate_data_frame
  for (nm in extra_covariates) {
    values <- extra[[nm]]
    if (nm %in% categorical_covariates) values <- factor(values)
    cov_df[[nm]] <- values
    n_na <- sum(is.na(values))
    if (is.factor(values)) {
      cat(sprintf("  %-16s factor with %d levels, %d NA\n",
                  nm, length(levels(values)), n_na))
    } else {
      cat(sprintf("  %-16s numeric [%.4g, %.4g], %d NA\n",
                  nm, min(values, na.rm = TRUE), max(values, na.rm = TRUE), n_na))
    }
    if (n_na > 0) stop(sprintf("%s has %d NA values", nm, n_na))
  }
  sceptre_object@covariate_data_frame <- cov_df
} else {
  cat("\nNo extra covariates requested; using sceptre's built-in set only.\n")
}

## ── 3. Build the pair sets ────────────────────────────────────────────────
cat("\nConstructing positive-control pairs...\n")
positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
cat("  positive-control pairs:", nrow(positive_control_pairs), "\n")

if (identical(opt$pair_type, "trans")) {
  cat("Constructing trans discovery pairs...\n")
  discovery_pairs <- construct_trans_pairs(
    sceptre_object         = sceptre_object,
    positive_control_pairs = positive_control_pairs,
    pairs_to_exclude       = "pc_pairs"
  )
} else if (identical(opt$pair_type, "cis")) {
  cat("Constructing cis discovery pairs (window =", opt$distance_threshold, "bp)...\n")
  discovery_pairs <- construct_cis_pairs(
    sceptre_object         = sceptre_object,
    positive_control_pairs = positive_control_pairs,
    distance_threshold     = opt$distance_threshold
  )
} else {
  stop(sprintf("--pair-type must be cis or trans, got %s", opt$pair_type))
}
cat("  discovery pairs:", nrow(discovery_pairs), "\n")

if (opt$protein_coding_only) {
  if (!nzchar(opt$gencode_gtf)) {
    stop("--protein-coding-only requires --gencode-gtf")
  }
  suppressPackageStartupMessages(library(rtracklayer))
  cat("\nRestricting discovery responses to protein-coding genes...\n")
  cat("  GTF:", opt$gencode_gtf, "\n")
  gtf <- import(opt$gencode_gtf)
  gtf_genes <- gtf[gtf$type == "gene"]
  pc_ids <- unique(sub("\\..*", "", gtf_genes$gene_id[gtf_genes$gene_type == "protein_coding"]))
  cat("  protein-coding gene IDs in the GTF:", length(pc_ids), "\n")
  n_before <- nrow(discovery_pairs)
  discovery_pairs <- discovery_pairs[
    sub("\\..*", "", discovery_pairs$response_id) %in% pc_ids, ]
  cat("  discovery pairs:", n_before, "->", nrow(discovery_pairs), "\n")
}

## ── 4. Build the formula explicitly ───────────────────────────────────────
## Base terms are sceptre's own defaults for high-MOI data. The +1 inside the
## gRNA log terms matters: those covariates contain zeros.
base_terms <- c(
  "log(response_n_umis)", "log(response_n_nonzero)",
  "log(grna_n_umis + 1)", "log(grna_n_nonzero + 1)", "batch"
)
formula_object <- as.formula(
  paste("~", paste(c(base_terms, extra_covariates), collapse = " + "))
)
cat("\nAssociation formula:\n  ", deparse(formula_object), "\n")
cat("(This governs the calibration, power and discovery analyses only --\n")
cat(" gRNA assignment builds its own default; see the script header.)\n")

## ── 5. Set the analysis parameters ────────────────────────────────────────
cat("\nSetting analysis parameters...\n")
sceptre_object <- set_analysis_parameters(
  sceptre_object            = sceptre_object,
  discovery_pairs           = discovery_pairs,
  positive_control_pairs    = positive_control_pairs,
  formula_object            = formula_object,
  side                      = opt$side,
  resampling_mechanism      = opt$resampling_mechanism,
  control_group             = opt$control_group,
  grna_integration_strategy = opt$grna_integration_strategy,
  multiple_testing_alpha    = opt$multiple_testing_alpha
)

print(sceptre_object)

## Confirm the formula actually survived, rather than trusting the call.
stored <- paste(deparse(sceptre_object@formula_object), collapse = "")
for (nm in extra_covariates) {
  if (!grepl(nm, stored, fixed = TRUE)) {
    stop(sprintf("covariate %s is absent from the stored formula: %s", nm, stored))
  }
}
cat("\nVerified: all", length(extra_covariates),
    "extra covariates are present in @formula_object.\n")

saveRDS(sceptre_object, opt$output)
cat("Wrote", opt$output, "\n")

if (nzchar(opt$summary)) {
  con <- file(opt$summary, open = "wt")
  sink(con); print(sceptre_object); sink(); close(con)
  cat("Wrote", opt$summary, "\n")
}
