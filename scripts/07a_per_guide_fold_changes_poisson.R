#!/usr/bin/env Rscript
## Step 7a -- per-gRNA Poisson log2 fold change.
##
## SCEPTRE reports one fold change per *target*, pooling that target's gRNAs
## (grna_integration_strategy = "union"). This recomputes the same estimator one
## gRNA at a time, which is what a per-guide knockdown figure needs.
##
## The estimator is sceptre's own (sceptre book §10.4, Eq. 10.15):
##
##   null model      Y_j ~ Pois(mu_j),  log(mu_j) = tau^T Z_j
##   fold change     xi_hat = log2( sum_{j in trt} Y_j / sum_{j in trt} mu_hat_j )
##
## The null GLM is fitted once per response gene across all cells and its fitted
## values cached, then reused for every gRNA -- which is what makes this
## tractable over ~18k gRNA-gene pairs.
##
## Two guide sets, selected with --guide-set:
##   targeting      every (target, response) pair in the sceptre object, i.e.
##                  positive controls plus discovery pairs
##   non-targeting  the non-targeting gRNAs against a supplied gene list, which
##                  gives the null distribution the targeting guides are
##                  compared against
##
## Cell ordering: the gRNA assignment matrix's columns are in sceptre cell
## order, which is the same order as the covariate data frame's rows. That is
## asserted below rather than assumed. The original script additionally read a
## `grna_binary_cols.csv` to re-attach column names; that is unnecessary because
## the matrix already carries them, and the estimator only ever uses column
## *positions*.

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
  library(data.table)
  library(Matrix)
})

opt <- parse_cli_args(list(
  sceptre_object      = list(required = TRUE, help = "sceptre object .rds"),
  response_odm        = list(required = TRUE, help = "gene.odm"),
  grna_odm            = list(required = TRUE, help = "grna.odm"),
  assignment_matrix   = list(required = TRUE,
    help = "grna_assignment_matrix.rds from the SCEPTRE pipeline"),
  guide_targets       = list(required = TRUE, help = "TSV with grna_id / grna_target"),
  ens2symbol          = list(required = TRUE, help = "TSV ensembl_id / gene_symbol"),
  guide_set           = list(default = "targeting", help = "targeting or non-targeting"),
  gene_list           = list(default = "",
    help = "Gene-symbol table restricting the responses (required for non-targeting)"),
  gene_list_column    = list(default = "Gene Name",
    help = "Column of --gene-list holding gene symbols"),
  nt_target_label     = list(default = "non-targeting",
    help = "grna_target value marking non-targeting gRNAs"),
  covariate_formula   = list(default = "~ batch + log(response_n_umis)",
    help = "Null-model covariates"),
  min_cells_per_guide = list(default = 5L, type = "integer",
    help = "Minimum treatment cells to report a fold change"),
  output              = list(required = TRUE, help = "Output TSV")
))

if (!opt$guide_set %in% c("targeting", "non-targeting")) {
  stop("--guide-set must be 'targeting' or 'non-targeting'")
}
if (identical(opt$guide_set, "non-targeting") && !nzchar(opt$gene_list)) {
  stop("--guide-set non-targeting requires --gene-list")
}

## ── Load ──────────────────────────────────────────────────────────────────
cat("Loading sceptre object...\n")
sceptre_obj <- read_ondisc_backed_sceptre_object(
  sceptre_object_fp    = opt$sceptre_object,
  response_odm_file_fp = opt$response_odm,
  grna_odm_file_fp     = opt$grna_odm
)

cov_df  <- sceptre_obj@covariate_data_frame
n_cells <- nrow(cov_df)
cat("Cells:", n_cells, "\n")

cat("Loading gRNA assignment matrix...\n")
guide_mat <- readRDS(opt$assignment_matrix)
cat("Assignment matrix:", nrow(guide_mat), "gRNAs x", ncol(guide_mat), "cells\n")

## The estimator indexes cells positionally, so a column-count mismatch would
## silently attribute counts to the wrong cells.
if (ncol(guide_mat) != n_cells) {
  stop(sprintf(
    "ABORT: assignment matrix has %d columns but the sceptre object has %d cells.",
    ncol(guide_mat), n_cells))
}
cat("Verified: assignment matrix column count matches the sceptre cell count\n")

## ── Null-model design matrix ──────────────────────────────────────────────
form <- as.formula(opt$covariate_formula)
needed <- all.vars(form)
missing <- setdiff(needed, colnames(cov_df))
if (length(missing) > 0) {
  stop(sprintf("covariate_data_frame lacks: %s (available: %s)",
               paste(missing, collapse = ", "),
               paste(colnames(cov_df), collapse = ", ")))
}
covariate_mat <- model.matrix(form, data = cov_df)
cat("Design matrix columns:", ncol(covariate_mat),
    "(", opt$covariate_formula, ")\n")

gene_odm <- sceptre_obj@response_matrix[[1]]
gene_ids <- gene_odm@feature_ids
cat("Response ODM:", length(gene_ids), "genes\n")

ens2sym <- fread(opt$ens2symbol)
sym_map <- setNames(ens2sym$gene_symbol, ens2sym$ensembl_id)
strip_version <- function(x) sub("\\..*", "", x)

meta <- fread(opt$guide_targets)
guides_in_mat <- rownames(guide_mat)
meta <- meta[meta$grna_id %in% guides_in_mat, ]
target_to_guides <- split(meta$grna_id, meta$grna_target)

## ── Build the (target, response) work list ────────────────────────────────
if (identical(opt$guide_set, "targeting")) {
  pos_pairs  <- sceptre_obj@positive_control_pairs
  disc_pairs <- sceptre_obj@discovery_pairs
  all_pairs <- unique(rbind(
    data.frame(grna_target = as.character(pos_pairs$grna_target),
               response_id = as.character(pos_pairs$response_id),
               is_pos_ctrl = TRUE, stringsAsFactors = FALSE),
    data.frame(grna_target = as.character(disc_pairs$grna_target),
               response_id = as.character(disc_pairs$response_id),
               is_pos_ctrl = FALSE, stringsAsFactors = FALSE)
  ))
  all_pairs <- all_pairs[all_pairs$grna_target != opt$nt_target_label, ]
  cat("Unique targeting pairs:", nrow(all_pairs),
      "across", length(unique(all_pairs$grna_target)), "targets\n")
} else {
  genes <- read.table(opt$gene_list, header = TRUE, sep = "\t",
                      quote = "\"", comment.char = "", check.names = FALSE)
  if (!opt$gene_list_column %in% colnames(genes)) {
    stop(sprintf("--gene-list-column '%s' not in %s (columns: %s)",
                 opt$gene_list_column, opt$gene_list,
                 paste(colnames(genes), collapse = ", ")))
  }
  symbols <- unique(trimws(as.character(genes[[opt$gene_list_column]])))
  cat("Gene list:", length(symbols), "symbols\n")

  ## symbol -> Ensembl, then match against the ODM's versioned IDs
  sym2ens <- setNames(ens2sym$ensembl_id, ens2sym$gene_symbol)
  wanted_ens <- unname(sym2ens[symbols])
  odm_stripped <- strip_version(gene_ids)
  response_ids <- gene_ids[odm_stripped %in% wanted_ens[!is.na(wanted_ens)]]
  cat("Matched to", length(response_ids), "responses in the ODM\n")
  unmatched <- symbols[is.na(sym2ens[symbols])]
  if (length(unmatched) > 0) {
    cat("  WARNING: no Ensembl ID for:", paste(unmatched, collapse = ", "), "\n")
  }

  all_pairs <- data.frame(
    grna_target = opt$nt_target_label,
    response_id = response_ids,
    is_pos_ctrl = FALSE,
    stringsAsFactors = FALSE
  )
  cat("Non-targeting pairs:", nrow(all_pairs), "\n")
}

## Drop anything not representable.
missing_resp <- setdiff(all_pairs$response_id, gene_ids)
if (length(missing_resp) > 0) {
  cat(sprintf("WARNING: dropping %d response_ids absent from the ODM\n",
              length(missing_resp)))
  all_pairs <- all_pairs[!all_pairs$response_id %in% missing_resp, ]
}
missing_tgt <- setdiff(all_pairs$grna_target, names(target_to_guides))
if (length(missing_tgt) > 0) {
  cat(sprintf("WARNING: dropping %d targets absent from the assignment matrix\n",
              length(missing_tgt)))
  all_pairs <- all_pairs[!all_pairs$grna_target %in% missing_tgt, ]
}
if (nrow(all_pairs) == 0L) stop("no pairs left to test")
cat(sprintf("Proceeding with %d pairs over %d response genes\n",
            nrow(all_pairs), length(unique(all_pairs$response_id))))

## ── Fit the null GLM once per response gene ───────────────────────────────
unique_responses <- unique(all_pairs$response_id)
cat(sprintf("\nFitting null Poisson GLMs for %d response genes...\n",
            length(unique_responses)))

fitted_cache <- list()
for (i in seq_along(unique_responses)) {
  resp_id  <- unique_responses[i]
  gene_idx <- which(gene_ids == resp_id)
  y <- as.numeric(gene_odm[gene_idx, ])

  if (sum(y) == 0) {
    fitted_cache[[resp_id]] <- rep(0, n_cells)
  } else {
    fit <- tryCatch(glm.fit(x = covariate_mat, y = y, family = poisson()),
                    error = function(e) NULL)
    ## Fall back to the marginal mean if the GLM will not converge, so one
    ## pathological gene cannot abort the whole run.
    fitted_cache[[resp_id]] <-
      if (is.null(fit)) rep(mean(y), n_cells) else fit$fitted.values
  }
  if (i %% 50 == 0) cat(sprintf("  fitted %d / %d\n", i, length(unique_responses)))
}
cat("Done fitting.\n")

## ── Per-gRNA fold changes ─────────────────────────────────────────────────
min_cells <- opt$min_cells_per_guide
targets <- intersect(unique(all_pairs$grna_target), names(target_to_guides))
results_list <- vector("list", nrow(all_pairs) * 10L)
k <- 0L

for (tgt in targets) {
  tgt_guides <- intersect(target_to_guides[[tgt]], guides_in_mat)
  if (length(tgt_guides) == 0L) next
  tgt_pairs <- all_pairs[all_pairs$grna_target == tgt, , drop = FALSE]

  for (i in seq_len(nrow(tgt_pairs))) {
    resp_id <- tgt_pairs$response_id[i]
    gene_idx <- which(gene_ids == resp_id)
    if (length(gene_idx) == 0L) next
    y      <- as.numeric(gene_odm[gene_idx, ])
    mu_hat <- fitted_cache[[resp_id]]

    for (g in tgt_guides) {
      trt_idx <- which(as.logical(guide_mat[g, ]))
      n_trt <- length(trt_idx)
      if (n_trt >= min_cells) {
        sum_y  <- sum(y[trt_idx])
        sum_mu <- sum(mu_hat[trt_idx])
        lfc <- if (sum_y > 0 && sum_mu > 0) log2(sum_y / sum_mu) else NA_real_
      } else {
        lfc <- NA_real_
      }
      k <- k + 1L
      results_list[[k]] <- list(
        grna_id          = g,
        grna_target      = tgt,
        response_id      = resp_id,
        gene_symbol      = unname(sym_map[strip_version(resp_id)]),
        is_pos_ctrl      = tgt_pairs$is_pos_ctrl[i],
        n_trt            = n_trt,
        log2_fold_change = round(lfc, 6L)
      )
    }
  }
  cat(sprintf("  %s: %d gRNAs x %d genes\n",
              tgt, length(tgt_guides), nrow(tgt_pairs)))
}

results <- rbindlist(results_list[seq_len(k)])
cat(sprintf("\nComputed fold changes for %d gRNA-gene pairs (%d non-NA)\n",
            nrow(results), sum(!is.na(results$log2_fold_change))))

setorder(results, grna_target, response_id, grna_id)
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
fwrite(results, opt$output, sep = "\t")
cat("Wrote", opt$output, "\n")
