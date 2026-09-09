#!/usr/bin/env Rscript
## Step 7b -- per-gRNA pseudobulk expression of the target genes.
##
## A model-free counterpart to the Poisson fold change in 07a. For each gene in
## the panel it reports CP10K pseudobulk expression of that gene in the cells
## carrying each of its TSS gRNAs, in the cells carrying each non-targeting
## gRNA, and in a per-gene non-targeting *pool*.
##
## Having both estimators matters: they share the assignment matrix and the
## count data but nothing else -- no GLM, no covariates beyond library size --
## so agreement between them is evidence the knockdown is real rather than a
## modelling artefact. `07d_plot_nt_vs_tss.py` reports the correlation.
##
## The per-gene NT pool is NT-carrying cells that carry *none* of that gene's
## TSS gRNAs. Under high MOI an NT cell usually also carries targeting gRNAs;
## those are overwhelmingly for other loci, and any for the gene being measured
## are excluded here so the reference is never contaminated by the very
## perturbation it is the reference for.
##
## Two CP10K flavours are reported and they are not the same statistic:
##   cp10k            pseudobulk -- pooled counts over pooled library size
##   mean_cell_cp10k  mean of per-cell CP10K
## The published figure uses `cp10k`.

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
  sceptre_object      = list(required = TRUE),
  response_odm        = list(required = TRUE),
  grna_odm            = list(required = TRUE),
  assignment_matrix   = list(required = TRUE),
  guide_targets       = list(required = TRUE),
  ens2symbol          = list(required = TRUE),
  gene_list           = list(required = TRUE,
    help = "Gene-symbol table naming the target panel"),
  gene_list_column    = list(default = "Gene Name"),
  nt_target_label     = list(default = "non-targeting"),
  tss_target_regex    = list(default = "^ENSG"),
  cp10k_scale         = list(default = 1e4, type = "numeric"),
  min_cells_per_guide = list(default = 5L, type = "integer"),
  output_per_guide    = list(required = TRUE),
  output_nt_pool      = list(required = TRUE)
))

SCALE <- opt$cp10k_scale
MIN_CELLS <- opt$min_cells_per_guide
strip_version <- function(x) sub("\\..*", "", x)

## ── Panel genes ───────────────────────────────────────────────────────────
genes <- read.table(opt$gene_list, header = TRUE, sep = "\t",
                    quote = "\"", comment.char = "", check.names = FALSE)
if (!opt$gene_list_column %in% colnames(genes)) {
  stop(sprintf("--gene-list-column '%s' not found in %s",
               opt$gene_list_column, opt$gene_list))
}
panel_syms <- unique(trimws(as.character(genes[[opt$gene_list_column]])))
cat("Panel genes:", length(panel_syms), "\n")

ens2sym <- fread(opt$ens2symbol)
sym_map <- setNames(ens2sym$gene_symbol, ens2sym$ensembl_id)

meta <- fread(opt$guide_targets)
tss_targets <- unique(meta$grna_target[grepl(opt$tss_target_regex, meta$grna_target)])
target_syms <- sym_map[strip_version(tss_targets)]
panel_targets <- tss_targets[!is.na(target_syms) & target_syms %in% panel_syms]
cat("Panel genes with a TSS target in the library:", length(panel_targets), "\n")

## ── Load ──────────────────────────────────────────────────────────────────
cat("Loading sceptre object...\n")
sceptre_obj <- read_ondisc_backed_sceptre_object(
  sceptre_object_fp    = opt$sceptre_object,
  response_odm_file_fp = opt$response_odm,
  grna_odm_file_fp     = opt$grna_odm
)

guide_mat <- readRDS(opt$assignment_matrix)
cov_df  <- sceptre_obj@covariate_data_frame
n_umis  <- as.numeric(cov_df$response_n_umis)
n_cells <- length(n_umis)
cat("Assignment matrix:", nrow(guide_mat), "gRNAs x", ncol(guide_mat), "cells\n")
if (ncol(guide_mat) != n_cells) {
  stop(sprintf("ABORT: assignment matrix has %d columns, sceptre object has %d cells",
               ncol(guide_mat), n_cells))
}

gene_odm <- sceptre_obj@response_matrix[[1]]
gene_ids <- gene_odm@feature_ids

guides_in_mat <- rownames(guide_mat)
ntc_guides <- intersect(
  meta$grna_id[meta$grna_target == opt$nt_target_label], guides_in_mat)
cat("Non-targeting gRNAs in the assignment matrix:", length(ntc_guides), "\n")
if (length(ntc_guides) == 0L) {
  stop(sprintf("no gRNAs with grna_target == '%s'", opt$nt_target_label))
}

nt_any <- colSums(guide_mat[ntc_guides, , drop = FALSE]) > 0L
cat("Cells carrying at least one non-targeting gRNA:", sum(nt_any), "\n")

panel_targets <- intersect(panel_targets, gene_ids)
cat("Panel responses present in the ODM:", length(panel_targets), "\n\n")
if (length(panel_targets) == 0L) stop("no panel genes found in the response ODM")

## ── Summarise ─────────────────────────────────────────────────────────────
summarise_cells <- function(idx, y) {
  list(n_cells         = length(idx),
       sum_counts      = sum(y[idx]),
       sum_umis        = sum(n_umis[idx]),
       cp10k           = SCALE * sum(y[idx]) / sum(n_umis[idx]),
       mean_cell_cp10k = mean(SCALE * y[idx] / n_umis[idx]),
       frac_nonzero    = mean(y[idx] > 0))
}

per_guide <- vector("list", length(panel_targets) * (length(ntc_guides) + 40L))
pool_rows <- vector("list", length(panel_targets))
k <- 0L

for (i in seq_along(panel_targets)) {
  resp_id <- panel_targets[i]
  sym <- unname(sym_map[strip_version(resp_id)])
  y <- as.numeric(gene_odm[which(gene_ids == resp_id), ])

  tss_guides <- intersect(meta$grna_id[meta$grna_target == resp_id], guides_in_mat)

  if (length(tss_guides)) {
    tss_any <- colSums(guide_mat[tss_guides, , drop = FALSE]) > 0L
  } else {
    tss_any <- rep(FALSE, n_cells)
  }
  nt_idx <- which(nt_any & !tss_any)
  pool <- summarise_cells(nt_idx, y)
  pool_rows[[i]] <- data.table(response_id = resp_id, gene_symbol = sym,
                               as.data.table(pool))

  for (grp in list(list(ids = tss_guides, series = "TSS"),
                   list(ids = ntc_guides, series = "NTC"))) {
    for (g in grp$ids) {
      idx <- which(as.logical(guide_mat[g, ]))
      if (length(idx) < MIN_CELLS) next
      s <- summarise_cells(idx, y)
      k <- k + 1L
      per_guide[[k]] <- data.table(
        grna_id     = g,
        series      = grp$series,
        response_id = resp_id,
        gene_symbol = sym,
        as.data.table(s),
        nt_pool_cp10k           = pool$cp10k,
        nt_pool_mean_cell_cp10k = pool$mean_cell_cp10k,
        nt_pool_n_cells         = pool$n_cells
      )
    }
  }
  cat(sprintf("  %-10s %-18s TSS gRNAs: %2d   NT pool cells: %6d   NT pool cp10k: %.2f\n",
              sym, resp_id, length(tss_guides), pool$n_cells, pool$cp10k))
}

per_guide_dt <- rbindlist(per_guide[seq_len(k)])
per_guide_dt[, log2fc_vs_nt_pool := log2(cp10k / nt_pool_cp10k)]
per_guide_dt[, pct_knockdown_vs_nt_pool := 100 * (1 - cp10k / nt_pool_cp10k)]
setorder(per_guide_dt, gene_symbol, series, grna_id)

pool_dt <- rbindlist(pool_rows)
setorder(pool_dt, gene_symbol)

for (path in c(opt$output_per_guide, opt$output_nt_pool)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}
fwrite(per_guide_dt, opt$output_per_guide, sep = "\t")
fwrite(pool_dt, opt$output_nt_pool, sep = "\t")

cat(sprintf("\nWrote %d per-gRNA rows to %s\n", nrow(per_guide_dt), opt$output_per_guide))
cat(sprintf("Wrote %d NT-pool rows to %s\n", nrow(pool_dt), opt$output_nt_pool))
cat(sprintf("  TSS rows: %d   NTC rows: %d\n",
            sum(per_guide_dt$series == "TSS"), sum(per_guide_dt$series == "NTC")))
