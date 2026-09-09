#!/usr/bin/env Rscript
## Step 6d -- export the SCEPTRE result RDS files as annotated TSVs.
##
## The Nextflow pipeline writes results_run_{power_check,discovery_analysis,
## calibration_check}.rds. This adds gene symbols and, for the power check,
## percent knockdown, and writes plain TSVs so downstream work does not need R.
##
## pct_knockdown = 100 * (1 - 2^log_2_fold_change)
##
## That transform is CONCAVE, so percentages must not be averaged across guides
## or targets and then compared with the percentage of a mean fold change. Where
## a mean is needed, average log2 fold changes and convert once at the end.
## `07f_plot_pct_knockdown_bar.py` asserts the Jensen inequality between the two
## for exactly this reason.

## Locate and source the shared base-R helpers next to this script.
this_script_dir_bootstrap <- function() {
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  getwd()
}
source(file.path(this_script_dir_bootstrap(), "r_utils.R"))


opt <- parse_cli_args(list(
  results_dir = list(required = TRUE,
    help = "SCEPTRE output directory containing results_run_*.rds"),
  ens2symbol  = list(required = TRUE, help = "TSV with ensembl_id / gene_symbol"),
  output_dir  = list(required = TRUE, help = "Where to write the TSVs"),
  prefix      = list(default = "sceptre", help = "Output filename prefix"),
  analyses    = list(default = "power_check,discovery_analysis,calibration_check",
    help = "Comma-separated analyses to export")
))

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

ens2symbol <- read.table(opt$ens2symbol, header = TRUE, sep = "\t",
                         quote = "", comment.char = "")
missing <- setdiff(c("ensembl_id", "gene_symbol"), colnames(ens2symbol))
if (length(missing) > 0) {
  stop(sprintf("%s is missing columns: %s", opt$ens2symbol,
               paste(missing, collapse = ", ")))
}
cat("Loaded", nrow(ens2symbol), "Ensembl -> symbol mappings\n\n")

analyses <- trimws(strsplit(opt$analyses, ",", fixed = TRUE)[[1]])

for (analysis in analyses) {
  rds <- file.path(opt$results_dir, sprintf("results_run_%s.rds", analysis))
  if (!file.exists(rds)) {
    cat("SKIP", analysis, "-- not found:", rds, "\n")
    next
  }
  cat("===", analysis, "===\n")
  df <- readRDS(rds)
  cat("  rows:", nrow(df), " columns:",
      paste(colnames(df), collapse = ", "), "\n")

  ## Strip the Ensembl version suffix so the mapping joins.
  df$ensembl_stripped <- sub("\\..*", "", df$response_id)
  df <- merge(df, ens2symbol, by.x = "ensembl_stripped", by.y = "ensembl_id",
              all.x = TRUE, sort = FALSE)
  n_unmapped <- sum(is.na(df$gene_symbol))
  cat("  unmapped responses:", n_unmapped, "/", nrow(df), "\n")

  if ("log_2_fold_change" %in% colnames(df)) {
    df$pct_knockdown <- 100 * (1 - 2^df$log_2_fold_change)
  }

  if (identical(analysis, "power_check") && "pass_qc" %in% colnames(df)) {
    passing <- df[df$pass_qc & !is.na(df$pct_knockdown), ]
    cat(sprintf("  on-target knockdown over %d pairs passing QC: mean %.1f%%, median %.1f%%\n",
                nrow(passing), mean(passing$pct_knockdown),
                median(passing$pct_knockdown)))
    cat(sprintf("  >= 50%% knockdown: %d / %d\n",
                sum(passing$pct_knockdown >= 50), nrow(passing)))
  }

  out <- file.path(opt$output_dir,
                   sprintf("%s_%s.tsv", opt$prefix, analysis))
  write.table(df, out, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  wrote", out, "\n\n")
}

cat("Done.\n")
