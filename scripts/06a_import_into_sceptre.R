#!/usr/bin/env Rscript
## Step 6a -- import the CellRanger matrices into an ondisc-backed sceptre object.
##
## Reads the feature-corrected CellRanger exports from stage 00 (one directory
## per sublibrary) plus the gRNA target table, and writes an on-disk-backed
## sceptre object. On-disk backing (ondisc) keeps the full 240k x 26k matrix off
## the heap, which is what makes the trans analysis tractable.
##
## import_data_from_cellranger() indexes features **positionally** across the
## directories, so every sublibrary must present an identical features.tsv.gz.
## Stage 00's validate_cellranger_features rule enforces that; if it did not
## hold, results would be silently wrong rather than erroring.
##
## Requires the sceptre/ondisc R library configured by sceptre.r_env.

## Locate and source the shared base-R helpers next to this script.
this_script_dir_bootstrap <- function() {
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  getwd()
}
source(file.path(this_script_dir_bootstrap(), "r_utils.R"))

suppressPackageStartupMessages({
  library(sceptre)
})

opt <- parse_cli_args(list(
  cellranger_dirs = list(required = TRUE,
    help = "Comma-separated filtered_feature_bc_matrix directories, in import order"),
  guide_targets   = list(required = TRUE,
    help = "TSV with grna_id / grna_target (and optional chr/start/end)"),
  moi             = list(default = "high", help = "high or low"),
  output_dir      = list(required = TRUE, help = "Directory for the ondisc backing files"),
  output          = list(required = TRUE, help = "Path of the sceptre object .rds to write")
))

directories <- strsplit(opt$cellranger_dirs, ",", fixed = TRUE)[[1]]
directories <- trimws(directories)
cat("Importing", length(directories), "sublibrary directories:\n")
for (d in directories) {
  if (!dir.exists(d)) stop(sprintf("directory does not exist: %s", d))
  cat("  ", d, "\n")
}

grna_target_data_frame <- read.table(opt$guide_targets, header = TRUE, sep = "\t")
required_cols <- c("grna_id", "grna_target")
missing <- setdiff(required_cols, colnames(grna_target_data_frame))
if (length(missing) > 0) {
  stop(sprintf("%s is missing columns: %s",
               opt$guide_targets, paste(missing, collapse = ", ")))
}
cat("\ngRNA target table:", nrow(grna_target_data_frame), "gRNAs across",
    length(unique(grna_target_data_frame$grna_target)), "targets\n")
cat("  non-targeting gRNAs:",
    sum(grna_target_data_frame$grna_target == "non-targeting"), "\n")

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

cat("\nRunning import_data_from_cellranger(moi =", opt$moi, ")...\n")
sceptre_object <- import_data_from_cellranger(
  directories            = directories,
  moi                    = opt$moi,
  grna_target_data_frame = grna_target_data_frame,
  use_ondisc             = TRUE,
  directory_to_write     = opt$output_dir
)

print(sceptre_object)

saveRDS(sceptre_object, opt$output)
cat("\nWrote", opt$output, "\n")
cat("ondisc backing files in", opt$output_dir, "\n")
