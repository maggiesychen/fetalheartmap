## Shared base-R helpers for the R step scripts.
##
## Deliberately dependency-free. The shared R library these scripts run against
## (sceptre.r_env.r_libs_user) has neither `optparse` nor `getopt`, and the
## `argparse` it does have is a reticulate wrapper around Python's argparse --
## which would make every R script depend on a working Python inside R. The
## SCEPTRE pipeline's own scripts use plain `commandArgs()` for the same reason.
##
## Usage:
##
##   source(file.path(dirname(this_file()), "r_utils.R"))
##   opt <- parse_cli_args(list(
##     input  = list(required = TRUE,  help = "input path"),
##     output = list(required = TRUE,  help = "output path"),
##     scale  = list(default = 1e4,    type = "numeric"),
##     strict = list(flag = TRUE,      help = "abort on the first warning")
##   ))
##   opt$input; opt$scale; opt$strict
##
## Option names are given with underscores and accepted on the command line
## with dashes, so `min_cells` is passed as `--min-cells`. `--help` prints the
## spec and exits 0.

parse_cli_args <- function(spec, args = commandArgs(trailingOnly = TRUE)) {
  to_flag <- function(name) paste0("--", gsub("_", "-", name))

  print_usage <- function() {
    script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
    cat("Usage:", if (length(script)) basename(script) else "script.R",
        "[options]\n\nOptions:\n")
    for (name in names(spec)) {
      entry <- spec[[name]]
      is_flag <- isTRUE(entry$flag)
      label <- if (is_flag) to_flag(name) else paste(to_flag(name), "VALUE")
      bits <- c()
      if (isTRUE(entry$required)) bits <- c(bits, "required")
      if (!is.null(entry$default)) {
        bits <- c(bits, paste0("default: ", paste(entry$default, collapse = ",")))
      }
      suffix <- if (length(bits)) paste0(" [", paste(bits, collapse = "; "), "]") else ""
      cat(sprintf("  %-32s %s%s\n", label,
                  if (is.null(entry$help)) "" else entry$help, suffix))
    }
    cat("  --help                           show this message and exit\n")
  }

  if (any(args %in% c("--help", "-h"))) {
    print_usage()
    quit(status = 0L)
  }

  ## Seed with declared defaults; flags default to FALSE.
  opt <- list()
  for (name in names(spec)) {
    entry <- spec[[name]]
    opt[[name]] <- if (isTRUE(entry$flag)) FALSE else entry$default
  }

  flag_to_name <- setNames(names(spec), vapply(names(spec), to_flag, character(1)))

  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) {
      print_usage()
      stop(sprintf("unexpected positional argument: %s", token), call. = FALSE)
    }
    ## Support both `--key value` and `--key=value`.
    if (grepl("=", token, fixed = TRUE)) {
      key <- sub("=.*$", "", token)
      inline <- sub("^[^=]*=", "", token)
    } else {
      key <- token
      inline <- NULL
    }
    if (!key %in% names(flag_to_name)) {
      print_usage()
      stop(sprintf("unknown option: %s", key), call. = FALSE)
    }
    name <- flag_to_name[[key]]
    entry <- spec[[name]]

    if (isTRUE(entry$flag)) {
      opt[[name]] <- TRUE
      if (!is.null(inline)) {
        opt[[name]] <- tolower(inline) %in% c("1", "true", "yes")
      }
      i <- i + 1L
      next
    }

    if (!is.null(inline)) {
      value <- inline
      i <- i + 1L
    } else {
      if (i + 1L > length(args)) {
        stop(sprintf("%s requires a value", key), call. = FALSE)
      }
      value <- args[[i + 1L]]
      i <- i + 2L
    }

    type <- if (is.null(entry$type)) "character" else entry$type
    opt[[name]] <- switch(
      type,
      numeric = as.numeric(value),
      integer = as.integer(value),
      character = value,
      stop(sprintf("unsupported type for %s: %s", name, type), call. = FALSE)
    )
    if (type %in% c("numeric", "integer") && is.na(opt[[name]])) {
      stop(sprintf("%s expects a %s, got %s", key, type, value), call. = FALSE)
    }
  }

  missing <- names(spec)[vapply(
    names(spec),
    function(n) isTRUE(spec[[n]]$required) &&
      (is.null(opt[[n]]) || !nzchar(as.character(opt[[n]])[1])),
    logical(1)
  )]
  if (length(missing)) {
    print_usage()
    stop(sprintf("missing required option(s): %s",
                 paste(vapply(missing, to_flag, character(1)), collapse = ", ")),
         call. = FALSE)
  }

  opt
}


split_csv <- function(x) {
  ## Comma-separated option value -> character vector; "" -> character(0).
  if (is.null(x) || !length(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1]])
}


this_script_dir <- function() {
  ## Directory of the running Rscript, so a script can source its siblings.
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  getwd()
}


strip_ensembl_version <- function(x) sub("\\..*", "", x)
