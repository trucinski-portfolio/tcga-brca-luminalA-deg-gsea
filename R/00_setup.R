# R/00_setup.R
# Stage 00: renv validation, directory scaffold, Xena data download + input validation
#
# Contrasts active: none (setup stage only)
#
# Inputs:  none (downloads raw data if absent)
# Outputs: data/raw/preprocessing_inputs/HiSeqV2_PANCAN.gz
#          data/raw/preprocessing_inputs/BRCA_clinicalMatrix.tsv
#          all project directories created

library(here)

# ---------------------------------------------------------------------------
# renv bootstrap
# ---------------------------------------------------------------------------

#' Bootstrap renv for reproducible package management
#'
#' Checks that renv is available and the project library is synchronized with
#' renv.lock. Installs renv itself if absent, then calls renv::restore() when
#' the library is out of sync.
#'
#' @return Invisible NULL. Called for side effects.
#'
#' @examples
#' bootstrap_renv()
bootstrap_renv <- function() {
  message("── Stage 00 | renv bootstrap ──────────────────────────────────────────")

  if (!requireNamespace("renv", quietly = TRUE)) {
    message("renv not found — installing from CRAN.")
    install.packages("renv", repos = "https://cloud.r-project.org")
  }

  lockfile <- here::here("renv.lock")
  if (!file.exists(lockfile)) {
    stop(
      "renv.lock not found at: ", lockfile, "\n",
      "Run renv::snapshot() from the repo root to create it."
    )
  }

  status <- renv::status(project = here::here())
  if (!isTRUE(status$synchronized)) {
    message("renv library out of sync — restoring from renv.lock ...")
    renv::restore(project = here::here(), prompt = FALSE)
  } else {
    message("renv library is synchronized.")
  }

  message("R version : ", getRversion())
  message("Library   : ", .libPaths()[1])
  message("Lockfile  : ", normalizePath(lockfile))
  message("── Stage 00 | renv bootstrap complete ─────────────────────────────────")
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Directory scaffold
# ---------------------------------------------------------------------------

#' Create all required project directories
#'
#' Creates every directory in the v3.0 project layout if it does not already
#' exist. Safe to re-run (idempotent).
#'
#' @return Invisible character vector of directory paths.
#'
#' @examples
#' create_project_dirs()
create_project_dirs <- function() {
  message("── Stage 00 | creating directories ────────────────────────────────────")

  dirs <- c(
    here::here("data", "raw", "preprocessing_inputs"),
    here::here("data", "processed"),
    here::here("results", "figures", "qc"),
    here::here("results", "figures", "clustering"),
    here::here("results", "figures", "deg"),
    here::here("results", "figures", "gsea"),
    here::here("results", "figures", "wgcna"),
    here::here("results", "figures", "survival"),
    here::here("results", "tables", "deg"),
    here::here("results", "tables", "fgsea"),
    here::here("results", "tables", "signatures"),
    here::here("results", "tables", "survival")
  )

  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
      message("  created: ", d)
    }
  }

  message("── Stage 00 | directories ready ────────────────────────────────────────")
  invisible(dirs)
}

# ---------------------------------------------------------------------------
# Xena download
# ---------------------------------------------------------------------------

#' Download a file from URL with retry logic
#'
#' Attempts to download a file up to \code{retries} times, pausing
#' \code{pause_sec} seconds between attempts. Skips download if the destination
#' file already exists and is non-empty.
#'
#' @param url Character. Source URL.
#' @param dest Character. Destination file path.
#' @param retries Integer. Number of download attempts. Default 3.
#' @param pause_sec Numeric. Seconds to wait between retries. Default 2.
#'
#' @return Invisible character scalar of the destination path.
#'
#' @examples
#' download_with_retry(
#'   url  = "https://example.com/file.gz",
#'   dest = here::here("data", "raw", "preprocessing_inputs", "file.gz")
#' )
download_with_retry <- function(url, dest, retries = 3L, pause_sec = 2) {
  if (file.exists(dest) && file.size(dest) > 0) {
    message("  already present (skipping): ", basename(dest))
    return(invisible(dest))
  }

  for (attempt in seq_len(retries)) {
    message("  downloading (attempt ", attempt, "/", retries, "): ", basename(dest))
    tryCatch(
      {
        utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE)
        if (file.exists(dest) && file.size(dest) > 0) {
          message("  download complete: ", basename(dest))
          return(invisible(dest))
        }
        message("  download produced empty file — retrying ...")
      },
      error = function(e) {
        message("  error on attempt ", attempt, ": ", conditionMessage(e))
      }
    )
    if (attempt < retries) Sys.sleep(pause_sec)
  }

  stop(
    "Failed to download after ", retries, " attempts: ", basename(dest), "\n",
    "URL: ", url
  )
}

#' Download TCGA-BRCA expression and clinical matrices from UCSC Xena
#'
#' Downloads the TOIL RSEM TPM expression matrix (HiSeqV2_PANCAN) and the
#' BRCA clinical matrix from the UCSC Xena data hub. Both files are required
#' by Stage 01 preprocessing. Downloads are skipped when files already exist.
#'
#' @return Invisible named character vector with paths to \code{expr} and
#'   \code{clinical} files.
#'
#' @examples
#' paths <- download_xena_inputs()
download_xena_inputs <- function() {
  message("── Stage 00 | downloading Xena inputs ─────────────────────────────────")

  inputs <- list(
    expr = list(
      url  = paste0(
        "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/",
        "TCGA.BRCA.sampleMap%2FHiSeqV2_PANCAN.gz"
      ),
      dest = here::here("data", "raw", "preprocessing_inputs", "HiSeqV2_PANCAN.gz")
    ),
    clinical = list(
      url  = paste0(
        "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/",
        "TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix"
      ),
      dest = here::here(
        "data", "raw", "preprocessing_inputs", "BRCA_clinicalMatrix.tsv"
      )
    )
  )

  paths <- vapply(inputs, function(x) {
    download_with_retry(url = x$url, dest = x$dest)
    x$dest
  }, character(1))

  message("── Stage 00 | Xena inputs downloaded ──────────────────────────────────")
  invisible(paths)
}

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

#' Validate that required raw input files are present and non-empty
#'
#' Checks each required input file and stops with a clear diagnostic message
#' if any file is absent or empty. Called after \code{download_xena_inputs()}
#' to provide a stage gate before Stage 01 can proceed.
#'
#' @return Invisible NULL. Stops with an error if validation fails.
#'
#' @examples
#' validate_inputs()
validate_inputs <- function() {
  message("── Stage 00 | validating inputs ────────────────────────────────────────")

  required <- c(
    here::here("data", "raw", "preprocessing_inputs", "HiSeqV2_PANCAN.gz"),
    here::here("data", "raw", "preprocessing_inputs", "BRCA_clinicalMatrix.tsv")
  )

  missing_files <- character(0)
  for (f in required) {
    if (!file.exists(f) || file.size(f) == 0) {
      missing_files <- c(missing_files, f)
      message("  MISSING: ", f)
    } else {
      sz_mb <- round(file.size(f) / 1024^2, 1)
      message("  OK (", sz_mb, " MB): ", basename(f))
    }
  }

  if (length(missing_files) > 0) {
    stop(
      length(missing_files), " required input file(s) missing or empty.\n",
      "Run download_xena_inputs() or re-run 00_setup.R to fetch them.\n",
      "Missing:\n  ", paste(missing_files, collapse = "\n  ")
    )
  }

  message("── Stage 00 | all inputs validated ────────────────────────────────────")
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

bootstrap_renv()
create_project_dirs()
download_xena_inputs()
validate_inputs()

message("✅ Stage 00 complete. Proceed to R/01_preprocessing.R")
