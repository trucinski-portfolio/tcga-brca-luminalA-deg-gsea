# R/00_setup.R
# Stage 00: renv verification, phenotype download, column inspection,
#           3-way sample filter, and metadata cache.
#
# Inputs:  none
# Outputs: data/processed/metadata_breast.rds  ($pheno data.frame, group-labeled)

library(here)
library(dplyr)
library(UCSCXenaTools)

TOIL_HUB      <- "toilHub"
PHENO_DATASET <- "TcgaTargetGTEX_phenotype"
EXPECTED_N    <- c(Tumor = 1100L, NAT = 113L, GTEx = 200L)

# ---------------------------------------------------------------------------
# renv
# ---------------------------------------------------------------------------

#' Verify renv is installed and restore the project library from renv.lock
#'
#' Calls \code{renv::restore()} which is a no-op when the library is already
#' synchronized, and installs any missing packages otherwise.
#'
#' @return Invisible NULL.
#' @examples
#' verify_renv()
verify_renv <- function() {
  message("── Stage 00 | renv ──────────────────────────────────────────────────────")
  if (!requireNamespace("renv", quietly = TRUE))
    stop("renv not installed. Run: install.packages('renv')")
  if (!file.exists(here::here("renv.lock")))
    stop("renv.lock not found. Run renv::snapshot() from the repo root.")

  renv::restore(prompt = FALSE)
  message("  R ", getRversion(), " | library: ", .libPaths()[1])
  message("── Stage 00 | renv OK ───────────────────────────────────────────────────")
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Phenotype download (includes hub connectivity check)
# ---------------------------------------------------------------------------

#' Download TcgaTargetGTEX_phenotype from the TOIL hub
#'
#' Follows the authoritative UCSCXenaTools pipeline from CLAUDE.md exactly:
#' \code{XenaGenerate → XenaFilter → XenaQuery → XenaDownload → XenaPrepare}.
#' Hub connectivity is verified implicitly — any network failure surfaces as a
#' clear error via \code{tryCatch}.
#'
#' @return data.frame of the full phenotype table.
#' @examples
#' pheno <- download_phenotype()
download_phenotype <- function() {
  message("── Stage 00 | phenotype download ───────────────────────────────────────")

  pheno <- tryCatch({
    pheno_dl <- XenaGenerate(subset = XenaHostNames == TOIL_HUB) |>
      XenaFilter(filterDatasets = PHENO_DATASET) |>
      XenaQuery() |>
      XenaDownload(destdir = here::here("data", "raw"))
    XenaPrepare(pheno_dl)
  }, error = function(e) {
    stop("Failed to reach toilHub or download phenotype — check connectivity.\n",
         "Detail: ", conditionMessage(e))
  })

  message("  rows: ", nrow(pheno), " | cols: ", ncol(pheno))
  message("── Stage 00 | phenotype downloaded ─────────────────────────────────────")
  pheno
}

# ---------------------------------------------------------------------------
# Column inspection
# ---------------------------------------------------------------------------

#' Print phenotype column names and detect the GTEx primary-site field
#'
#' R may import \code{_primary_site} (leading underscore) as \code{X_primary_site}.
#' Prints all column names so the confirmed name can be recorded in CLAUDE.md,
#' then returns whichever variant is present in this download.
#'
#' @param pheno data.frame. Full phenotype table from \code{download_phenotype()}.
#' @return Character scalar — the confirmed GTEx site column name.
#' @examples
#' site_col <- inspect_columns(pheno)
inspect_columns <- function(pheno) {
  message("── Stage 00 | column inspection ────────────────────────────────────────")
  message("  colnames(pheno):\n  ", paste(colnames(pheno), collapse = "\n  "))

  candidates <- c("X_primary_site", "_primary_site")
  site_col   <- intersect(candidates, colnames(pheno))[1]

  if (is.na(site_col))
    stop("Neither 'X_primary_site' nor '_primary_site' found in phenotype. ",
         "Update CLAUDE.md with the correct GTEx site column name.")
  if (sum(candidates %in% colnames(pheno)) > 1)
    warning("Both column name variants present — using '", site_col, "'. Verify.")

  message("  GTEx site column confirmed: '", site_col, "'")
  message("  >> Record this in CLAUDE.md under 'Confirmed GTEx field name'")
  message("── Stage 00 | column inspection complete ───────────────────────────────")
  site_col
}

# ---------------------------------------------------------------------------
# 3-way filter + group labeling
# ---------------------------------------------------------------------------

#' Apply the 3-way OR filter and assign Tumor / NAT / GTEx group labels
#'
#' Implements the exact OR filter from CLAUDE.md using confirmed column names:
#' \itemize{
#'   \item TCGA arm — \code{`primary disease or tissue` == "breast invasive carcinoma"}
#'         (captures Tumor barcode 01 and NAT barcode 11 in one pass)
#'   \item GTEx arm — \code{`_study` == "GTEX" & `_primary_site` == "Breast"}
#' }
#' Groups are assigned by TCGA barcode positions 14–15 for TCGA samples and
#' by study origin for GTEx. Unclassified rows are dropped with a warning.
#'
#' @param pheno    data.frame. Full phenotype table.
#' @param site_col Character. Confirmed GTEx site column from \code{inspect_columns()}.
#' @return data.frame with column \code{group} (\code{"Tumor"}, \code{"NAT"},
#'   \code{"GTEx"}) — no "Other" rows.
#' @examples
#' breast <- filter_breast(pheno, site_col)
filter_breast <- function(pheno, site_col) {
  message("── Stage 00 | 3-way breast filter ──────────────────────────────────────")

  breast <- pheno |>
    filter(
      `primary disease or tissue` == "Breast Invasive Carcinoma" |
      (`_study` == "GTEX" & `_primary_site` == "Breast")
    ) |>
    mutate(group = case_when(
      substr(sample, 14, 15) == "01" ~ "Tumor",
      substr(sample, 14, 15) == "11" ~ "NAT",
      `_study` == "GTEX"             ~ "GTEx",
      TRUE                           ~ "Other"
    ))

  n_other <- sum(breast$group == "Other")
  if (n_other > 0) {
    warning(n_other, " sample(s) could not be classified and will be dropped.")
    breast <- filter(breast, group != "Other")
  }

  counts <- table(breast$group)
  for (grp in names(EXPECTED_N)) {
    n <- as.integer(counts[grp])
    message(sprintf("  %-8s: %d", grp, if (is.na(n)) 0L else n))
    if (is.na(n) || n < EXPECTED_N[[grp]] * 0.5)
      warning("Group '", grp, "' n=", if (is.na(n)) 0L else n,
              " is below 50% of expected (~", EXPECTED_N[[grp]], "). Verify data.")
  }

  message("── Stage 00 | filter complete ───────────────────────────────────────────")
  breast
}

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

#' Save the labeled breast metadata to data/processed/metadata_breast.rds
#'
#' Writes the Stage 00 output consumed by \code{R/01_preprocessing.R}.
#'
#' @param breast data.frame. Labeled output of \code{filter_breast()}.
#' @return Invisible character path.
#' @examples
#' save_metadata(breast)
save_metadata <- function(breast) {
  message("── Stage 00 | saving ────────────────────────────────────────────────────")
  out <- here::here("data", "processed", "metadata_breast.rds")
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(pheno = breast), file = out)
  message("  saved: ", out)
  message("── Stage 00 | done ──────────────────────────────────────────────────────")
  invisible(out)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

verify_renv()
pheno    <- download_phenotype()
site_col <- inspect_columns(pheno)
breast   <- filter_breast(pheno, site_col)
rm(pheno)
save_metadata(breast)

message("✅ Stage 00 complete. Proceed to R/01_preprocessing.R")
