# R/05_fgsea.R
# Stage 05: Fast Gene Set Enrichment Analysis (fgsea) per PAM50 subtype.
#           Gene sets: MSigDB Hallmark + KEGG Legacy (via msigdbr).
#           Ranking metric: limma moderated t-statistic from {Subtype}_vs_GTEx.
#           One .csv per subtype, all pathways retained (no pre-filter).
#
# Inputs:  results/tables/deg/{Subtype}_vs_GTEx.csv   (4 files)
#          results/tables/deg/{Subtype}_vs_NAT.csv    (4 files)
# Outputs: results/tables/fgsea/{Subtype}_fgsea.csv   (4 files)
#            Columns: contrast, collection, pathway, NES, pval, padj,
#                     size, leadingEdge

library(here)
library(dplyr)
library(fgsea)
library(msigdbr)

set.seed(42)

PAM50_LEVELS <- c("LumA", "LumB", "Her2", "Basal")
MIN_SIZE     <- 15L
MAX_SIZE     <- 500L
N_PERM       <- 1000L
FDR_CUTOFF   <- 0.05

# Contrasts to run GSEA on for each subtype.
# GTEx = absolute dysregulation vs healthy; NAT = field biology reference.
CONTRASTS <- c("vs_GTEx", "vs_NAT")

# ---------------------------------------------------------------------------
# Stage gate
# ---------------------------------------------------------------------------

#' Stop with a clear error if Stage 03 DEG outputs are missing
#'
#' @return Invisible NULL.
#' @examples
#' assert_stage03_complete()
assert_stage03_complete <- function() {
  message("── Stage 05 | stage gate ────────────────────────────────────────────────")

  deg_dir  <- here::here("results", "tables", "deg")
  required <- c(
    paste0(PAM50_LEVELS, "_vs_GTEx.csv"),
    paste0(PAM50_LEVELS, "_vs_NAT.csv")
  )
  missing <- required[!file.exists(file.path(deg_dir, required))]
  if (length(missing) > 0)
    stop("Missing DEG files: ", paste(missing, collapse = ", "),
         "\nRun R/03_deg.R first.")

  message("  DEG files found for all subtypes × both contrasts.")
  message("── Stage 05 | gate passed ───────────────────────────────────────────────")
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Gene sets
# ---------------------------------------------------------------------------

#' Load MSigDB Hallmark and KEGG Legacy gene sets as a named list
#'
#' Uses \code{msigdbr} (local, no network after first download) to retrieve
#' human gene sets. Returns a single named list suitable for \code{fgsea()}.
#' Pathway names are prefixed with their collection (\code{HALLMARK_} or
#' \code{KEGG_}) to avoid name collisions.
#'
#' @return Named list of character vectors. Names = pathway IDs,
#'   values = HGNC gene symbols.
#' @examples
#' gene_sets <- load_gene_sets()
load_gene_sets <- function() {
  message("── Stage 05 | loading gene sets ─────────────────────────────────────────")

  hallmark <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H") |>
    dplyr::select(gs_name, gene_symbol) |>
    dplyr::group_by(gs_name) |>
    dplyr::summarise(genes = list(unique(gene_symbol)), .groups = "drop")

  kegg <- msigdbr::msigdbr(species = "Homo sapiens",
                             collection    = "C2",
                             subcollection = "CP:KEGG_LEGACY") |>
    dplyr::select(gs_name, gene_symbol) |>
    dplyr::group_by(gs_name) |>
    dplyr::summarise(genes = list(unique(gene_symbol)), .groups = "drop")

  sets <- dplyr::bind_rows(hallmark, kegg)
  gene_sets <- setNames(sets$genes, sets$gs_name)

  message("  Hallmark pathways : ", nrow(hallmark))
  message("  KEGG Legacy pathways: ", nrow(kegg))
  message("  Total gene sets   : ", length(gene_sets))
  message("── Stage 05 | gene sets loaded ──────────────────────────────────────────")
  gene_sets
}

# ---------------------------------------------------------------------------
# Rank vector
# ---------------------------------------------------------------------------

#' Build a named t-statistic rank vector from a DEG result CSV
#'
#' Reads the specified DEG CSV, drops genes with NA t-statistics, and
#' resolves duplicate gene names by keeping the row with the largest
#' absolute t-statistic. Returns a named numeric vector sorted in
#' descending order, ready for \code{fgsea()}.
#'
#' @param subtype  Character. One of PAM50_LEVELS.
#' @param contrast Character. One of \code{"vs_GTEx"} or \code{"vs_NAT"}.
#' @return Named numeric vector. Names = HGNC symbols, values = t-statistics.
#' @examples
#' ranks <- build_rank_vector("LumA", "vs_GTEx")
build_rank_vector <- function(subtype, contrast) {
  path <- here::here("results", "tables", "deg",
                     paste0(subtype, "_", contrast, ".csv"))
  deg  <- read.csv(path, stringsAsFactors = FALSE)

  # Drop NAs and resolve duplicate gene symbols
  deg  <- deg[!is.na(deg$t), ]
  deg  <- deg[order(-abs(deg$t)), ]
  deg  <- deg[!duplicated(deg$gene), ]

  ranks <- setNames(deg$t, deg$gene)
  sort(ranks, decreasing = TRUE)
}

# ---------------------------------------------------------------------------
# Run fgsea
# ---------------------------------------------------------------------------

#' Run fgsea for one subtype × one contrast and return tidy results
#'
#' Runs \code{fgsea::fgsea()} with \code{nPermSimple = N_PERM}.
#' The \code{leadingEdge} list column is collapsed to a semicolon-separated
#' string for CSV compatibility. All pathways are returned (no FDR pre-filter);
#' the \code{padj} column should be used for filtering downstream.
#'
#' @param gene_sets Named list of gene symbol vectors.
#' @param subtype   Character. One of PAM50_LEVELS.
#' @param contrast  Character. One of \code{"vs_GTEx"} or \code{"vs_NAT"}.
#' @return data.frame: contrast, collection, pathway, NES, pval, padj,
#'   size, leadingEdge.
#' @examples
#' res <- run_fgsea_contrast(gene_sets, "LumA", "vs_GTEx")
run_fgsea_contrast <- function(gene_sets, subtype, contrast) {
  message("  ", subtype, " ", contrast, " ...")

  ranks <- build_rank_vector(subtype, contrast)

  set.seed(42)
  res <- fgsea::fgsea(
    pathways    = gene_sets,
    stats       = ranks,
    minSize     = MIN_SIZE,
    maxSize     = MAX_SIZE,
    nPermSimple = N_PERM,
    eps         = 0        # accurate p-values for pathways with p < 1e-50
  )

  n_sig <- sum(res$padj < FDR_CUTOFF, na.rm = TRUE)
  message("    tested: ", nrow(res), " pathways | significant (padj<",
          FDR_CUTOFF, "): ", n_sig)
  if (n_sig < 5L)
    warning(subtype, " ", contrast, ": only ", n_sig,
            " significant pathways — check gene set overlap.")

  # Determine collection from pathway name prefix
  res$collection <- dplyr::case_when(
    startsWith(res$pathway, "HALLMARK_") ~ "Hallmark",
    startsWith(res$pathway, "KEGG_")     ~ "KEGG_Legacy",
    TRUE                                 ~ "Other"
  )

  # Collapse leadingEdge list to semicolon string for CSV output
  res$leadingEdge <- vapply(res$leadingEdge,
                             function(x) paste(x, collapse = ";"),
                             character(1L))

  data.frame(
    contrast    = paste0(subtype, "_", contrast),
    collection  = res$collection,
    pathway     = res$pathway,
    NES         = res$NES,
    pval        = res$pval,
    padj        = res$padj,
    size        = res$size,
    leadingEdge = res$leadingEdge,
    stringsAsFactors = FALSE,
    row.names = NULL
  ) |>
    dplyr::arrange(padj, desc(abs(NES)))
}

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

#' Write fgsea results for one subtype to results/tables/fgsea/{Subtype}_fgsea.csv
#'
#' Combines results from all contrasts for one subtype into a single CSV.
#' Both GTEx and NAT contrasts are included as rows, distinguished by the
#' \code{contrast} column.
#'
#' @param results data.frame. Combined output of \code{run_fgsea_contrast()}
#'   calls for one subtype.
#' @param subtype Character. Used as filename stem.
#' @return Invisible character path.
#' @examples
#' save_fgsea(results, "LumA")
save_fgsea <- function(results, subtype) {
  out_dir <- here::here("results", "tables", "fgsea")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, paste0(subtype, "_fgsea.csv"))
  write.csv(results, file = path, row.names = FALSE, quote = FALSE)

  n_sig <- sum(results$padj < FDR_CUTOFF, na.rm = TRUE)
  message("  saved: ", basename(path),
          " (", nrow(results), " pathway-rows, ",
          n_sig, " at padj<", FDR_CUTOFF, ")")
  invisible(path)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

assert_stage03_complete()

gene_sets <- load_gene_sets()

message("── Stage 05 | running fgsea ─────────────────────────────────────────────")

for (st in PAM50_LEVELS) {
  message("  subtype: ", st)
  contrast_results <- lapply(CONTRASTS, function(ct) {
    run_fgsea_contrast(gene_sets, st, ct)
  })
  combined <- dplyr::bind_rows(contrast_results)
  save_fgsea(combined, st)
}

# Verify 4 files written
fgsea_files <- list.files(here::here("results", "tables", "fgsea"),
                           pattern = "_fgsea\\.csv$")
message("── Stage 05 | output summary ────────────────────────────────────────────")
message("  fgsea files written: ", length(fgsea_files), " / 4 expected")
stopifnot("expected 4 fgsea files" = length(fgsea_files) == 4L)

message("✅ Stage 05 complete. Proceed to R/06_wgcna.R")
