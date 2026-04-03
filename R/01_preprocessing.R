# R/01_preprocessing.R
# Stage 01: Load metadata, download + subset expression matrix, map Ensembl → HGNC,
#           assign 3-way group labels and PAM50, build contrast list, save cohort.rds.
#
# Inputs:  data/processed/metadata_breast.rds   ($pheno from Stage 00)
# Outputs: data/processed/cohort.rds
#            $expr      — matrix (HGNC symbols × sample IDs), log2(TPM+0.001)
#            $meta      — data.frame (sample_id, group, PAM50, sample_type)
#            $contrasts — named character vector, length == 18

library(here)
library(dplyr)
library(stringr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(UCSCXenaTools)

TOIL_HUB      <- "toilHub"
EXPR_DATASET  <- "TcgaTargetGtex_rsem_gene_tpm"
PAM50_LEVELS  <- c("LumA", "LumB", "Her2", "Basal")
GROUP_LEVELS  <- c("Tumor", "NAT", "GTEx")

# ---------------------------------------------------------------------------
# Stage gate
# ---------------------------------------------------------------------------

#' Stop with a clear error if Stage 00 output is missing
#'
#' @return Invisible NULL.
#' @examples
#' assert_stage00_complete()
assert_stage00_complete <- function() {
  message("── Stage 01 | stage gate ────────────────────────────────────────────────")
  path <- here::here("data", "processed", "metadata_breast.rds")
  if (!file.exists(path))
    stop("Stage 00 output missing: ", path, "\nRun R/00_setup.R first.")
  message("  metadata_breast.rds found.")
  message("── Stage 01 | gate passed ───────────────────────────────────────────────")
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Load metadata
# ---------------------------------------------------------------------------

#' Load Stage 00 phenotype output and extract breast sample IDs
#'
#' Reads \code{metadata_breast.rds} and returns the \code{$pheno} data.frame
#' produced by Stage 00, which already contains all three sample groups
#' (Tumor barcode-01, NAT barcode-11, GTEx Healthy) and the \code{group}
#' column.
#'
#' @return data.frame with columns: sample, group, detailed_category,
#'   primary disease or tissue, _primary_site, _sample_type, _gender, _study.
#' @examples
#' pheno <- load_metadata()
load_metadata <- function() {
  message("── Stage 01 | loading metadata ─────────────────────────────────────────")
  obj   <- readRDS(here::here("data", "processed", "metadata_breast.rds"))
  pheno <- obj$pheno
  message("  samples: ", nrow(pheno),
          " | groups: ", paste(sort(unique(pheno$group)), collapse = ", "))
  message("── Stage 01 | metadata loaded ───────────────────────────────────────────")
  pheno
}

# ---------------------------------------------------------------------------
# Expression download
# ---------------------------------------------------------------------------

#' Download TcgaTargetGtex_rsem_gene_tpm and subset to breast samples
#'
#' Follows the authoritative UCSCXenaTools pipeline from CLAUDE.md exactly.
#' The full matrix (~19 K samples) is downloaded once and cached in
#' \code{data/raw/}; subsequent runs reuse the cached file. After
#' \code{XenaPrepare()} the matrix is subsetted to \code{breast_samples} and
#' returned with Ensembl IDs as row names.
#'
#' Data are already in log2(TPM + 0.001) — do NOT re-log.
#'
#' @param breast_samples Character vector of sample IDs to retain.
#' @return Numeric matrix, Ensembl IDs × breast sample IDs.
#' @examples
#' expr <- download_expression(breast_samples)
download_expression <- function(breast_samples) {
  message("── Stage 01 | downloading expression matrix ─────────────────────────────")
  message("  dataset  : ", EXPR_DATASET)
  message("  samples to keep: ", length(breast_samples))

  # The full matrix has ~19 K columns; vroom's default 128 KB line buffer is too
  # small. Raise to 256 MB before XenaPrepare() reads the file.
  Sys.setenv("VROOM_CONNECTION_SIZE" = "268435456")

  expr_raw <- tryCatch({
    dl <- XenaGenerate(subset = XenaHostNames == TOIL_HUB) |>
      XenaFilter(filterDatasets = EXPR_DATASET) |>
      XenaQuery() |>
      XenaDownload(destdir = here::here("data", "raw"))
    XenaPrepare(dl)
  }, error = function(e) {
    stop("Expression download failed — check connectivity.\nDetail: ",
         conditionMessage(e))
  })

  gene_col <- names(expr_raw)[1]
  genes    <- expr_raw[[gene_col]]
  mat      <- as.matrix(expr_raw[, -1, drop = FALSE])
  rownames(mat) <- genes
  storage.mode(mat) <- "double"

  # Subset columns to breast samples only
  keep <- intersect(colnames(mat), breast_samples)
  message("  matched breast samples in matrix: ", length(keep),
          " / ", length(breast_samples))
  if (length(keep) == 0)
    stop("No breast sample IDs matched expression matrix columns. ",
         "Check that sample IDs are consistent between phenotype and expression data.")
  mat <- mat[, keep, drop = FALSE]

  message("  expression matrix: ", nrow(mat), " genes × ", ncol(mat), " samples")
  message("── Stage 01 | expression downloaded ────────────────────────────────────")
  mat
}

# ---------------------------------------------------------------------------
# Ensembl → HGNC mapping
# ---------------------------------------------------------------------------

#' Map Ensembl gene IDs to HGNC symbols via org.Hs.eg.db
#'
#' Uses the local \code{org.Hs.eg.db} Bioconductor annotation package —
#' no network calls required. TOIL row IDs carry a version suffix
#' (e.g. \code{ENSG00000000003.14}); the suffix is stripped before lookup
#' and the original versioned ID is used as the join key. Genes with no
#' HGNC symbol are dropped. When multiple Ensembl IDs map to the same
#' HGNC symbol, expression values are collapsed to the row-wise mean.
#'
#' @param expr Numeric matrix. Rows are Ensembl IDs (versioned or bare),
#'   columns are sample IDs.
#' @return Numeric matrix. Rows are unique HGNC symbols, columns unchanged.
#' @examples
#' expr_hgnc <- map_ensembl_to_hgnc(expr)
map_ensembl_to_hgnc <- function(expr) {
  message("── Stage 01 | Ensembl → HGNC mapping ───────────────────────────────────")

  ensembl_ids      <- rownames(expr)
  ensembl_ids_base <- sub("\\.\\d+$", "", ensembl_ids)   # strip version suffix

  message("  mapping ", length(ensembl_ids), " Ensembl IDs via org.Hs.eg.db ...")

  symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys      = ensembl_ids_base,
    column    = "SYMBOL",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )

  # Build gene_map keyed on the original (versioned) IDs so the join below works
  gene_map <- data.frame(
    ensembl_gene_id = ensembl_ids,
    hgnc_symbol     = unname(symbols),
    stringsAsFactors = FALSE
  )

  gene_map <- gene_map[!is.na(gene_map$hgnc_symbol) & nchar(gene_map$hgnc_symbol) > 0, ]
  gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]
  message("  mapped: ", nrow(gene_map), " Ensembl→HGNC pairs")

  # Inner-join: keep only Ensembl IDs that have a mapping
  keep_ids <- intersect(ensembl_ids, gene_map$ensembl_gene_id)
  message("  Ensembl IDs with HGNC symbol: ", length(keep_ids),
          " / ", length(ensembl_ids))
  expr <- expr[keep_ids, , drop = FALSE]

  lookup <- stats::setNames(gene_map$hgnc_symbol, gene_map$ensembl_gene_id)
  rownames(expr) <- lookup[rownames(expr)]

  dup_count <- sum(duplicated(rownames(expr)))
  if (dup_count > 0) {
    message("  collapsing ", dup_count, " duplicate HGNC symbols by mean ...")
    expr <- rowsum(expr, rownames(expr)) / as.vector(table(rownames(expr)))
  }

  message("  final gene count (unique HGNC symbols): ", nrow(expr))
  message("── Stage 01 | mapping complete ─────────────────────────────────────────")
  expr
}

# ---------------------------------------------------------------------------
# Build metadata table
# ---------------------------------------------------------------------------

#' Fetch PAM50 subtype calls from the TCGA BRCA clinical matrix
#'
#' The TOIL combined phenotype (\code{TcgaTargetGTEX_phenotype}) does not
#' embed PAM50 subtypes for TCGA BRCA samples. The authoritative source is the
#' TCGA BRCA clinical matrix on \code{tcgaHub}, which contains the
#' \code{PAM50Call_RNAseq} column. Result is cached to
#' \code{data/processed/pam50_calls.rds} so subsequent runs skip the download.
#'
#' @return data.frame with columns: sample (15-character TCGA barcode), PAM50.
#' @examples
#' pam50 <- fetch_pam50()
fetch_pam50 <- function() {
  cache_path <- here::here("data", "processed", "pam50_calls.rds")
  if (file.exists(cache_path)) {
    message("  using cached PAM50 calls: ", cache_path)
    return(readRDS(cache_path))
  }

  message("  downloading TCGA BRCA clinical matrix for PAM50 ...")
  clin_raw <- tryCatch({
    dl <- XenaGenerate(subset = XenaHostNames == "tcgaHub") |>
      XenaFilter(filterDatasets = "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix") |>
      XenaQuery() |>
      XenaDownload(destdir = here::here("data", "raw"))
    XenaPrepare(dl)
  }, error = function(e) {
    stop("PAM50 clinical matrix download failed: ", conditionMessage(e))
  })

  # sampleID column holds 15-char TCGA barcodes
  pam50 <- clin_raw[, c("sampleID", "PAM50Call_RNAseq")] |>
    dplyr::rename(sample = sampleID, PAM50 = PAM50Call_RNAseq) |>
    dplyr::filter(!is.na(PAM50))

  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(pam50, cache_path)
  message("  PAM50 calls cached: ", cache_path, " (", nrow(pam50), " samples)")
  pam50
}

#' Assign PAM50 labels and build the final per-sample metadata table
#'
#' PAM50 calls are fetched from the TCGA BRCA clinical matrix (tcgaHub) and
#' joined on the 15-character sample barcode. Normal-like tumors are excluded
#' from the Tumor arm only — NAT and GTEx are never filtered by PAM50.
#' Control groups receive PAM50 = NA.
#'
#' @param pheno data.frame. Labeled output from \code{load_metadata()}.
#' @param expr_samples Character vector. Column names of the final expression matrix
#'   (used to inner-join metadata to expression).
#' @return data.frame with columns: sample_id, group, PAM50, sample_type.
#' @examples
#' meta <- build_meta(pheno, colnames(expr))
build_meta <- function(pheno, expr_samples) {
  message("── Stage 01 | building metadata ─────────────────────────────────────────")

  pam50 <- fetch_pam50()

  # TCGA barcodes in the expression matrix are 15-char (TCGA-xx-xxxx-01A etc.)
  # The clinical matrix uses 15-char barcodes too — join on those directly.
  meta <- pheno |>
    dplyr::filter(sample %in% expr_samples) |>
    dplyr::mutate(
      sample_id   = sample,
      sample_type = .data[["_sample_type"]]
    ) |>
    dplyr::left_join(pam50, by = "sample") |>
    dplyr::mutate(
      # Control groups (NAT, GTEx) keep PAM50 = NA; Tumor gets the joined call
      PAM50 = dplyr::if_else(group == "Tumor", PAM50, NA_character_)
    ) |>
    # Exclude Normal-like from Tumor arm only; never touch NAT or GTEx
    dplyr::filter(
      group != "Tumor" |
        (group == "Tumor" & PAM50 %in% PAM50_LEVELS)
    ) |>
    dplyr::select(sample_id, group, PAM50, sample_type)

  message("  sample counts after PAM50 filter:")
  print(table(meta$group, useNA = "ifany"))
  message("  PAM50 distribution (Tumor arm):")
  print(table(meta$PAM50[meta$group == "Tumor"], useNA = "ifany"))

  message("── Stage 01 | metadata built ────────────────────────────────────────────")
  meta
}

# ---------------------------------------------------------------------------
# Build contrasts
# ---------------------------------------------------------------------------

#' Build the 18-element named contrast vector matching CLAUDE.md
#'
#' Returns a named character vector of length 18. Names are the contrast
#' identifiers used as filenames in Stage 03. The four Within-subtype entries
#' are placeholders; Stage 02 populates them with sub-subcluster definitions
#' after UMAP clustering.
#'
#' @return Named character vector, length == 18.
#' @examples
#' contrasts <- build_contrasts()
build_contrasts <- function() {
  message("── Stage 01 | building contrasts ────────────────────────────────────────")

  contrast_names <- c(
    # Layer 1: Subtype vs NAT (4)
    "LumA_vs_NAT", "LumB_vs_NAT", "Her2_vs_NAT", "Basal_vs_NAT",
    # Layer 2: Subtype vs GTEx (4)
    "LumA_vs_GTEx", "LumB_vs_GTEx", "Her2_vs_GTEx", "Basal_vs_GTEx",
    # Layer 3: Subtype vs Subtype (5)
    "LumA_vs_Basal", "LumA_vs_LumB", "LumA_vs_Her2",
    "Her2_vs_Basal", "LumB_vs_Basal",
    # Layer 4: NAT vs GTEx (1)
    "NAT_vs_GTEx",
    # Layer 5: Within-subtype placeholders (4) — populated by Stage 02
    "Within_LumA", "Within_LumB", "Within_Her2", "Within_Basal"
  )

  stopifnot("contrast count must be 18" = length(contrast_names) == 18L)

  message("  contrasts defined: ", length(contrast_names))
  for (nm in contrast_names) message("    ", nm)
  message("── Stage 01 | contrasts built ───────────────────────────────────────────")
  contrast_names
}

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

#' Save cohort.rds to data/processed/
#'
#' Writes a named list with \code{$expr} (HGNC matrix), \code{$meta}
#' (data.frame), and \code{$contrasts} (length-18 named character vector).
#' This file is the required input for all downstream stages.
#'
#' @param expr     Numeric matrix. HGNC symbols × sample IDs.
#' @param meta     data.frame. sample_id, group, PAM50, sample_type.
#' @param contrasts Named character vector of length 18.
#' @return Invisible character path.
#' @examples
#' save_cohort(expr, meta, contrasts)
save_cohort <- function(expr, meta, contrasts) {
  message("── Stage 01 | saving cohort.rds ─────────────────────────────────────────")

  # Align matrix columns to meta row order (inner join on sample_id)
  shared <- intersect(colnames(expr), meta$sample_id)
  expr   <- expr[, shared, drop = FALSE]
  meta   <- meta[match(shared, meta$sample_id), , drop = FALSE]

  stopifnot(
    "expr cols must match meta rows" = identical(colnames(expr), meta$sample_id),
    "contrasts must be length 18"    = length(contrasts) == 18L
  )

  out <- here::here("data", "processed", "cohort.rds")
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(expr = expr, meta = meta, contrasts = contrasts), file = out)
  message("  saved: ", out)
  message("  $expr      : ", nrow(expr), " genes × ", ncol(expr), " samples")
  message("  $meta      : ", nrow(meta), " rows")
  message("  $contrasts : ", length(contrasts), " contrasts")
  message("── Stage 01 | save complete ──────────────────────────────────────────────")
  invisible(out)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

assert_stage00_complete()

pheno          <- load_metadata()
breast_samples <- pheno$sample

expr           <- download_expression(breast_samples)
expr           <- map_ensembl_to_hgnc(expr)
meta           <- build_meta(pheno, colnames(expr))
contrasts      <- build_contrasts()
save_cohort(expr, meta, contrasts)

message("✅ Stage 01 complete. Proceed to R/02_clustering.R")
