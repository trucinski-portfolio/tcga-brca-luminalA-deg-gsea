# R/03_deg.R
# Stage 03: Differential expression analysis via limma-trend.
#           One .csv per contrast, 18 files total.
#
# Inputs:  data/processed/cohort.rds
#          data/processed/umap_clusters.rds
# Outputs: results/tables/deg/{contrast_name}.csv   (18 files)
#            Columns: gene, logFC, AveExpr, t, P.Value, adj.P.Val, B
#
# Model design:
#   Main model  — 6-level group factor (LumA/LumB/Her2/Basal/NAT/GTEx)
#                 model.matrix(~0 + group), makeContrasts for 14 contrasts
#   Within model — 2-level subcluster factor per PAM50 subtype (4 contrasts)
#                 fitted independently on the subtype-specific expression subset

library(here)
library(dplyr)
library(limma)

set.seed(42)

PAM50_LEVELS <- c("LumA", "LumB", "Her2", "Basal")
FDR_CUTOFF   <- 0.05

# ---------------------------------------------------------------------------
# Stage gate
# ---------------------------------------------------------------------------

#' Stop with a clear error if Stage 01 or Stage 02 outputs are missing
#'
#' @return Invisible NULL.
#' @examples
#' assert_stage02_complete()
assert_stage02_complete <- function() {
  message("── Stage 03 | stage gate ────────────────────────────────────────────────")
  cohort_path  <- here::here("data", "processed", "cohort.rds")
  umap_path    <- here::here("data", "processed", "umap_clusters.rds")
  if (!file.exists(cohort_path))
    stop("Stage 01 output missing: ", cohort_path)
  if (!file.exists(umap_path))
    stop("Stage 02 output missing: ", umap_path)
  message("  cohort.rds found.")
  message("  umap_clusters.rds found.")
  message("── Stage 03 | gate passed ───────────────────────────────────────────────")
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

#' Load cohort.rds and umap_clusters.rds
#'
#' @return Named list: \code{$expr}, \code{$meta}, \code{$contrasts},
#'   \code{$umap} (data.frame with subcluster column).
#' @examples
#' inputs <- load_inputs()
load_inputs <- function() {
  message("── Stage 03 | loading inputs ────────────────────────────────────────────")
  cohort <- readRDS(here::here("data", "processed", "cohort.rds"))
  umap   <- readRDS(here::here("data", "processed", "umap_clusters.rds"))

  message("  $expr      : ", nrow(cohort$expr), " genes × ", ncol(cohort$expr), " samples")
  message("  $meta      : ", nrow(cohort$meta), " rows")
  message("  $umap      : ", nrow(umap$umap), " rows")
  message("── Stage 03 | inputs loaded ─────────────────────────────────────────────")

  list(
    expr      = cohort$expr,
    meta      = cohort$meta,
    contrasts = cohort$contrasts,
    umap      = umap$umap
  )
}

# ---------------------------------------------------------------------------
# Main limma model (14 contrasts)
# ---------------------------------------------------------------------------

#' Build a 6-level group label for the main limma model
#'
#' Tumor samples receive their PAM50 subtype as the label. NAT and GTEx
#' samples retain their group label. Returns a factor with levels ordered
#' as: LumA, LumB, Her2, Basal, NAT, GTEx.
#'
#' @param meta data.frame. Columns: sample_id, group, PAM50.
#' @return Named factor aligned to \code{meta$sample_id}.
#' @examples
#' group_label <- build_group_label(meta)
build_group_label <- function(meta) {
  label <- dplyr::case_when(
    meta$group == "Tumor" ~ meta$PAM50,
    TRUE                  ~ meta$group
  )
  factor(label, levels = c(PAM50_LEVELS, "NAT", "GTEx"))
}

#' Fit the main limma-trend model across all 6 groups
#'
#' Uses a ~ 0 + group design (no intercept) so every group has a coefficient
#' and any pairwise contrast can be formed via \code{makeContrasts}.
#' \code{eBayes(trend = TRUE)} is used because the data are in
#' log2(TPM+0.001) space rather than log-counts.
#'
#' @param expr        Numeric matrix. HGNC symbols × sample IDs.
#' @param group_label Factor. Length == \code{ncol(expr)}, 6 levels.
#' @return MArrayLM object with empirical Bayes statistics computed.
#' @examples
#' fit <- fit_main_model(expr, group_label)
fit_main_model <- function(expr, group_label) {
  message("── Stage 03 | fitting main limma model ──────────────────────────────────")
  design <- model.matrix(~ 0 + group_label)
  colnames(design) <- levels(group_label)

  fit  <- limma::lmFit(expr, design)
  fit  <- limma::eBayes(fit, trend = TRUE)
  message("  design: ", nrow(design), " samples × ", ncol(design), " coefficients")
  message("── Stage 03 | main model fitted ─────────────────────────────────────────")
  fit
}

#' Extract limma topTable results for one contrast as a tidy data.frame
#'
#' Applies \code{makeContrasts} for the specified contrast string, re-fits
#' with \code{contrasts.fit}, re-applies \code{eBayes(trend=TRUE)}, and
#' returns all genes (no FDR pre-filter — downstream stages filter as needed).
#'
#' @param fit          MArrayLM from \code{fit_main_model()}.
#' @param contrast_str Character. e.g. \code{"LumA - NAT"}.
#' @param design       model.matrix used in \code{fit_main_model()}.
#' @return data.frame: gene, logFC, AveExpr, t, P.Value, adj.P.Val, B.
#' @examples
#' res <- extract_contrast(fit, "LumA - NAT", design)
extract_contrast <- function(fit, contrast_str, design) {
  cm   <- limma::makeContrasts(contrasts = contrast_str,
                                levels    = colnames(design))
  fit2 <- limma::contrasts.fit(fit, cm)
  fit2 <- limma::eBayes(fit2, trend = TRUE)
  tt   <- limma::topTable(fit2, coef = 1, number = Inf, sort.by = "none")
  data.frame(
    gene      = rownames(tt),
    logFC     = tt$logFC,
    AveExpr   = tt$AveExpr,
    t         = tt$t,
    P.Value   = tt$P.Value,
    adj.P.Val = tt$adj.P.Val,
    B         = tt$B,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# ---------------------------------------------------------------------------
# Within-subtype limma model (4 contrasts)
# ---------------------------------------------------------------------------

#' Fit a limma-trend model for within-subtype sub-subcluster contrast
#'
#' Subsets expression and metadata to one PAM50 subtype, builds a 2-level
#' factor from the k-means sub-subcluster labels assigned in Stage 02, and
#' fits a limma-trend model. Returns NULL (with a warning) if fewer than
#' two populated sub-subclusters exist for this subtype.
#'
#' @param expr    Numeric matrix. HGNC symbols × all sample IDs.
#' @param meta    data.frame. Columns: sample_id, group, PAM50.
#' @param umap_df data.frame. Columns: sample_id, subcluster (from Stage 02).
#' @param subtype Character. One of PAM50_LEVELS.
#' @return data.frame (gene, logFC, AveExpr, t, P.Value, adj.P.Val, B)
#'   or NULL if contrast cannot be formed.
#' @examples
#' res <- run_within_contrast(expr, meta, umap_df, "LumA")
run_within_contrast <- function(expr, meta, umap_df, subtype) {
  message("  within-subtype: ", subtype)

  # Samples for this subtype with a subcluster assignment
  sub_umap <- umap_df[!is.na(umap_df$subcluster) &
                        startsWith(umap_df$subcluster, subtype), ]

  if (nrow(sub_umap) == 0) {
    warning("  no subcluster assignments for ", subtype, " — skipping.")
    return(NULL)
  }

  clusters <- unique(sub_umap$subcluster)
  if (length(clusters) < 2) {
    warning("  fewer than 2 sub-subclusters for ", subtype, " — skipping.")
    return(NULL)
  }

  # Align to expression matrix columns
  keep_ids <- intersect(sub_umap$sample_id, colnames(expr))
  sub_expr <- expr[, keep_ids, drop = FALSE]
  sub_clust <- sub_umap$subcluster[match(keep_ids, sub_umap$sample_id)]
  sub_clust <- factor(sub_clust, levels = sort(unique(sub_clust)))

  message("    ", levels(sub_clust)[1], " n=", sum(sub_clust == levels(sub_clust)[1]),
          " vs ", levels(sub_clust)[2], " n=", sum(sub_clust == levels(sub_clust)[2]))

  design <- model.matrix(~ 0 + sub_clust)
  colnames(design) <- levels(sub_clust)

  contrast_str <- paste0(levels(sub_clust)[2], " - ", levels(sub_clust)[1])
  cm   <- limma::makeContrasts(contrasts = contrast_str, levels = colnames(design))
  fit  <- limma::lmFit(sub_expr, design)
  fit2 <- limma::contrasts.fit(fit, cm)
  fit2 <- limma::eBayes(fit2, trend = TRUE)

  tt <- limma::topTable(fit2, coef = 1, number = Inf, sort.by = "none")
  data.frame(
    gene      = rownames(tt),
    logFC     = tt$logFC,
    AveExpr   = tt$AveExpr,
    t         = tt$t,
    P.Value   = tt$P.Value,
    adj.P.Val = tt$adj.P.Val,
    B         = tt$B,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

#' Write one contrast result to results/tables/deg/{name}.csv
#'
#' Prints a summary of significant genes (adj.P.Val < FDR_CUTOFF) and saves
#' all genes (unfiltered) so downstream stages can apply their own thresholds.
#'
#' @param res  data.frame. Output of \code{extract_contrast()} or
#'   \code{run_within_contrast()}.
#' @param name Character. Contrast name used as filename stem.
#' @return Invisible character path.
#' @examples
#' save_contrast(res, "LumA_vs_NAT")
save_contrast <- function(res, name) {
  out_dir <- here::here("results", "tables", "deg")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, paste0(name, ".csv"))
  write.csv(res, file = path, row.names = FALSE, quote = FALSE)

  n_sig <- sum(res$adj.P.Val < FDR_CUTOFF, na.rm = TRUE)
  message("  ", name, ": ", nrow(res), " genes tested, ",
          n_sig, " DEGs (adj.P.Val < ", FDR_CUTOFF, ")")
  invisible(path)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

assert_stage02_complete()

inputs      <- load_inputs()
expr        <- inputs$expr
meta        <- inputs$meta
umap_df     <- inputs$umap

group_label <- build_group_label(meta)
stopifnot("group_label length must match ncol(expr)" =
            length(group_label) == ncol(expr))

fit    <- fit_main_model(expr, group_label)
design <- model.matrix(~ 0 + group_label)
colnames(design) <- levels(group_label)

out_dir <- here::here("results", "tables", "deg")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("── Stage 03 | extracting contrasts ──────────────────────────────────────")

# Layer 1: Subtype vs NAT (4)
for (st in PAM50_LEVELS) {
  res <- extract_contrast(fit, paste0(st, " - NAT"), design)
  save_contrast(res, paste0(st, "_vs_NAT"))
}

# Layer 2: Subtype vs GTEx (4)
for (st in PAM50_LEVELS) {
  res <- extract_contrast(fit, paste0(st, " - GTEx"), design)
  save_contrast(res, paste0(st, "_vs_GTEx"))
}

# Layer 3: Subtype vs Subtype (5)
subtype_pairs <- list(
  c("LumA", "Basal"), c("LumA", "LumB"), c("LumA", "Her2"),
  c("Her2", "Basal"), c("LumB", "Basal")
)
for (pair in subtype_pairs) {
  res <- extract_contrast(fit, paste0(pair[1], " - ", pair[2]), design)
  save_contrast(res, paste0(pair[1], "_vs_", pair[2]))
}

# Layer 4: NAT vs GTEx (1)
res <- extract_contrast(fit, "NAT - GTEx", design)
save_contrast(res, "NAT_vs_GTEx")

# Layer 5: Within-subtype (4)
message("── Stage 03 | within-subtype contrasts ──────────────────────────────────")
for (st in PAM50_LEVELS) {
  res <- run_within_contrast(expr, meta, umap_df, st)
  if (!is.null(res)) {
    save_contrast(res, paste0("Within_", st))
  } else {
    message("  Within_", st, ": skipped (insufficient subclusters)")
  }
}

# Verify all 18 expected contrast files are present
expected_files <- paste0(c(
  paste0(PAM50_LEVELS, "_vs_NAT"),
  paste0(PAM50_LEVELS, "_vs_GTEx"),
  "LumA_vs_Basal", "LumA_vs_LumB", "LumA_vs_Her2",
  "Her2_vs_Basal", "LumB_vs_Basal",
  "NAT_vs_GTEx",
  paste0("Within_", PAM50_LEVELS)
), ".csv")

missing_files <- expected_files[!file.exists(file.path(out_dir, expected_files))]
message("── Stage 03 | output summary ────────────────────────────────────────────")
message("  expected contrast files present: ",
        length(expected_files) - length(missing_files), " / 18")
if (length(missing_files) > 0)
  stop("Missing DEG files: ", paste(missing_files, collapse = ", "))
stopifnot("expected 18 contrast files" = length(expected_files) == 18L)

message("✅ Stage 03 complete. Proceed to R/04_feature_selection.R")
