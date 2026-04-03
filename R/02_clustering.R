# R/02_clustering.R
# Stage 02: UMAP dimensionality reduction + within-subtype sub-subcluster detection.
#           All plots are printed to screen AND saved as .png + .pdf.
#
# Inputs:  data/processed/cohort.rds
# Outputs: data/processed/umap_clusters.rds
#            $umap          — data.frame: sample_id, UMAP1, UMAP2, group, PAM50, subcluster
#            $n_per_subtype — named integer: sub-subcluster counts per PAM50 subtype
#          results/figures/umap_by_group.{png,pdf}
#          results/figures/umap_by_subtype.{png,pdf}
#          results/figures/umap_subclusters_{subtype}.{png,pdf}   (4 files × 2 formats)

library(here)
library(dplyr)
library(ggplot2)
library(uwot)

set.seed(42)

N_HVG        <- 5000L    # highly variable genes fed to PCA
N_PCA        <- 50L      # PCA components fed to UMAP
K_SUBCLUSTER <- 2L       # k-means k within each PAM50 subtype
PAM50_LEVELS <- c("LumA", "LumB", "Her2", "Basal")

subtype_colors <- c(
  LumA   = "#2166AC",
  LumB   = "#92C5DE",
  Her2   = "#D6604D",
  Basal  = "#1A1A1A",
  Normal = "#4DAC26"
)

group_colors <- c(
  Tumor = "#D73027",
  NAT   = "#FC8D59",
  GTEx  = "#4575B4"
)

# ---------------------------------------------------------------------------
# Stage gate
# ---------------------------------------------------------------------------

#' Stop with a clear error if Stage 01 output is missing
#'
#' @return Invisible NULL.
#' @examples
#' assert_stage01_complete()
assert_stage01_complete <- function() {
  message("── Stage 02 | stage gate ────────────────────────────────────────────────")
  path <- here::here("data", "processed", "cohort.rds")
  if (!file.exists(path))
    stop("Stage 01 output missing: ", path, "\nRun R/01_preprocessing.R first.")
  message("  cohort.rds found.")
  message("── Stage 02 | gate passed ───────────────────────────────────────────────")
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

#' Load cohort.rds and return its three components
#'
#' @return Named list: \code{$expr} (matrix), \code{$meta} (data.frame),
#'   \code{$contrasts} (character vector).
#' @examples
#' cohort <- load_cohort()
load_cohort <- function() {
  message("── Stage 02 | loading cohort ────────────────────────────────────────────")
  cohort <- readRDS(here::here("data", "processed", "cohort.rds"))
  message("  $expr      : ", nrow(cohort$expr), " genes × ", ncol(cohort$expr), " samples")
  message("  $meta      : ", nrow(cohort$meta), " rows")
  message("── Stage 02 | cohort loaded ─────────────────────────────────────────────")
  cohort
}

# ---------------------------------------------------------------------------
# Highly variable genes
# ---------------------------------------------------------------------------

#' Select top N most variable genes across all samples
#'
#' Variance is computed per gene across all 1,109 samples (Tumor + NAT + GTEx).
#' Using all samples here maximises separation between groups in UMAP space,
#' which is the goal of the visualisation stage.
#'
#' @param expr Numeric matrix. HGNC symbols × sample IDs, log2(TPM+0.001).
#' @param n    Integer. Number of top-variance genes to retain. Default 5000.
#' @return Numeric matrix subset to \code{n} rows.
#' @examples
#' hvg <- select_hvg(cohort$expr)
select_hvg <- function(expr, n = N_HVG) {
  message("── Stage 02 | selecting HVGs ────────────────────────────────────────────")
  gene_var <- matrixStats::rowVars(expr)
  top_idx  <- order(gene_var, decreasing = TRUE)[seq_len(min(n, nrow(expr)))]
  hvg      <- expr[top_idx, , drop = FALSE]
  message("  HVGs selected: ", nrow(hvg))
  message("── Stage 02 | HVG selection complete ───────────────────────────────────")
  hvg
}

# ---------------------------------------------------------------------------
# PCA
# ---------------------------------------------------------------------------

#' Run PCA on transposed HVG matrix and return top components
#'
#' Samples are rows, genes are columns — the standard orientation for
#' \code{prcomp}. Returns the rotated scores matrix (samples × PCs).
#'
#' @param hvg Numeric matrix. Genes × samples (HVG-filtered).
#' @param n   Integer. Number of PCs to retain. Default 50.
#' @return Numeric matrix. Samples × PCs.
#' @examples
#' pcs <- run_pca(hvg)
run_pca <- function(hvg, n = N_PCA) {
  message("── Stage 02 | running PCA ───────────────────────────────────────────────")
  pca  <- prcomp(t(hvg), center = TRUE, scale. = FALSE, rank. = n)
  pcs  <- pca$x
  var_exp <- summary(pca)$importance[2, seq_len(min(5L, n))]
  message("  variance explained (PC1–5): ",
          paste0(round(var_exp * 100, 1), "%", collapse = ", "))
  message("── Stage 02 | PCA complete ──────────────────────────────────────────────")
  pcs
}

# ---------------------------------------------------------------------------
# UMAP
# ---------------------------------------------------------------------------

#' Run UMAP on PCA scores
#'
#' Uses \code{uwot::umap} with reproducible seed. Returns a two-column
#' data.frame aligned to the rows of \code{pcs}.
#'
#' @param pcs Numeric matrix. Samples × PCs from \code{run_pca()}.
#' @return data.frame with columns UMAP1, UMAP2. Row names = sample IDs.
#' @examples
#' umap_coords <- run_umap(pcs)
run_umap <- function(pcs) {
  message("── Stage 02 | running UMAP ──────────────────────────────────────────────")
  coords <- uwot::umap(pcs, n_neighbors = 30L, min_dist = 0.3,
                       n_components = 2L, seed = 42L)
  df <- data.frame(UMAP1 = coords[, 1], UMAP2 = coords[, 2],
                   row.names = rownames(pcs))
  message("  UMAP complete: ", nrow(df), " samples embedded")
  message("── Stage 02 | UMAP complete ─────────────────────────────────────────────")
  df
}

# ---------------------------------------------------------------------------
# Sub-subcluster detection
# ---------------------------------------------------------------------------

#' Detect within-subtype sub-subclusters via k-means on UMAP coordinates
#'
#' For each PAM50 subtype in scope, runs k-means (k = \code{k}) on the
#' two UMAP dimensions. Sub-subcluster labels are stored as
#' \code{"{Subtype}_C{k}"} (e.g. \code{"LumA_C1"}). Control samples
#' (NAT, GTEx) receive \code{NA}.
#'
#' @param umap_df data.frame. Columns: UMAP1, UMAP2, group, PAM50.
#'   Row names = sample IDs.
#' @param k Integer. Number of clusters per subtype. Default 2.
#' @return Character vector of sub-subcluster labels, aligned to rows of
#'   \code{umap_df}.
#' @examples
#' umap_df$subcluster <- detect_subclusters(umap_df)
detect_subclusters <- function(umap_df, k = K_SUBCLUSTER) {
  message("── Stage 02 | detecting sub-subclusters (k = ", k, ") ───────────────────")
  labels <- rep(NA_character_, nrow(umap_df))

  for (subtype in PAM50_LEVELS) {
    idx <- which(umap_df$PAM50 == subtype & !is.na(umap_df$PAM50))
    if (length(idx) < k) {
      message("  ", subtype, ": too few samples (", length(idx), ") — skipping")
      next
    }
    km <- kmeans(umap_df[idx, c("UMAP1", "UMAP2")],
                 centers = k, nstart = 25L, iter.max = 100L)
    labels[idx] <- paste0(subtype, "_C", km$cluster)
    message("  ", subtype, ": ", paste(table(km$cluster), collapse = " | "),
            " samples per cluster")
  }

  message("── Stage 02 | sub-subcluster detection complete ─────────────────────────")
  labels
}

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

#' UMAP coloured by sample group (Tumor / NAT / GTEx)
#'
#' @param umap_df data.frame. Columns: UMAP1, UMAP2, group, PAM50.
#' @return ggplot object.
#' @examples
#' p <- plot_umap_group(umap_df)
plot_umap_group <- function(umap_df) {
  caption <- paste0("UMAP (", N_HVG, " HVGs, ", N_PCA, " PCs, n_neighbors=30). ",
                    "n = ", nrow(umap_df), " samples.")
  ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = group)) +
    geom_point(size = 0.6, alpha = 0.7) +
    scale_colour_manual(values = group_colors, name = "Group") +
    labs(title = "UMAP — Sample Groups",
         subtitle = "Tumor (PAM50-classified) · NAT · GTEx Healthy Breast",
         x = "UMAP 1", y = "UMAP 2",
         caption = caption) +
    theme_bw(base_size = 12) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

#' UMAP coloured by PAM50 subtype (Tumor samples only highlighted)
#'
#' @param umap_df data.frame. Columns: UMAP1, UMAP2, group, PAM50.
#' @return ggplot object.
#' @examples
#' p <- plot_umap_subtype(umap_df)
plot_umap_subtype <- function(umap_df) {
  # Controls plotted in grey underneath; Tumor subtypes on top in subtype_colors
  controls <- dplyr::filter(umap_df, group != "Tumor")
  tumors   <- dplyr::filter(umap_df, group == "Tumor")

  caption <- paste0("PAM50 subtypes: LumA n=", sum(tumors$PAM50 == "LumA", na.rm=TRUE),
                    ", LumB n=", sum(tumors$PAM50 == "LumB", na.rm=TRUE),
                    ", Her2 n=", sum(tumors$PAM50 == "Her2", na.rm=TRUE),
                    ", Basal n=", sum(tumors$PAM50 == "Basal", na.rm=TRUE), ".")

  ggplot() +
    geom_point(data = controls, aes(x = UMAP1, y = UMAP2),
               colour = "grey80", size = 0.5, alpha = 0.5) +
    geom_point(data = tumors, aes(x = UMAP1, y = UMAP2, colour = PAM50),
               size = 0.8, alpha = 0.8) +
    scale_colour_manual(values = subtype_colors, name = "PAM50") +
    labs(title = "UMAP — PAM50 Subtypes",
         subtitle = "Controls (NAT, GTEx) shown in grey",
         x = "UMAP 1", y = "UMAP 2",
         caption = caption) +
    theme_bw(base_size = 12) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

#' UMAP coloured by within-subtype sub-subcluster for one PAM50 subtype
#'
#' @param umap_df  data.frame. Columns: UMAP1, UMAP2, group, PAM50, subcluster.
#' @param subtype  Character. One of PAM50_LEVELS.
#' @return ggplot object.
#' @examples
#' p <- plot_umap_subcluster(umap_df, "LumA")
plot_umap_subcluster <- function(umap_df, subtype) {
  sub_df    <- dplyr::filter(umap_df, PAM50 == subtype)
  other_df  <- dplyr::filter(umap_df, PAM50 != subtype | is.na(PAM50))
  n_clusters <- length(unique(stats::na.omit(sub_df$subcluster)))
  cluster_colors <- setNames(
    scales::hue_pal()(n_clusters),
    sort(unique(stats::na.omit(sub_df$subcluster)))
  )

  ggplot() +
    geom_point(data = other_df, aes(x = UMAP1, y = UMAP2),
               colour = "grey85", size = 0.4, alpha = 0.4) +
    geom_point(data = sub_df, aes(x = UMAP1, y = UMAP2, colour = subcluster),
               size = 1.0, alpha = 0.9) +
    scale_colour_manual(values = cluster_colors, name = "Sub-cluster",
                        na.value = "grey85") +
    labs(title = paste0("UMAP — Within-Subtype Sub-Clusters: ", subtype),
         subtitle = paste0(n_clusters, " clusters (k-means, k=", K_SUBCLUSTER,
                           "), n=", nrow(sub_df), " samples"),
         x = "UMAP 1", y = "UMAP 2") +
    theme_bw(base_size = 12) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

# ---------------------------------------------------------------------------
# Save helper
# ---------------------------------------------------------------------------

#' Print a ggplot to screen and save as .png (300 dpi) and .pdf
#'
#' @param p        ggplot object.
#' @param stem     Character. File stem (no extension), relative to
#'   \code{results/figures/}.
#' @param width    Numeric. Width in inches. Default 10.
#' @param height   Numeric. Height in inches. Default 8.
#' @return Invisible NULL.
#' @examples
#' save_plot(p, "umap_by_group")
save_plot <- function(p, stem, width = 10, height = 8) {
  print(p)   # display to screen during interactive or Rscript runs

  dir.create(here::here("results", "figures"), recursive = TRUE,
             showWarnings = FALSE)
  base <- here::here("results", "figures", stem)

  ggplot2::ggsave(paste0(base, ".png"), plot = p,
                  width = width, height = height, dpi = 300)
  ggplot2::ggsave(paste0(base, ".pdf"), plot = p,
                  width = width, height = height)

  message("  saved: ", stem, ".png + .pdf")
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Save cohort output
# ---------------------------------------------------------------------------

#' Save umap_clusters.rds to data/processed/
#'
#' @param umap_df      data.frame. Full UMAP table with subcluster column.
#' @param n_per_subtype Named integer vector.
#' @return Invisible character path.
#' @examples
#' save_clusters(umap_df, n_per_subtype)
save_clusters <- function(umap_df, n_per_subtype) {
  message("── Stage 02 | saving umap_clusters.rds ──────────────────────────────────")
  out <- here::here("data", "processed", "umap_clusters.rds")
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(umap = umap_df, n_per_subtype = n_per_subtype), file = out)
  message("  saved: ", out)
  message("  $umap          : ", nrow(umap_df), " rows")
  message("  $n_per_subtype :")
  print(n_per_subtype)
  message("── Stage 02 | save complete ──────────────────────────────────────────────")
  invisible(out)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

assert_stage01_complete()

cohort <- load_cohort()
expr   <- cohort$expr
meta   <- cohort$meta

# Check matrixStats is available (used in select_hvg)
if (!requireNamespace("matrixStats", quietly = TRUE))
  stop("Package 'matrixStats' required. Install with: install.packages('matrixStats')")

hvg  <- select_hvg(expr)
pcs  <- run_pca(hvg)
umap_coords <- run_umap(pcs)

# Join UMAP coords with metadata
umap_df <- umap_coords |>
  tibble::rownames_to_column("sample_id") |>
  dplyr::left_join(meta, by = "sample_id")

# Detect within-subtype sub-subclusters
umap_df$subcluster <- detect_subclusters(umap_df)

# Sub-subcluster counts per subtype
n_per_subtype <- table(umap_df$subcluster) |>
  as.integer() |>
  setNames(names(table(umap_df$subcluster)))

# Plots -------------------------------------------------------------------
message("── Stage 02 | generating plots ──────────────────────────────────────────")

p_group   <- plot_umap_group(umap_df)
save_plot(p_group, "umap_by_group")

p_subtype <- plot_umap_subtype(umap_df)
save_plot(p_subtype, "umap_by_subtype")

for (st in PAM50_LEVELS) {
  p_sub <- plot_umap_subcluster(umap_df, st)
  save_plot(p_sub, paste0("umap_subclusters_", st))
}

message("── Stage 02 | all plots saved ───────────────────────────────────────────")

save_clusters(umap_df, n_per_subtype)

message("✅ Stage 02 complete. Proceed to R/03_deg.R")
