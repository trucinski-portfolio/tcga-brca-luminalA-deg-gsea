# R/06_wgcna.R
# Stage 06: WGCNA co-expression network analysis on Tumor samples.
#           Detects modules, correlates with PAM50 traits, saves heatmap.
#
# Inputs:  data/processed/cohort.rds
# Outputs: data/processed/wgcna_modules.rds
#            $module_labels     — named integer vector (gene → module number)
#            $module_colors     — named character vector (gene → colour label)
#            $MEs               — data.frame: module eigengenes (samples × modules)
#            $module_trait_cor  — matrix: Pearson r (modules × PAM50 traits)
#            $module_trait_pval — matrix: p-values
#            $power             — integer: selected soft-threshold power
#            $gene_count        — named integer: genes per module
#          results/figures/wgcna_soft_threshold.{png,pdf}
#          results/figures/wgcna_module_trait_heatmap.{png,pdf}
#
# Note: WGCNA on 5000 genes × 817 samples takes ~15-30 min.
#       Run in a fresh R session with ≥8 GB RAM available.

library(here)
library(dplyr)
library(WGCNA)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(matrixStats)

# WGCNA replaces stats::cor with its own faster implementation — allow this.
options(stringsAsFactors = FALSE)
WGCNA::enableWGCNAThreads()  # use all available cores

set.seed(42)

N_HVG          <- 5000L   # top MAD genes fed to WGCNA
NETWORK_TYPE   <- "signed"
MIN_MODULE     <- 30L     # minimum genes per module
MERGE_CUT      <- 0.25    # dendrogram cut height for module merging
POWER_RANGE    <- 1:30    # candidate soft-threshold powers to test
R2_TARGET      <- 0.85    # minimum scale-free topology R²

PAM50_LEVELS <- c("LumA", "LumB", "Her2", "Basal")

subtype_colors <- c(
  LumA  = "#2166AC",
  LumB  = "#92C5DE",
  Her2  = "#D6604D",
  Basal = "#1A1A1A"
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
  message("── Stage 06 | stage gate ────────────────────────────────────────────────")
  path <- here::here("data", "processed", "cohort.rds")
  if (!file.exists(path))
    stop("cohort.rds missing — run R/01_preprocessing.R first.")
  message("  cohort.rds found.")
  message("── Stage 06 | gate passed ───────────────────────────────────────────────")
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Load and prepare
# ---------------------------------------------------------------------------

#' Load cohort and return Tumor-only expression in WGCNA orientation
#'
#' Subsets to Tumor samples, selects top \code{n} genes by MAD, and
#' transposes to samples × genes (WGCNA convention). Low-variance genes
#' are excluded to reduce noise in the adjacency matrix.
#'
#' @param n Integer. Number of top-MAD genes to retain. Default 5000.
#' @return Named list: \code{$expr_t} (matrix, samples × genes),
#'   \code{$meta} (Tumor-only data.frame).
#' @examples
#' dat <- prepare_expr()
prepare_expr <- function(n = N_HVG) {
  message("── Stage 06 | preparing expression ──────────────────────────────────────")
  cohort <- readRDS(here::here("data", "processed", "cohort.rds"))
  expr   <- cohort$expr
  meta   <- cohort$meta

  # Tumor samples only
  tumor_meta <- meta[meta$group == "Tumor" & !is.na(meta$PAM50), ]
  tumor_expr <- expr[, tumor_meta$sample_id, drop = FALSE]
  message("  Tumor samples: ", ncol(tumor_expr))

  # Top n genes by MAD
  mad_vals  <- matrixStats::rowMads(tumor_expr)
  top_genes <- order(mad_vals, decreasing = TRUE)[seq_len(min(n, nrow(tumor_expr)))]
  tumor_expr <- tumor_expr[top_genes, , drop = FALSE]
  message("  Genes selected (top MAD): ", nrow(tumor_expr))

  # Check for genes with zero variance (WGCNA will error on these)
  gsg <- WGCNA::goodSamplesGenes(t(tumor_expr), verbose = 0)
  if (!gsg$allOK) {
    tumor_expr <- tumor_expr[gsg$goodGenes, gsg$goodSamples, drop = FALSE]
    message("  Removed low-quality genes/samples. Remaining: ",
            nrow(tumor_expr), " genes × ", ncol(tumor_expr), " samples")
  }

  message("  Final matrix: ", nrow(tumor_expr), " genes × ", ncol(tumor_expr),
          " samples")
  message("── Stage 06 | expression prepared ──────────────────────────────────────")
  list(expr_t = t(tumor_expr), meta = tumor_meta)
}

# ---------------------------------------------------------------------------
# Soft threshold
# ---------------------------------------------------------------------------

#' Select soft-threshold power and produce diagnostic plot
#'
#' Tests powers \code{POWER_RANGE} and selects the lowest power where the
#' scale-free topology fit R² ≥ \code{R2_TARGET}. Falls back to the power
#' with the highest R² if the target is not met. Prints the plot to screen
#' and saves as .png + .pdf.
#'
#' @param expr_t Numeric matrix. Samples × genes.
#' @return Integer. Selected soft-threshold power.
#' @examples
#' power <- pick_soft_threshold(expr_t)
pick_soft_threshold <- function(expr_t) {
  message("── Stage 06 | picking soft threshold ────────────────────────────────────")
  message("  testing powers: ", min(POWER_RANGE), "–", max(POWER_RANGE),
          " (target R² ≥ ", R2_TARGET, ") ...")

  sft <- WGCNA::pickSoftThreshold(
    expr_t,
    powerVector  = POWER_RANGE,
    networkType  = NETWORK_TYPE,
    verbose      = 0
  )

  fit_df <- data.frame(
    power = sft$fitIndices$Power,
    r2    = -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
    k     = sft$fitIndices$mean.k.
  )

  # Select power
  meets_target <- fit_df$power[fit_df$r2 >= R2_TARGET]
  power <- if (length(meets_target) > 0) min(meets_target) else
    fit_df$power[which.max(fit_df$r2)]

  message("  selected power: ", power,
          " (R² = ", round(fit_df$r2[fit_df$power == power], 3), ")")

  # Plot
  p1 <- ggplot(fit_df, aes(x = power, y = r2)) +
    geom_line(colour = "steelblue") +
    geom_point(colour = "steelblue", size = 2) +
    geom_hline(yintercept = R2_TARGET, linetype = "dashed", colour = "red") +
    geom_vline(xintercept = power, linetype = "dotted", colour = "darkgreen") +
    annotate("text", x = power + 0.5, y = 0.1,
             label = paste0("power = ", power), hjust = 0, colour = "darkgreen") +
    scale_x_continuous(breaks = POWER_RANGE[seq(1, length(POWER_RANGE), by = 2)]) +
    labs(title  = "WGCNA Soft-Threshold Power Selection",
         subtitle = paste0(NETWORK_TYPE, " network | target R² = ", R2_TARGET),
         x = "Soft-Threshold Power",
         y = "Scale-Free Topology R²") +
    theme_bw(base_size = 12)

  print(p1)
  dir.create(here::here("results", "figures"), recursive = TRUE,
             showWarnings = FALSE)
  ggplot2::ggsave(here::here("results", "figures", "wgcna_soft_threshold.png"),
                  plot = p1, width = 8, height = 5, dpi = 300)
  ggplot2::ggsave(here::here("results", "figures", "wgcna_soft_threshold.pdf"),
                  plot = p1, width = 8, height = 5)
  message("  saved: wgcna_soft_threshold.png + .pdf")
  message("── Stage 06 | soft threshold selected ───────────────────────────────────")
  power
}

# ---------------------------------------------------------------------------
# Build network
# ---------------------------------------------------------------------------

#' Run blockwiseModules to detect co-expression modules
#'
#' Uses a signed network with the selected soft-threshold power.
#' \code{blockwiseModules} handles large datasets by processing in blocks;
#' with 5,000 genes and \code{maxBlockSize = 6000} the full matrix fits in
#' one block. The grey module (label 0) contains unassigned genes.
#'
#' @param expr_t Numeric matrix. Samples × genes.
#' @param power  Integer. Soft-threshold power from \code{pick_soft_threshold()}.
#' @return WGCNA module object (list) from \code{blockwiseModules()}.
#' @examples
#' net <- build_network(expr_t, power)
build_network <- function(expr_t, power) {
  message("── Stage 06 | building network (power = ", power, ") ─────────────────────")
  message("  this step takes ~15–30 min ...")

  net <- WGCNA::blockwiseModules(
    expr_t,
    power             = power,
    networkType       = NETWORK_TYPE,
    TOMType           = NETWORK_TYPE,
    minModuleSize     = MIN_MODULE,
    mergeCutHeight    = MERGE_CUT,
    maxBlockSize      = 6000L,
    numericLabels     = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs          = FALSE,
    verbose           = 2
  )

  n_modules <- length(unique(net$colors)) - 1L  # exclude grey (0)
  message("  modules detected (excl. grey): ", n_modules)
  message("  genes per module (top 10):")
  mod_sizes <- sort(table(net$colors), decreasing = TRUE)
  print(head(mod_sizes, 10))

  message("── Stage 06 | network built ─────────────────────────────────────────────")
  net
}

# ---------------------------------------------------------------------------
# Module–trait correlation
# ---------------------------------------------------------------------------

#' Correlate module eigengenes with PAM50 binary trait indicators
#'
#' Creates one binary trait column per PAM50 subtype (1 = member, 0 = other)
#' and correlates each module's eigengene with each trait using Pearson
#' correlation. Returns correlation and p-value matrices.
#'
#' @param net    blockwiseModules result from \code{build_network()}.
#' @param meta   data.frame. Tumor-only metadata with PAM50 column.
#' @return Named list: \code{$cor} (modules × traits), \code{$pval},
#'   \code{$MEs} (module eigengenes data.frame).
#' @examples
#' mt <- correlate_modules_traits(net, meta)
correlate_modules_traits <- function(net, meta) {
  message("── Stage 06 | module–trait correlations ─────────────────────────────────")

  MEs <- WGCNA::orderMEs(net$MEs)

  # Binary trait matrix: samples × PAM50 subtypes
  traits <- sapply(PAM50_LEVELS, function(st) {
    as.integer(meta$PAM50 == st)
  })
  rownames(traits) <- meta$sample_id

  # Pearson correlation + p-values
  cor_res  <- WGCNA::cor(MEs, traits, use = "p")
  pval_res <- WGCNA::corPvalueStudent(cor_res, nrow(MEs))

  n_sig <- sum(pval_res < 0.05, na.rm = TRUE)
  message("  modules × traits: ", nrow(cor_res), " × ", ncol(cor_res))
  message("  significant associations (p<0.05): ", n_sig)

  message("── Stage 06 | correlations complete ─────────────────────────────────────")
  list(cor = cor_res, pval = pval_res, MEs = MEs)
}

# ---------------------------------------------------------------------------
# Heatmap
# ---------------------------------------------------------------------------

#' Draw and save module–trait correlation heatmap
#'
#' Rows = modules, columns = PAM50 subtypes. Cell colour encodes Pearson r
#' (blue = negative, red = positive). Significance asterisks are overlaid:
#' *** p<0.001, ** p<0.01, * p<0.05. Printed to screen and saved as
#' .png + .pdf.
#'
#' @param cor_mat  Matrix. Module × trait Pearson correlations.
#' @param pval_mat Matrix. Corresponding p-values.
#' @return Invisible NULL.
#' @examples
#' plot_module_trait_heatmap(mt$cor, mt$pval)
plot_module_trait_heatmap <- function(cor_mat, pval_mat) {
  message("── Stage 06 | plotting module–trait heatmap ─────────────────────────────")

  sig_text <- matrix("", nrow = nrow(pval_mat), ncol = ncol(pval_mat))
  sig_text[pval_mat < 0.001] <- "***"
  sig_text[pval_mat >= 0.001 & pval_mat < 0.01]  <- "**"
  sig_text[pval_mat >= 0.01  & pval_mat < 0.05]  <- "*"

  col_fun <- circlize::colorRamp2(c(-1, 0, 1),
                                   c("#4575B4", "white", "#D73027"))

  col_anno <- ComplexHeatmap::HeatmapAnnotation(
    Subtype = PAM50_LEVELS,
    col     = list(Subtype = subtype_colors),
    show_annotation_name = FALSE
  )

  ht <- ComplexHeatmap::Heatmap(
    cor_mat,
    name               = "Pearson r",
    col                = col_fun,
    top_annotation     = col_anno,
    cell_fun           = function(j, i, x, y, width, height, fill) {
      grid::grid.text(sig_text[i, j], x, y,
                      gp = grid::gpar(fontsize = 8, col = "grey20"))
    },
    column_names_gp    = grid::gpar(col = subtype_colors[PAM50_LEVELS],
                                     fontsize = 11, fontface = "bold"),
    row_names_gp       = grid::gpar(fontsize = 8),
    cluster_columns    = FALSE,
    show_row_dend      = TRUE,
    row_title          = "Co-expression Modules",
    column_title       = "PAM50 Subtype",
    heatmap_legend_param = list(
      title    = "Pearson r",
      at       = c(-1, -0.5, 0, 0.5, 1),
      labels   = c("-1", "-0.5", "0", "0.5", "1")
    )
  )

  draw(ht)   # print to screen

  out_base <- here::here("results", "figures", "wgcna_module_trait_heatmap")
  png(paste0(out_base, ".png"), width = 2400, height = 2000, res = 300)
  draw(ht)
  dev.off()
  pdf(paste0(out_base, ".pdf"), width = 10, height = 9)
  draw(ht)
  dev.off()
  message("  saved: wgcna_module_trait_heatmap.png + .pdf")
  message("── Stage 06 | heatmap saved ──────────────────────────────────────────────")
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

#' Save wgcna_modules.rds to data/processed/
#'
#' @param net     blockwiseModules result.
#' @param mt      Module-trait list from \code{correlate_modules_traits()}.
#' @param power   Integer. Selected soft-threshold power.
#' @return Invisible character path.
#' @examples
#' save_modules(net, mt, power)
save_modules <- function(net, mt, power) {
  message("── Stage 06 | saving wgcna_modules.rds ──────────────────────────────────")
  out <- here::here("data", "processed", "wgcna_modules.rds")
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

  saveRDS(list(
    module_labels     = net$colors,
    module_colors     = WGCNA::labels2colors(net$colors),
    MEs               = mt$MEs,
    module_trait_cor  = mt$cor,
    module_trait_pval = mt$pval,
    power             = power,
    gene_count        = table(net$colors)
  ), file = out)

  message("  saved: ", out)
  message("  modules (excl. grey): ",
          length(unique(net$colors[net$colors != 0])))
  message("── Stage 06 | save complete ──────────────────────────────────────────────")
  invisible(out)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

assert_stage01_complete()

dat    <- prepare_expr()
expr_t <- dat$expr_t
meta   <- dat$meta

power  <- pick_soft_threshold(expr_t)
net    <- build_network(expr_t, power)
mt     <- correlate_modules_traits(net, meta)

plot_module_trait_heatmap(mt$cor, mt$pval)
save_modules(net, mt, power)

message("✅ Stage 06 complete. Proceed to R/07_survival.R")
