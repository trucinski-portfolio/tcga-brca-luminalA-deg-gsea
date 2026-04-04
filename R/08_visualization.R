# R/08_visualization.R
# Stage 08: Publication-quality figures.
#           Assembles all major results into a consistent figure set.
#
# Inputs:  data/processed/cohort.rds
#          data/processed/umap_clusters.rds
#          data/processed/wgcna_modules.rds
#          results/tables/deg/*.csv           (18 files)
#          results/tables/signatures/*.csv    (4 files)
#          results/tables/fgsea/*.csv         (4 files)
#          results/tables/survival/*.csv      (4 files)
#
# Outputs (results/figures/):
#   fig1_umap_overview.{png,pdf}
#   fig2_deg_counts.{png,pdf}
#   fig3_volcano_{LumA,LumB,Her2,Basal}_vs_GTEx.{png,pdf}
#   fig4_deg_heatmap.{png,pdf}
#   fig5_signature_heatmap.{png,pdf}
#   fig6_fgsea_dotplot.{png,pdf}
#   fig7_cox_forest.{png,pdf}

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(matrixStats)

set.seed(42)

FIG_DIR    <- here::here("results", "figures")
PAM50_LEVELS <- c("LumA", "LumB", "Her2", "Basal")

subtype_colors <- c(
  LumA  = "#2166AC",
  LumB  = "#92C5DE",
  Her2  = "#D6604D",
  Basal = "#1A1A1A"
)

group_colors <- c(
  Tumor = "#E41A1C",
  NAT   = "#FF7F00",
  GTEx  = "#4DAC26"
)

# ---------------------------------------------------------------------------
# Stage gate
# ---------------------------------------------------------------------------

#' Verify all Stage 08 inputs exist before running any plotting code
#'
#' @return Invisible NULL. Stops with a clear error if any input is missing.
#' @examples
#' assert_all_stages_complete()
assert_all_stages_complete <- function() {
  message("── Stage 08 | stage gate ────────────────────────────────────────────────")

  required_rds <- c(
    "data/processed/cohort.rds",
    "data/processed/umap_clusters.rds",
    "data/processed/wgcna_modules.rds"
  )
  missing_rds <- required_rds[!file.exists(here::here(required_rds))]
  if (length(missing_rds) > 0)
    stop("Missing .rds files: ", paste(missing_rds, collapse = ", "))

  deg_n  <- length(list.files(here::here("results", "tables", "deg"),
                               pattern = "[.]csv$"))
  sig_n  <- length(list.files(here::here("results", "tables", "signatures"),
                               pattern = "[.]csv$"))
  gsea_n <- length(list.files(here::here("results", "tables", "fgsea"),
                               pattern = "[.]csv$"))
  surv_n <- length(list.files(here::here("results", "tables", "survival"),
                               pattern = "[.]csv$"))

  checks <- c("DEG CSVs (need 18)"       = deg_n  == 18L,
               "Signature CSVs (need 4)"  = sig_n  == 4L,
               "fgsea CSVs (need 4)"      = gsea_n == 4L,
               "Survival CSVs (need 4)"   = surv_n == 4L)
  failed <- names(checks)[!checks]
  if (length(failed) > 0)
    stop("Stage gate failed: ", paste(failed, collapse = "; "))

  message("  all inputs present.")
  message("── Stage 08 | gate passed ───────────────────────────────────────────────")
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Save a ggplot as both .png (300 dpi) and .pdf, printing to screen first
#'
#' @param p        ggplot object.
#' @param stem     Character. Filename without extension (placed in FIG_DIR).
#' @param width    Numeric. Width in inches. Default 10.
#' @param height   Numeric. Height in inches. Default 8.
#' @return Invisible NULL.
#' @examples
#' save_fig(p, "fig1_umap_overview")
save_fig <- function(p, stem, width = 10, height = 8) {
  print(p)
  ggsave(file.path(FIG_DIR, paste0(stem, ".png")),
         plot = p, width = width, height = height, dpi = 300)
  ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")),
         plot = p, width = width, height = height)
  message("  saved: ", stem, ".png + .pdf")
  invisible(NULL)
}

#' Save a ComplexHeatmap as both .png and .pdf, drawing to screen first
#'
#' @param ht       Heatmap or HeatmapList object.
#' @param stem     Character. Filename without extension.
#' @param width_px Integer. PNG width in pixels. Default 2400.
#' @param height_px Integer. PNG height in pixels. Default 2000.
#' @param width_in Numeric. PDF width in inches. Default 10.
#' @param height_in Numeric. PDF height in inches. Default 8.
#' @return Invisible NULL.
#' @examples
#' save_heatmap(ht, "fig4_deg_heatmap")
save_heatmap <- function(ht, stem,
                         width_px = 2400, height_px = 2000,
                         width_in = 10,  height_in = 8) {
  draw(ht)
  png(file.path(FIG_DIR, paste0(stem, ".png")),
      width = width_px, height = height_px, res = 300)
  draw(ht)
  dev.off()
  pdf(file.path(FIG_DIR, paste0(stem, ".pdf")),
      width = width_in, height = height_in)
  draw(ht)
  dev.off()
  message("  saved: ", stem, ".png + .pdf")
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Fig 1 — UMAP overview (group + subtype panels)
# ---------------------------------------------------------------------------

#' Figure 1: two-panel UMAP — left colored by group, right by PAM50 subtype
#'
#' @param umap_df  data.frame. Columns: sample_id, UMAP1, UMAP2, subtype, group.
#' @return ggplot object.
#' @examples
#' p <- plot_umap_overview(umap_df)
plot_umap_overview <- function(umap_df) {
  message("  Fig 1: UMAP overview ...")

  df <- umap_df |>
    dplyr::mutate(
      subtype_label = dplyr::case_when(
        !is.na(PAM50) & PAM50 %in% PAM50_LEVELS ~ PAM50,
        TRUE                                     ~ "Control"
      )
    )

  # Panel A — by group
  pa <- ggplot(df, aes(UMAP1, UMAP2, colour = group)) +
    geom_point(size = 0.7, alpha = 0.7) +
    scale_colour_manual(values = group_colors, name = "Group") +
    labs(title = "A   Sample groups", x = "UMAP 1", y = "UMAP 2") +
    theme_bw(base_size = 11) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"))

  # Panel B — by PAM50 subtype (controls grey)
  subtype_pal <- c(subtype_colors, Control = "grey80")
  pb <- ggplot(df |> dplyr::arrange(subtype_label == "Control"),
               aes(UMAP1, UMAP2, colour = subtype_label)) +
    geom_point(size = 0.7, alpha = 0.7) +
    scale_colour_manual(values = subtype_pal, name = "PAM50") +
    labs(title = "B   PAM50 subtype", x = "UMAP 1", y = "UMAP 2") +
    theme_bw(base_size = 11) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"))

  p <- cowplot::plot_grid(pa, pb, nrow = 1, rel_widths = c(1, 1))
  p
}

# ---------------------------------------------------------------------------
# Fig 2 — DEG count bar chart
# ---------------------------------------------------------------------------

#' Figure 2: bar chart of significant DEG counts per contrast (FDR < 0.05, |logFC| >= 1)
#'
#' Reads all 18 DEG CSVs, applies thresholds, and plots a horizontal bar chart
#' faceted by comparison type (vs GTEx, vs NAT, subtype vs subtype).
#'
#' @param deg_dir Character. Path to results/tables/deg/.
#' @return ggplot object.
#' @examples
#' p <- plot_deg_counts(here::here("results", "tables", "deg"))
plot_deg_counts <- function(deg_dir) {
  message("  Fig 2: DEG count bar chart ...")

  files <- list.files(deg_dir, pattern = "[.]csv$", full.names = TRUE)

  counts <- lapply(files, function(f) {
    d        <- read.csv(f, stringsAsFactors = FALSE)
    contrast <- sub("[.]csv$", "", basename(f))
    n_up     <- sum(d$adj.P.Val < 0.05 & d$logFC >  1, na.rm = TRUE)
    n_dn     <- sum(d$adj.P.Val < 0.05 & d$logFC < -1, na.rm = TRUE)
    data.frame(contrast = contrast, up = n_up, down = -n_dn,
               stringsAsFactors = FALSE)
  }) |> dplyr::bind_rows()

  # Classify contrast type for faceting
  counts <- counts |>
    dplyr::mutate(
      type = dplyr::case_when(
        grepl("_vs_GTEx$",   contrast) ~ "vs GTEx",
        grepl("_vs_NAT$",    contrast) ~ "vs NAT",
        grepl("NAT_vs_GTEx", contrast) ~ "NAT vs GTEx",
        grepl("_C[12]",      contrast) ~ "Within-subtype",
        TRUE                           ~ "Subtype vs Subtype"
      )
    )

  long <- counts |>
    tidyr::pivot_longer(c(up, down), names_to = "direction", values_to = "n") |>
    dplyr::mutate(
      label     = ifelse(direction == "up", "Up", "Down"),
      fill_col  = ifelse(direction == "up", "#D6604D", "#4393C3")
    )

  p <- ggplot(long, aes(x = reorder(contrast, abs(n)), y = n, fill = label)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = c(Up = "#D6604D", Down = "#4393C3"),
                      name = "Direction") +
    facet_wrap(~ type, scales = "free_y", ncol = 1) +
    coord_flip() +
    labs(
      title   = "Significant DEGs per contrast",
      subtitle = "FDR < 0.05, |log\u2082FC| \u2265 1",
      x = NULL, y = "Gene count"
    ) +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey92"),
          plot.title = element_text(face = "bold"))
  p
}

# ---------------------------------------------------------------------------
# Fig 3 — Volcano plots (subtype vs GTEx)
# ---------------------------------------------------------------------------

#' Figure 3: volcano plot for one subtype vs GTEx contrast
#'
#' Labels the top 15 genes by significance. Highlights genes in the
#' subtype's ElasticNet signature.
#'
#' @param deg_path Character. Path to the DEG CSV.
#' @param sig_path Character. Path to the signature CSV.
#' @param subtype  Character. Subtype name for title.
#' @return ggplot object.
#' @examples
#' p <- plot_volcano(deg_path, sig_path, "LumA")
plot_volcano <- function(deg_path, sig_path, subtype) {
  deg <- read.csv(deg_path, stringsAsFactors = FALSE) |>
    dplyr::filter(!is.na(adj.P.Val), !is.na(logFC))
  sig_genes <- read.csv(sig_path,  stringsAsFactors = FALSE)$gene

  deg <- deg |>
    dplyr::mutate(
      neg_log10_p = -log10(adj.P.Val + 1e-300),
      status = dplyr::case_when(
        adj.P.Val < 0.05 & logFC >  1 & gene %in% sig_genes ~ "Signature (up)",
        adj.P.Val < 0.05 & logFC < -1 & gene %in% sig_genes ~ "Signature (down)",
        adj.P.Val < 0.05 & logFC >  1                        ~ "Up",
        adj.P.Val < 0.05 & logFC < -1                        ~ "Down",
        TRUE                                                  ~ "NS"
      )
    )

  top_genes <- deg |>
    dplyr::filter(status != "NS") |>
    dplyr::slice_max(neg_log10_p, n = 15)

  vol_colors <- c(
    "Up"             = "#D6604D",
    "Down"           = "#4393C3",
    "Signature (up)" = "#B2182B",
    "Signature (down)" = "#2166AC",
    "NS"             = "grey75"
  )

  p <- ggplot(deg, aes(logFC, neg_log10_p, colour = status)) +
    geom_point(size = 0.6, alpha = 0.6) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed",
               linewidth = 0.4, colour = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               linewidth = 0.4, colour = "grey40") +
    ggrepel::geom_text_repel(
      data = top_genes, aes(label = gene),
      size = 2.5, max.overlaps = 20, colour = "black",
      segment.size = 0.3
    ) +
    scale_colour_manual(values = vol_colors, name = NULL) +
    labs(
      title    = paste0(subtype, " vs GTEx Healthy"),
      subtitle = "FDR < 0.05 threshold shown; labelled: top 15 by significance",
      x        = "log\u2082 fold-change",
      y        = "-log\u2081\u2080(adjusted p-value)"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", colour = subtype_colors[[subtype]]),
          legend.position = "bottom")
  p
}

# ---------------------------------------------------------------------------
# Fig 4 — Top DEG heatmap (ComplexHeatmap)
# ---------------------------------------------------------------------------

#' Figure 4: heatmap of top 20 up/down DEGs per subtype vs GTEx, tumors only
#'
#' Selects top 20 up and 20 down genes per subtype by t-statistic from the
#' vs_GTEx contrasts, takes the union, and draws a ComplexHeatmap with
#' PAM50 column annotation.
#'
#' @param cohort   List. $expr and $meta from cohort.rds.
#' @param deg_dir  Character. Path to results/tables/deg/.
#' @return Heatmap object.
#' @examples
#' ht <- build_deg_heatmap(cohort, here::here("results", "tables", "deg"))
build_deg_heatmap <- function(cohort, deg_dir) {
  message("  Fig 4: DEG heatmap ...")

  expr <- cohort$expr
  meta <- cohort$meta |>
    dplyr::filter(group == "Tumor", PAM50 %in% PAM50_LEVELS)

  # Collect top genes per subtype
  top_genes <- lapply(PAM50_LEVELS, function(st) {
    f   <- file.path(deg_dir, paste0(st, "_vs_GTEx.csv"))
    deg <- read.csv(f, stringsAsFactors = FALSE)
    up  <- deg |> dplyr::arrange(desc(t)) |> head(20) |> dplyr::pull(gene)
    dn  <- deg |> dplyr::arrange(t)       |> head(20) |> dplyr::pull(gene)
    c(up, dn)
  }) |> unlist() |> unique()
  top_genes <- intersect(top_genes, rownames(expr))

  mat <- expr[top_genes, meta$sample_id]
  # Z-score per gene
  mat <- t(scale(t(mat)))
  mat[mat >  3] <-  3
  mat[mat < -3] <- -3

  col_ann <- HeatmapAnnotation(
    PAM50 = meta$PAM50,
    col   = list(PAM50 = subtype_colors),
    annotation_legend_param = list(PAM50 = list(title = "PAM50"))
  )

  col_fun <- colorRamp2(c(-3, 0, 3),
                        c("#2166AC", "white", "#D6604D"))

  ht <- Heatmap(
    mat,
    name                    = "Z-score",
    col                     = col_fun,
    top_annotation          = col_ann,
    show_column_names       = FALSE,
    show_row_names          = TRUE,
    row_names_gp            = grid::gpar(fontsize = 7),
    cluster_rows            = TRUE,
    cluster_columns         = TRUE,
    column_split            = meta$PAM50,
    column_title_gp         = grid::gpar(fontsize = 10, fontface = "bold"),
    column_title            = "Top 20 up/down DEGs per subtype vs GTEx Healthy",
    heatmap_legend_param    = list(title = "Z-score")
  )
  ht
}

# ---------------------------------------------------------------------------
# Fig 5 — Signature gene heatmap
# ---------------------------------------------------------------------------

#' Figure 5: heatmap of all ElasticNet signature genes across tumor samples
#'
#' Unions the 4 subtype signatures, extracts expression for tumor samples,
#' and draws a ComplexHeatmap annotated by PAM50 and mean_coef direction.
#'
#' @param cohort   List. $expr and $meta.
#' @param sig_dir  Character. Path to results/tables/signatures/.
#' @return Heatmap object.
#' @examples
#' ht <- build_signature_heatmap(cohort, here::here("results","tables","signatures"))
build_signature_heatmap <- function(cohort, sig_dir) {
  message("  Fig 5: signature heatmap ...")

  expr <- cohort$expr
  meta <- cohort$meta |>
    dplyr::filter(group == "Tumor", PAM50 %in% PAM50_LEVELS)

  sig_list <- lapply(PAM50_LEVELS, function(st) {
    d <- read.csv(file.path(sig_dir, paste0(st, "_signature.csv")),
                  stringsAsFactors = FALSE)
    d$subtype <- st
    d
  }) |> dplyr::bind_rows()

  # For genes shared across subtypes, keep the one with highest |mean_coef|
  sig_dedup <- sig_list |>
    dplyr::group_by(gene) |>
    dplyr::slice_max(abs(mean_coef), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  genes <- intersect(sig_dedup$gene, rownames(expr))
  mat   <- expr[genes, meta$sample_id]
  mat   <- t(scale(t(mat)))
  mat[mat >  3] <-  3
  mat[mat < -3] <- -3

  coef_info <- sig_dedup |>
    dplyr::filter(gene %in% genes) |>
    dplyr::arrange(match(gene, genes))

  row_ann <- rowAnnotation(
    Direction = ifelse(coef_info$mean_coef > 0, "Up", "Down"),
    col = list(Direction = c(Up = "#D6604D", Down = "#4393C3")),
    annotation_legend_param = list(Direction = list(title = "Coef"))
  )

  col_ann <- HeatmapAnnotation(
    PAM50 = meta$PAM50,
    col   = list(PAM50 = subtype_colors)
  )

  col_fun <- colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#D6604D"))

  ht <- Heatmap(
    mat,
    name              = "Z-score",
    col               = col_fun,
    top_annotation    = col_ann,
    right_annotation  = row_ann,
    show_column_names = FALSE,
    show_row_names    = TRUE,
    row_names_gp      = grid::gpar(fontsize = 7),
    cluster_rows      = TRUE,
    cluster_columns   = TRUE,
    column_split      = meta$PAM50,
    column_title_gp   = grid::gpar(fontsize = 10, fontface = "bold"),
    column_title      = "ElasticNet signature genes across PAM50 tumor subtypes"
  )
  ht
}

# ---------------------------------------------------------------------------
# Fig 6 — fgsea NES dot plot
# ---------------------------------------------------------------------------

#' Figure 6: dot plot of top fgsea pathways per subtype (vs GTEx, Hallmark only)
#'
#' Shows top 10 pathways by |NES| per subtype at padj < 0.05. Point size
#' encodes -log10(padj), color encodes NES direction.
#'
#' @param gsea_dir Character. Path to results/tables/fgsea/.
#' @return ggplot object.
#' @examples
#' p <- plot_fgsea_dotplot(here::here("results","tables","fgsea"))
plot_fgsea_dotplot <- function(gsea_dir) {
  message("  Fig 6: fgsea dot plot ...")

  files  <- list.files(gsea_dir, pattern = "[.]csv$", full.names = TRUE)
  gsea   <- lapply(files, read.csv, stringsAsFactors = FALSE) |>
    dplyr::bind_rows()

  # subtype is encoded in the contrast column e.g. "LumA_vs_GTEx"
  # collection value is "Hallmark" (full word, not "H")
  gsea <- gsea |>
    dplyr::mutate(
      subtype = sub("_vs_.*$", "", contrast)
    )

  top <- gsea |>
    dplyr::filter(
      grepl("vs_GTEx", contrast),
      collection == "Hallmark",
      padj < 0.05
    ) |>
    dplyr::group_by(subtype) |>
    dplyr::slice_max(abs(NES), n = 10, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      pathway_short = sub("^HALLMARK_", "", pathway),
      pathway_short = gsub("_", " ", pathway_short),
      direction     = ifelse(NES > 0, "Activated", "Suppressed")
    )

  if (nrow(top) == 0) {
    message("    no Hallmark pathways at padj<0.05 vs GTEx — relaxing to padj<0.25")
    top <- gsea |>
      dplyr::filter(grepl("vs_GTEx", contrast), collection == "Hallmark",
                    padj < 0.25) |>
      dplyr::group_by(subtype) |>
      dplyr::slice_max(abs(NES), n = 10, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        pathway_short = sub("^HALLMARK_", "", pathway),
        pathway_short = gsub("_", " ", pathway_short),
        direction     = ifelse(NES > 0, "Activated", "Suppressed")
      )
  }

  p <- ggplot(top,
              aes(x = factor(subtype, levels = PAM50_LEVELS),
                  y = reorder(pathway_short, NES),
                  colour = NES,
                  size   = -log10(padj + 1e-10))) +
    geom_point() +
    scale_colour_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                           midpoint = 0, name = "NES") +
    scale_size_continuous(name = "-log\u2081\u2080(padj)", range = c(1.5, 6)) +
    labs(
      title    = "Top Hallmark pathways per subtype vs GTEx Healthy",
      subtitle = "Top 10 by |NES|, padj < 0.05",
      x = "PAM50 subtype", y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x  = element_text(colour = subtype_colors[PAM50_LEVELS],
                                  face = "bold"),
      plot.title   = element_text(face = "bold"),
      panel.grid.major.x = element_blank()
    )
  p
}

# ---------------------------------------------------------------------------
# Fig 7 — Cox HR forest plot
# ---------------------------------------------------------------------------

#' Figure 7: forest plot of Cox HR per subtype
#'
#' Reads the 4 subtype Cox CSVs and produces a forest plot showing HR with
#' 95% CI. Reference line at HR = 1. Point color by subtype.
#'
#' @param surv_dir Character. Path to results/tables/survival/.
#' @return ggplot object.
#' @examples
#' p <- plot_cox_forest(here::here("results","tables","survival"))
plot_cox_forest <- function(surv_dir) {
  message("  Fig 7: Cox forest plot ...")

  files <- list.files(surv_dir, pattern = "[.]csv$", full.names = TRUE)
  cox   <- lapply(files, read.csv, stringsAsFactors = FALSE) |>
    dplyr::bind_rows() |>
    dplyr::mutate(subtype = factor(subtype, levels = rev(PAM50_LEVELS)))

  p <- ggplot(cox, aes(x = gene_score_hr,
                       y = subtype,
                       colour = subtype)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               linewidth = 0.5, colour = "grey40") +
    geom_errorbar(aes(xmin = gene_score_lower95,
                      xmax = gene_score_upper95),
                  width = 0.2, linewidth = 0.8,
                  orientation = "y") +
    geom_point(size = 4) +
    geom_text(aes(label = paste0("p=", signif(gene_score_pval, 2))),
              hjust = -0.2, size = 3.2, colour = "black") +
    scale_colour_manual(values = subtype_colors, guide = "none") +
    scale_x_log10() +
    labs(
      title    = "Cox proportional-hazards: signature score vs Overall Survival",
      subtitle = "HR per 1 SD increase in mean signature expression (log\u2082-TPM)",
      x        = "Hazard Ratio (95% CI, log scale)",
      y        = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title  = element_text(face = "bold"),
          axis.text.y = element_text(
            colour = subtype_colors[as.character(rev(PAM50_LEVELS))],
            face   = "bold",
            size   = 12
          ))
  p
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

assert_all_stages_complete()

message("── Stage 08 | loading data ──────────────────────────────────────────────")
cohort     <- readRDS(here::here("data", "processed", "cohort.rds"))
umap_obj   <- readRDS(here::here("data", "processed", "umap_clusters.rds"))
umap_df    <- umap_obj$umap   # already contains group, PAM50, subcluster
message("  cohort loaded: ", nrow(cohort$meta), " samples")
message("  UMAP loaded:   ", nrow(umap_df), " rows")
message("── Stage 08 | data loaded ───────────────────────────────────────────────")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Check for cowplot (needed for Fig 1 panel assembly)
if (!requireNamespace("cowplot", quietly = TRUE))
  stop("Package 'cowplot' required. Install with renv::install('cowplot').")
if (!requireNamespace("ggrepel", quietly = TRUE))
  stop("Package 'ggrepel' required. Install with renv::install('ggrepel').")

message("── Stage 08 | Fig 1: UMAP overview ─────────────────────────────────────")
p1 <- plot_umap_overview(umap_df)
save_fig(p1, "fig1_umap_overview", width = 14, height = 6)

message("── Stage 08 | Fig 2: DEG counts ─────────────────────────────────────────")
p2 <- plot_deg_counts(here::here("results", "tables", "deg"))
save_fig(p2, "fig2_deg_counts", width = 10, height = 12)

message("── Stage 08 | Fig 3: Volcano plots ──────────────────────────────────────")
deg_dir <- here::here("results", "tables", "deg")
sig_dir <- here::here("results", "tables", "signatures")
for (st in PAM50_LEVELS) {
  message("  volcano: ", st, " vs GTEx")
  p3 <- plot_volcano(
    deg_path = file.path(deg_dir, paste0(st, "_vs_GTEx.csv")),
    sig_path = file.path(sig_dir, paste0(st, "_signature.csv")),
    subtype  = st
  )
  save_fig(p3, paste0("fig3_volcano_", st, "_vs_GTEx"), width = 9, height = 7)
}

message("── Stage 08 | Fig 4: DEG heatmap ────────────────────────────────────────")
ht4 <- build_deg_heatmap(cohort, here::here("results", "tables", "deg"))
save_heatmap(ht4, "fig4_deg_heatmap",
             width_px = 3600, height_px = 3000,
             width_in = 14,   height_in = 12)

message("── Stage 08 | Fig 5: Signature heatmap ──────────────────────────────────")
ht5 <- build_signature_heatmap(cohort,
                                here::here("results", "tables", "signatures"))
save_heatmap(ht5, "fig5_signature_heatmap",
             width_px = 3000, height_px = 2800,
             width_in = 12,   height_in = 11)

message("── Stage 08 | Fig 6: fgsea dot plot ─────────────────────────────────────")
p6 <- plot_fgsea_dotplot(here::here("results", "tables", "fgsea"))
save_fig(p6, "fig6_fgsea_dotplot", width = 12, height = 10)

message("── Stage 08 | Fig 7: Cox forest plot ────────────────────────────────────")
p7 <- plot_cox_forest(here::here("results", "tables", "survival"))
save_fig(p7, "fig7_cox_forest", width = 9, height = 5)

message("── Stage 08 | output summary ────────────────────────────────────────────")
figs <- list.files(FIG_DIR, pattern = "[.](png|pdf)$")
message("  total figure files: ", length(figs))
message("✅ Stage 08 complete. Pipeline finished.")
