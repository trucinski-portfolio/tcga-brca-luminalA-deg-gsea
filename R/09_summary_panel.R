# R/09_summary_panel.R
# Summary panel: assembles existing pipeline figure .pngs into a single
# at-a-glance composite figure using cowplot::draw_image (no analysis re-run).
#
# Inputs:  results/figures/*.png  (must exist — run 08_visualization.R first)
# Outputs: results/figures/summary_panel.png  (300 dpi)
#          results/figures/summary_panel.pdf  (vector)

library(here)
library(cowplot)
library(ggplot2)

FIG_DIR <- here::here("results", "figures")

# ---------------------------------------------------------------------------
# Stage gate
# ---------------------------------------------------------------------------

#' Stop with a clear error if required input figures are missing
#'
#' @return Invisible NULL.
#' @examples
#' assert_figures_exist()
assert_figures_exist <- function() {
  message("── Summary panel | stage gate ───────────────────────────────────────────")
  required <- c(
    "fig1_umap_overview.png",
    "fig4_deg_heatmap.png",
    "fig6_fgsea_dotplot.png",
    "fig5_signature_heatmap.png",
    "fig7_cox_forest.png",
    "wgcna_module_trait_heatmap.png"
  )
  missing <- required[!file.exists(file.path(FIG_DIR, required))]
  if (length(missing) > 0)
    stop("Missing figures: ", paste(missing, collapse = ", "),
         "\nRun R/08_visualization.R (and R/06_wgcna.R) first.")

  if (!requireNamespace("magick", quietly = TRUE))
    stop("Package 'magick' required for cowplot::draw_image().\n",
         "Install with: renv::install('magick')")

  message("  all required figures found.")
  message("── Summary panel | gate passed ──────────────────────────────────────────")
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Build panel
# ---------------------------------------------------------------------------

#' Load a .png as a cowplot ggdraw image object
#'
#' @param filename Character. Filename within FIG_DIR (without path).
#' @return ggplot object suitable for cowplot::plot_grid.
#' @examples
#' img <- load_img("fig1_umap_overview.png")
load_img <- function(filename) {
  path <- file.path(FIG_DIR, filename)
  cowplot::ggdraw() + cowplot::draw_image(path)
}

#' Build the 2-column × 3-row composite summary panel
#'
#' Layout:
#'   Row 1:  UMAP overview            |  DEG heatmap
#'   Row 2:  fgsea dotplot            |  Signature heatmap
#'   Row 3:  Cox forest plot          |  WGCNA module-trait heatmap
#'
#' @return ggplot object (cowplot grid).
#' @examples
#' panel <- build_summary_panel()
build_summary_panel <- function() {
  message("── Summary panel | loading figures ──────────────────────────────────────")

  img_umap   <- load_img("fig1_umap_overview.png")
  img_deg    <- load_img("fig4_deg_heatmap.png")
  img_fgsea  <- load_img("fig6_fgsea_dotplot.png")
  img_sig    <- load_img("fig5_signature_heatmap.png")
  img_cox    <- load_img("fig7_cox_forest.png")
  img_wgcna  <- load_img("wgcna_module_trait_heatmap.png")

  message("  all 6 figures loaded.")

  # Row labels as small ggplot text panels
  label <- function(txt) {
    cowplot::ggdraw() +
      cowplot::draw_label(txt, fontface = "bold", size = 11,
                          x = 0.02, hjust = 0, colour = "grey30")
  }

  # Assemble: each row is [label | figure | label | figure]
  # Use rel_heights to give label rows minimal height
  row1 <- cowplot::plot_grid(
    label("A  UMAP — sample groups and PAM50 subtypes"),
    label("B  Top DEGs per subtype vs GTEx (z-scored)"),
    nrow = 1, rel_widths = c(1, 1)
  )
  row1_figs <- cowplot::plot_grid(img_umap, img_deg, nrow = 1, rel_widths = c(1.3, 1))

  row2 <- cowplot::plot_grid(
    label("C  Hallmark pathway enrichment (fgsea NES)"),
    label("D  ElasticNet signature genes across tumors"),
    nrow = 1, rel_widths = c(1, 1)
  )
  row2_figs <- cowplot::plot_grid(img_fgsea, img_sig, nrow = 1)

  row3 <- cowplot::plot_grid(
    label("E  Cox PH: signature score vs Overall Survival"),
    label("F  WGCNA module–trait correlations"),
    nrow = 1, rel_widths = c(1, 1)
  )
  row3_figs <- cowplot::plot_grid(img_cox, img_wgcna, nrow = 1, rel_widths = c(0.8, 1))

  # Title strip
  title_strip <- cowplot::ggdraw() +
    cowplot::draw_label(
      "TCGA-BRCA PAM50 Subtype Analysis \u2014 Pipeline Summary",
      fontface = "bold", size = 14, colour = "grey15"
    )

  # Stack everything
  panel <- cowplot::plot_grid(
    title_strip,
    row1,    row1_figs,
    row2,    row2_figs,
    row3,    row3_figs,
    ncol        = 1,
    rel_heights = c(0.04, 0.03, 1, 0.03, 1, 0.03, 0.7)
  )

  message("── Summary panel | panel assembled ──────────────────────────────────────")
  panel
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

assert_figures_exist()

panel <- build_summary_panel()

message("── Summary panel | saving ───────────────────────────────────────────────")

print(panel)

ggsave(
  filename = file.path(FIG_DIR, "summary_panel.png"),
  plot     = panel,
  width    = 22, height = 26, dpi = 300, bg = "white"
)

ggsave(
  filename = file.path(FIG_DIR, "summary_panel.pdf"),
  plot     = panel,
  width    = 22, height = 26, bg = "white"
)

message("  saved: summary_panel.png + summary_panel.pdf")
message("✅ Summary panel complete.")
