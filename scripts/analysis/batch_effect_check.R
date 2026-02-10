# Batch Effect Assessment Script
# Checks for tissue_source_site (TSS) batch effects via PCA

# Setup
source("../setup/r_bootstrap.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# Load data
meta_fp <- file.path(REPO_ROOT, "data/processed/preprocessing_outputs/metadata_LumA_IDC_Tumor_vs_AllNormals.tsv")
expr_fp <- file.path(REPO_ROOT, "data/processed/preprocessing_outputs/expr_LumA_IDC_Tumor_vs_AllNormals.tsv")
clin_fp <- file.path(REPO_ROOT, "data/raw/preprocessing_inputs/BRCA_clinicalMatrix.tsv")

meta <- read.delim(meta_fp, stringsAsFactors = FALSE, check.names = FALSE)
expr <- read.delim(expr_fp, row.names = 1, check.names = FALSE)
clin <- read.delim(clin_fp, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)

cat("Expression matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")

# Extract tissue_source_site from clinical data
# First harmonize clinical index to match our samples
harmonize_id <- function(x) {
  gsub("\\.", "-", substr(x, 1, 16))
}

clin_ids <- harmonize_id(rownames(clin))
rownames(clin) <- clin_ids

# Get TSS for our samples
common <- intersect(meta$Sample, rownames(clin))
cat("Samples with TSS info:", length(common), "of", nrow(meta), "\n")

meta$TSS <- NA
meta$TSS[meta$Sample %in% common] <- clin[meta$Sample[meta$Sample %in% common], "tissue_source_site"]

# Also extract TSS from barcode (characters 6-7, e.g., "A2" from TCGA-A2-XXXX)
meta$TSS_from_barcode <- substr(meta$Sample, 6, 7)

cat("\nTSS distribution (from barcode):\n")
print(table(meta$TSS_from_barcode, useNA = "ifany"))

cat("\nTSS by Group:\n")
print(table(meta$TSS_from_barcode, meta$Group))

# Run PCA on expression data
expr_mat <- as.matrix(expr)

# Filter to high-variance genes for cleaner PCA
gene_vars <- apply(expr_mat, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:5000])
expr_pca_input <- t(expr_mat[top_genes, ])  # samples x genes

pca <- prcomp(expr_pca_input, scale. = TRUE, center = TRUE)

# Create PCA data frame
pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3],
  PC4 = pca$x[, 4]
)
pca_df <- merge(pca_df, meta, by = "Sample")

# Variance explained
var_explained <- summary(pca)$importance[2, 1:4] * 100

cat("\nVariance explained by top PCs:\n")
cat("  PC1:", round(var_explained[1], 1), "%\n")
cat("  PC2:", round(var_explained[2], 1), "%\n")
cat("  PC3:", round(var_explained[3], 1), "%\n")
cat("  PC4:", round(var_explained[4], 1), "%\n")

# Output directory
fig_dir <- file.path(REPO_ROOT, "results/figures/qc")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Plot 1: PCA colored by Group (expected strong separation)
p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(
    title = "PCA: Tumor vs Normal (Expected Separation)",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Normal" = "steelblue", "Tumor" = "firebrick"))

ggsave(file.path(fig_dir, "PCA_by_Group.png"), p1, width = 8, height = 6, dpi = 150)

# Plot 2: PCA colored by TSS (batch check)
# Only show top 10 most common TSS to avoid legend clutter
tss_counts <- sort(table(pca_df$TSS_from_barcode), decreasing = TRUE)
top_tss <- names(tss_counts)[1:min(10, length(tss_counts))]
pca_df$TSS_plot <- ifelse(pca_df$TSS_from_barcode %in% top_tss,
                           pca_df$TSS_from_barcode, "Other")

p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = TSS_plot)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(
    title = "PCA: Colored by Tissue Source Site (Batch Check)",
    subtitle = "Look for TSS-specific clustering within tumor/normal groups",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "PCA_by_TSS.png"), p2, width = 10, height = 6, dpi = 150)

# Plot 3: PCA PC3 vs PC4 by TSS (sometimes batch effects appear in later PCs)
p3 <- ggplot(pca_df, aes(x = PC3, y = PC4, color = TSS_plot)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(
    title = "PCA: PC3 vs PC4 by TSS",
    subtitle = "Batch effects sometimes visible in later components",
    x = paste0("PC3 (", round(var_explained[3], 1), "%)"),
    y = paste0("PC4 (", round(var_explained[4], 1), "%)")
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "PCA_PC3_PC4_by_TSS.png"), p3, width = 10, height = 6, dpi = 150)

# Plot 4: Faceted by Group to see within-group TSS patterns
p4 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = TSS_plot)) +
  geom_point(alpha = 0.7, size = 2) +
  facet_wrap(~Group) +
  labs(
    title = "PCA by TSS: Separated by Tumor/Normal",
    subtitle = "TSS clustering WITHIN groups would indicate batch effect",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "PCA_by_TSS_faceted.png"), p4, width = 12, height = 5, dpi = 150)

# Statistical test: PERMANOVA-like assessment
# Test if TSS explains variance after accounting for Group
cat("\n=== Batch Effect Assessment ===\n")

# Within tumors only: does TSS explain PC1/PC2 variance?
tumor_pca <- pca_df[pca_df$Group == "Tumor", ]

# Simple ANOVA on PC1 within tumors by TSS
if (length(unique(tumor_pca$TSS_from_barcode)) > 1) {
  aov_pc1 <- aov(PC1 ~ TSS_from_barcode, data = tumor_pca)
  cat("\nANOVA: PC1 ~ TSS (Tumors only):\n")
  print(summary(aov_pc1))

  aov_pc2 <- aov(PC2 ~ TSS_from_barcode, data = tumor_pca)
  cat("\nANOVA: PC2 ~ TSS (Tumors only):\n")
  print(summary(aov_pc2))

  # R-squared for TSS effect
  ss_tss <- sum(summary(aov_pc1)[[1]]$"Sum Sq"[1])
  ss_total <- sum(summary(aov_pc1)[[1]]$"Sum Sq")
  r2_tss <- ss_tss / ss_total
  cat("\nR² for TSS effect on PC1 (tumors):", round(r2_tss, 4), "\n")
}

cat("\nPlots saved to:", fig_dir, "\n")
cat("\n=== Summary ===\n")
cat("If TSS explains >5% of variance in top PCs within groups, consider batch correction.\n")
cat("If TSS clusters are mixed/overlapping, batch effect is likely minimal.\n")
