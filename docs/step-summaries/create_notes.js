const { Document, Packer, Paragraph, TextRun } = require("docx");
const fs = require("fs");
const path = require("path");

const OUT = __dirname;

function p(text, opts = {}) {
  return new Paragraph({
    spacing: { after: 160 },
    children: [new TextRun({ text, font: "Arial", size: 22, ...opts })],
  });
}

function blank() {
  return new Paragraph({ spacing: { after: 160 }, children: [] });
}

function bold(text) {
  return p(text, { bold: true });
}

function code(text) {
  return new Paragraph({
    spacing: { after: 100 },
    children: [new TextRun({ text, font: "Courier New", size: 18 })],
  });
}

async function write(filename, children) {
  const doc = new Document({
    sections: [{
      properties: {
        page: {
          size: { width: 12240, height: 15840 },
          margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 },
        },
      },
      children,
    }],
  });
  const buf = await Packer.toBuffer(doc);
  fs.writeFileSync(path.join(OUT, filename), buf);
  console.log("wrote", filename);
}

// ─── 01 preprocessing ────────────────────────────────────────────────────────
write("01_preprocessing_notes.docx", [
  bold("01_preprocessing.R"),
  blank(),
  bold("Purpose"),
  p("Downloads the TOIL combined expression matrix, maps gene IDs from Ensembl to HGNC symbols, fetches PAM50 subtype calls, and builds the master cohort object used by every downstream stage."),
  blank(),
  bold("Produces"),
  p("data/processed/cohort.rds — a list with three elements:"),
  p("  $expr: numeric matrix, 28,344 HGNC symbols x 1,109 sample IDs, values in log2(TPM+0.001). Do not re-log."),
  p("  $meta: data.frame with columns sample_id, group (Tumor/NAT/GTEx), PAM50, sample_type"),
  p("  $contrasts: character vector of 18 contrast names"),
  blank(),
  bold("Key constants"),
  code("TOIL_HUB      <- \"toilHub\""),
  code("EXPR_DATASET  <- \"TcgaTargetGtex_rsem_gene_tpm\""),
  code("PAM50_LEVELS  <- c(\"LumA\", \"LumB\", \"Her2\", \"Basal\")"),
  code("GROUP_LEVELS  <- c(\"Tumor\", \"NAT\", \"GTEx\")"),
  blank(),
  bold("assert_stage00_complete()"),
  p("Checks that data/processed/metadata_breast.rds exists. Hard stops if missing."),
  blank(),
  bold("load_metadata()"),
  p("Reads metadata_breast.rds, extracts $pheno. The group column (Tumor/NAT/GTEx) was already assigned by Stage 00's 3-way filter."),
  blank(),
  bold("download_expression()"),
  p("Uses the authoritative UCSCXenaTools pipeline: XenaGenerate -> XenaFilter -> XenaQuery -> XenaDownload -> XenaPrepare. Sets VROOM_CONNECTION_SIZE = 268435456 (256 MB) before XenaPrepare because the matrix has ~19,000 columns per row and exceeds vroom's default 128 KB line buffer. The .gz file is cached in data/raw/ so re-runs skip the download. Immediately subsets columns to the 1,384 breast samples."),
  blank(),
  bold("map_ensembl_to_hgnc()"),
  p("Maps Ensembl IDs to HGNC symbols using AnnotationDbi::mapIds(org.Hs.eg.db). Fully local — no network required. TOIL row IDs include version suffixes like ENSG00000000003.14; these are stripped with sub(\"\\\\.\\\\d+$\", \"\", ids) before lookup. Unmapped genes are dropped. The 57 cases where multiple Ensembl IDs collapse to the same HGNC symbol are resolved by row-wise mean using rowsum / table."),
  blank(),
  bold("fetch_pam50()"),
  p("Critical: the TOIL phenotype does NOT contain PAM50 calls. detailed_category just says \"Breast Invasive Carcinoma\" for every tumor. PAM50 must be fetched separately from tcgaHub (not toilHub), dataset TCGA.BRCA.sampleMap/BRCA_clinicalMatrix, column PAM50Call_RNAseq. Values: LumA, LumB, Her2, Basal, Normal. Cached to data/processed/pam50_calls.rds. Join is on the 15-character TCGA barcode."),
  blank(),
  bold("build_meta()"),
  p("Joins pheno with PAM50 calls. Normal-like tumors (PAM50 == \"Normal\") are excluded from the Tumor arm only. NAT and GTEx are never filtered by PAM50 — they get PAM50 = NA. Final counts: Tumor 817 (LumA 420, LumB 192, Her2 66, Basal 139), NAT 113, GTEx 179."),
  blank(),
  bold("build_contrasts()"),
  p("Returns a plain character vector of 18 contrast names. The 4 Within-subtype names are placeholders that Stage 02 k-means fills in."),
  blank(),
  bold("save_cohort()"),
  p("Aligns expr columns to meta rows via intersect(), verifies identity with stopifnot(identical(...)) as a regression guard, then saves."),
]);

// ─── 02 clustering ───────────────────────────────────────────────────────────
write("02_clustering_notes.docx", [
  bold("02_clustering.R"),
  blank(),
  bold("Purpose"),
  p("UMAP dimensionality reduction and within-subtype sub-subcluster detection. Produces all UMAP visualizations and the cluster assignments that Stage 03 uses for its Within-subtype contrasts."),
  blank(),
  bold("Produces"),
  p("data/processed/umap_clusters.rds"),
  p("  $umap: data.frame with columns sample_id, UMAP1, UMAP2, subtype, subcluster"),
  p("  $n_per_subtype: named integer vector"),
  p("results/figures/ — 6 figure pairs (.png 300dpi + .pdf):"),
  p("  umap_by_group, umap_by_subtype, umap_subclusters_LumA/LumB/Her2/Basal"),
  blank(),
  bold("Key constants"),
  code("N_HVG        <- 5000L   # highly variable genes fed to PCA"),
  code("N_PCA        <- 50L     # PCA components fed to UMAP"),
  code("K_SUBCLUSTER <- 2L      # k-means k within each PAM50 subtype"),
  blank(),
  bold("select_hvg()"),
  p("Picks top 5,000 genes by variance across all 1,109 samples using matrixStats::rowVars. Uses all sample groups (not just Tumor) to maximise between-group separation in UMAP space."),
  blank(),
  bold("run_pca()"),
  p("prcomp(t(hvg), center=TRUE, scale.=FALSE, rank.=50). Reports percentage variance explained by PC1-5. Returns the sample x 50 PC score matrix."),
  blank(),
  bold("run_umap()"),
  p("uwot::umap(pcs, n_neighbors=30, min_dist=0.3, n_components=2, seed=42). n_neighbors=30 is standard for ~1,000-sample datasets. min_dist=0.3 allows moderate cluster separation without over-compressing local structure."),
  blank(),
  bold("detect_subclusters()"),
  p("For each PAM50 subtype, runs k-means (k=2) on the 2D UMAP coordinates of just that subtype's samples. Labels are formatted as LumA_C1, LumA_C2, etc. Controls (NAT, GTEx) get NA. These labels feed Stage 03's Within-subtype limma contrasts."),
  blank(),
  bold("save_plot(p, stem)"),
  p("The display-and-save pattern: print(p) fires to screen, then ggsave saves .png (300 dpi) and .pdf. Every plot goes through here."),
  blank(),
  bold("Plot layering strategy"),
  p("Controls are drawn in grey first (geom_point), then tumor subtypes on top in subtype_colors. For subcluster plots, all other subtypes are greyed out so the focal subtype's two clusters are the only colored points."),
]);

// ─── 03 deg ──────────────────────────────────────────────────────────────────
write("03_deg_notes.docx", [
  bold("03_deg.R"),
  blank(),
  bold("Purpose"),
  p("Differential expression analysis for all 18 contrasts using limma-trend. The workhorse output that feeds Stages 04, 05, 07, and 08."),
  blank(),
  bold("Produces"),
  p("18 CSV files in results/tables/deg/. Columns: gene, logFC, AveExpr, t, P.Value, adj.P.Val, B. All 28,344 genes written per file — no pre-filtering. Downstream stages apply their own thresholds."),
  blank(),
  bold("Model design"),
  p("Single 6-level factor model across all sample groups:"),
  code("group_label <- factor(label, levels = c(\"LumA\",\"LumB\",\"Her2\",\"Basal\",\"NAT\",\"GTEx\"))"),
  code("design <- model.matrix(~ 0 + group_label)"),
  code("fit    <- lmFit(expr, design)"),
  code("fit    <- eBayes(fit, trend = TRUE)"),
  blank(),
  bold("Why trend=TRUE"),
  p("eBayes(trend=TRUE) is appropriate for log2-TPM data. It fits a mean-variance trend line to the prior variance, making the moderated t-statistics valid for this data type. Never use trend=FALSE or DESeq2 on this data."),
  blank(),
  bold("extract_contrast()"),
  p("For each of the 14 main contrasts (subtype vs NAT, subtype vs GTEx, subtype vs subtype, NAT vs GTEx), calls makeContrasts(contrasts=\"LumA - NAT\", levels=colnames(design)), then contrasts.fit() and eBayes(trend=TRUE) again. topTable(number=Inf, sort.by=\"none\") returns all genes in original order."),
  blank(),
  bold("run_within_contrast()"),
  p("For the 4 Within-subtype contrasts, subsets to just that PAM50 subtype's samples, builds a 2-level factor from Stage 02 k-means subcluster labels (LumA_C1 vs LumA_C2), and fits a fresh lmFit on that subset. Returns NULL with a warning if fewer than 2 subclusters exist."),
  blank(),
  bold("Verification"),
  p("Checks that all 18 expected filenames exist by name (not just that the count is 18), so stale files from previous runs don't cause false passes."),
]);

// ─── 04 feature selection ────────────────────────────────────────────────────
write("04_feature_selection_notes.docx", [
  bold("04_feature_selection.R"),
  blank(),
  bold("Purpose"),
  p("Bootstrap ElasticNet stability selection to find minimal gene signatures for each PAM50 subtype. Answers: which genes robustly distinguish this subtype from all other tumors across many simulated cohorts?"),
  blank(),
  bold("Produces"),
  p("4 CSV files in results/tables/signatures/ — LumA_signature.csv, LumB_signature.csv, Her2_signature.csv, Basal_signature.csv"),
  p("Columns: gene, selection_frequency, mean_coef"),
  p("Final counts: LumA 30 genes, LumB 27, Her2 36, Basal 34"),
  blank(),
  bold("Key constants"),
  code("ALPHA        <- 0.5     # ElasticNet (0=Ridge, 1=Lasso)"),
  code("N_FOLDS      <- 10L     # cross-validation folds"),
  code("N_BOOT       <- 100L    # bootstrap iterations"),
  code("FREQ_CUTOFF  <- 0.80    # selection frequency threshold"),
  blank(),
  bold("Why 0.80?"),
  p("Based on Meinshausen & Buhlmann (2010), who derive 0.80 as the principled threshold for controlling false positives in stability selection. Applied uniformly across all subtypes — methodologically cleaner for a thesis defense than mixed thresholds."),
  blank(),
  bold("Step 1 — initial cv.glmnet"),
  p("cv.glmnet runs once on all 28,344 genes with 10-fold cross-validation to find lambda.1se. This is the only expensive call."),
  blank(),
  bold("Step 2 — reduce feature space"),
  p("Reduces from 28,344 genes to ~100-300 by keeping only genes active at lambda.min (guaranteed superset of lambda.1se genes, so no signal is lost). Uses direct beta matrix access:"),
  code("lmin_idx     <- which.min(abs(cv_fit$lambda - cv_fit$lambda.min))"),
  code("beta_min_col <- cv_fit$glmnet.fit$beta[, lmin_idx]"),
  code("active_genes <- rownames(cv_fit$glmnet.fit$beta)[beta_min_col != 0]"),
  p("Important: as.numeric(coef(cv_fit, s=\"lambda.min\")) was found to return all zeros for sparse dgCMatrix objects — always use the direct beta matrix indexing above, never coef()."),
  blank(),
  bold("Step 3 — 100 bootstrap iterations"),
  p("Each iteration: sample with replacement, fit glmnet at fixed lambda.1se on the reduced feature space, record which genes had non-zero coefficients. selection_frequency = fraction of 100 runs where the gene was selected. mean_coef captures average effect direction and magnitude."),
  blank(),
  bold("Resume behavior"),
  p("Script skips already-completed subtypes by checking if the output CSV exists. Partial runs can resume without re-running completed subtypes."),
  blank(),
  bold("Why Her2 and Basal have broader signatures"),
  p("Her2 has only 66 samples and Basal 139. Smaller positive class means the ElasticNet needs more genes to maintain discriminative power at the same lambda.1se, which produces slightly wider signatures. This is expected and not a bug."),
]);

// ─── 05 fgsea ────────────────────────────────────────────────────────────────
write("05_fgsea_notes.docx", [
  bold("05_fgsea.R"),
  blank(),
  bold("Purpose"),
  p("Gene set enrichment analysis for each PAM50 subtype. Answers: what biological processes are transcriptionally activated or suppressed in this subtype relative to healthy tissue and adjacent normal?"),
  blank(),
  bold("Produces"),
  p("4 CSV files in results/tables/fgsea/ — LumA_fgsea.csv, etc."),
  p("Columns: contrast, collection, pathway, NES, pval, padj, size, leadingEdge"),
  p("Each file contains results for BOTH contrasts (vs_GTEx and vs_NAT) for that subtype, distinguished by the contrast column."),
  p("Results: 58-93 significant pathways per subtype per contrast at padj < 0.05."),
  blank(),
  bold("Gene sets"),
  p("MSigDB Hallmark (50 curated oncology pathways) + KEGG Legacy (186 metabolic/signaling pathways) = 236 total. Loaded via msigdbr which caches locally after the first run — no network needed after that."),
  p("API note: as of msigdbr v10.0.0, use collection= and subcollection= parameters. The older category= and subcategory= arguments are deprecated and will error."),
  blank(),
  bold("Ranking metric"),
  p("limma moderated t-statistic — NOT log2FC. The t-statistic captures both effect size and precision. Low-variance genes with real effects rank higher than noisy genes with big fold-changes. build_rank_vector() drops NA t-stats, resolves duplicates by highest |t|, sorts descending."),
  blank(),
  bold("fgsea call"),
  code("fgsea(pathways=gene_sets, stats=ranks,"),
  code("      minSize=15, maxSize=500,"),
  code("      nPermSimple=1000, eps=0)"),
  p("eps=0 ensures accurate p-value estimation for pathways with p < 1e-50, which occurs given the large sample sizes and strong cancer-vs-normal signal. Without eps=0, those pathways get reported as p=0."),
  blank(),
  bold("Why both contrasts?"),
  p("vs_GTEx = absolute dysregulation relative to truly healthy breast tissue. Most clinically interpretable. vs_NAT = dysregulation relative to cancer-adjacent tissue, capturing what is changed beyond the field effect. Together they characterize the full biology."),
  blank(),
  bold("leadingEdge handling"),
  p("The leadingEdge list column (genes driving the enrichment score) is collapsed to semicolon-separated strings for CSV compatibility: paste(x, collapse=\";\")."),
  blank(),
  bold("Tie rate note"),
  p("About 11% of the rank vector values are tied. This is normal for limma t-statistics on log2-TPM data and has negligible effect across 1,000 permutations. fgsea prints a warning about it — this is expected, not an error."),
]);

// ─── 06 wgcna ────────────────────────────────────────────────────────────────
write("06_wgcna_notes.docx", [
  bold("06_wgcna.R"),
  blank(),
  bold("Purpose"),
  p("Weighted Gene Co-expression Network Analysis. Finds clusters of genes that move together across tumor samples (modules), then tests which modules are correlated with PAM50 subtype identity."),
  blank(),
  bold("Produces"),
  p("data/processed/wgcna_modules.rds — module gene assignments, eigengenes, module-trait correlation matrix, power used"),
  p("results/figures/wgcna_soft_threshold.png + .pdf — R2-vs-power diagnostic curve"),
  p("results/figures/wgcna_module_trait_heatmap.png + .pdf — modules x PAM50 traits heatmap"),
  blank(),
  bold("Confirmed output"),
  p("5 modules (excluding grey), power = 8 (R2 = 0.907), 20/24 significant module-trait associations at p < 0.05. Module sizes: grey=2061 (unassigned), M1=775, M2=609, M3=598, M4=562, M5=395."),
  blank(),
  bold("Key constants"),
  code("N_HVG        <- 5000L    # top MAD genes"),
  code("NETWORK_TYPE <- \"signed\" # preserves co-expression direction"),
  code("MIN_MODULE   <- 30L      # minimum genes per module"),
  code("MERGE_CUT    <- 0.25     # dendrogram merge height"),
  code("R2_TARGET    <- 0.85     # scale-free topology target"),
  blank(),
  bold("Why signed network?"),
  p("Keeps positively and negatively co-expressed genes in separate modules. Biologically more interpretable than unsigned networks that lump both directions together."),
  blank(),
  bold("Why MAD instead of variance for gene selection?"),
  p("MAD (median absolute deviation) is more robust to outlier samples than variance. Computed with matrixStats::rowMads()."),
  blank(),
  bold("prepare_expr()"),
  p("Subsets to 817 Tumor samples, selects top 5,000 genes by MAD, then runs WGCNA::goodSamplesGenes() to remove any zero-variance genes or outlier samples before network construction."),
  blank(),
  bold("pick_soft_threshold()"),
  p("Tests powers 1-30 to find the lowest power where the network approximates a scale-free topology (R2 >= 0.85). This is a WGCNA biological validity requirement. Produces and saves the R2-vs-power diagnostic immediately so you can see it while the network builds."),
  blank(),
  bold("build_network()"),
  p("blockwiseModules() with maxBlockSize=6000, which fits all 5,000 genes in one block. The grey module (label 0) is always the garbage collector for genes that do not fit any module."),
  blank(),
  bold("correlate_modules_traits()"),
  p("Creates binary trait indicators (LumA=1/0, LumB=1/0, etc.) and Pearson-correlates each module's eigengene with each trait. A module eigengene is the first principal component of that module — one summary expression value per sample per module."),
  blank(),
  bold("plot_module_trait_heatmap()"),
  p("ComplexHeatmap with blue-white-red color scale (-1 to +1), significance asterisks in cells (* p<0.05, ** p<0.01, *** p<0.001), PAM50 subtype colors as column annotation. draw(ht) prints to screen; then png() and pdf() device blocks save the same plot."),
]);

// ─── 07 survival ─────────────────────────────────────────────────────────────
write("07_survival_notes.docx", [
  bold("07_survival.R"),
  blank(),
  bold("Purpose"),
  p("Cox proportional-hazards survival analysis. For each PAM50 subtype, scores each tumor sample on its ElasticNet signature (mean log2-TPM of the signature genes), then tests whether high vs. low scorers have significantly different Overall Survival."),
  blank(),
  bold("Produces"),
  p("results/tables/survival/LumA_cox.csv, LumB_cox.csv, Her2_cox.csv, Basal_cox.csv"),
  p("Columns: subtype, gene_score_hr, gene_score_lower95, gene_score_upper95, gene_score_pval, n_events, n_total"),
  p("results/figures/survival_km_{Subtype}.png + .pdf — Kaplan-Meier curves split at median signature score"),
  blank(),
  bold("fetch_survival()"),
  p("Downloads the TCGA BRCA clinical matrix from tcgaHub (same dataset used for PAM50 in Stage 01). Extracts columns sampleID, OS (event: 0=alive, 1=dead), OS.time (days). Filters to os_time > 0. Caches to data/processed/survival.rds."),
  blank(),
  bold("compute_signature_score()"),
  p("For a given set of samples, takes the mean log2-TPM across all signature genes present in the expression matrix. Any signature genes absent from the matrix are reported and skipped. Returns a named numeric vector (sample -> score)."),
  blank(),
  bold("run_cox()"),
  p("Joins signature scores with survival data, then scales the score to SD units using scale(). Fits:"),
  code("coxph(Surv(os_time, os_event) ~ score_scaled, data = df)"),
  p("Scaling to SD units makes the HR interpretable as \"hazard ratio per 1 standard deviation increase in signature expression.\" Extracts HR, 95% CI, and p-value from summary(fit)$conf.int and $coefficients. Warns if fewer than 5 events (estimates unreliable at that point)."),
  blank(),
  bold("plot_km()"),
  p("Dichotomises samples into High/Low groups at the within-subtype median score. Uses survminer::ggsurvplot() with a risk table, confidence intervals, and log-rank p-value displayed on the plot. The High group color uses that subtype's entry in subtype_colors for visual consistency across the whole pipeline."),
  p("Important: survminer objects require print() inside the png()/pdf() device block — ggsave does not work for ggsurvplot output."),
  blank(),
  bold("Stage gate"),
  p("Checks both cohort.rds and all 4 signature CSVs before running anything. Hard stops with a clear error message if either is missing."),
]);

console.log("All 7 documents written.");
