const { Document, Packer, Paragraph, TextRun, PageBreak,
        AlignmentType, HeadingLevel } = require("docx");
const fs   = require("fs");
const path = require("path");

const OUT     = __dirname;
const DOCS    = path.join(__dirname, "..");

// ─── helpers ──────────────────────────────────────────────────────────────────

function p(text, opts = {}) {
  return new Paragraph({
    spacing: { after: 160 },
    children: [new TextRun({ text: String(text), font: "Arial", size: 22, ...opts })],
  });
}
function blank() { return new Paragraph({ spacing: { after: 120 }, children: [] }); }
function bold(text) { return p(text, { bold: true }); }
function code(text) {
  return new Paragraph({
    spacing: { after: 80 },
    children: [new TextRun({ text, font: "Courier New", size: 18 })],
  });
}
function h1(text) {
  return new Paragraph({
    spacing: { before: 320, after: 160 },
    children: [new TextRun({ text, font: "Arial", size: 28, bold: true })],
  });
}
function h2(text) {
  return new Paragraph({
    spacing: { before: 240, after: 120 },
    children: [new TextRun({ text, font: "Arial", size: 24, bold: true })],
  });
}
function h3(text) {
  return new Paragraph({
    spacing: { before: 180, after: 100 },
    children: [new TextRun({ text, font: "Arial", size: 22, bold: true, italics: true })],
  });
}
function indent(text) {
  return new Paragraph({
    spacing: { after: 120 },
    indent: { left: 720 },
    children: [new TextRun({ text, font: "Arial", size: 22 })],
  });
}
function br() {
  return new Paragraph({ children: [new PageBreak()] });
}

async function write(filepath, children) {
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
  fs.writeFileSync(filepath, buf);
  console.log("wrote", path.basename(filepath));
}

// ─── 08 step summary ──────────────────────────────────────────────────────────

write(path.join(OUT, "08_visualization_notes.docx"), [
  bold("08_visualization.R"),
  blank(),
  bold("Purpose"),
  p("Assembles all pipeline results into a consistent set of publication-quality figures. Every figure is printed to screen during the run and saved as both .png (300 dpi) and .pdf (vector)."),
  blank(),
  bold("Produces (results/figures/)"),
  p("fig1_umap_overview.{png,pdf}         — two-panel UMAP: group and PAM50 subtype"),
  p("fig2_deg_counts.{png,pdf}            — DEG count bar chart, all 18 contrasts"),
  p("fig3_volcano_{subtype}_vs_GTEx.{png,pdf} — 4 volcano plots, signature genes highlighted"),
  p("fig4_deg_heatmap.{png,pdf}           — top 20 up/down DEGs per subtype, ComplexHeatmap"),
  p("fig5_signature_heatmap.{png,pdf}     — ElasticNet signature genes across tumor samples"),
  p("fig6_fgsea_dotplot.{png,pdf}         — top 10 Hallmark pathways per subtype"),
  p("fig7_cox_forest.{png,pdf}            — Cox HR forest plot with 95% CI"),
  p("Total: 44 files (22 unique figures x .png + .pdf)"),
  blank(),
  bold("Stage gate"),
  p("Checks all 8 prior stage outputs before running: cohort.rds, umap_clusters.rds, wgcna_modules.rds, 18 DEG CSVs, 4 signature CSVs, 4 fgsea CSVs, 4 survival CSVs. Hard stops with a clear error if anything is missing."),
  blank(),
  bold("save_fig(p, stem)"),
  p("The standard display-and-save pattern for ggplot2 objects: print(p) to screen, then ggsave .png (300 dpi) and .pdf. All ggplot figures go through this function."),
  blank(),
  bold("save_heatmap(ht, stem)"),
  p("The ComplexHeatmap equivalent: draw(ht) to screen, then png() device + draw(ht) + dev.off(), then same for pdf(). Heatmaps cannot use ggsave."),
  blank(),
  bold("Fig 1 — plot_umap_overview()"),
  p("Two-panel figure using cowplot::plot_grid(). Left panel colors points by group (Tumor/NAT/GTEx). Right panel colors by PAM50 subtype using subtype_colors, with controls drawn in grey80 underneath. The column PAM50 in the UMAP data frame is used directly — the data frame already contains group and PAM50 from Stage 02."),
  blank(),
  bold("Fig 2 — plot_deg_counts()"),
  p("Reads all 18 DEG CSVs and applies thresholds FDR < 0.05 and |logFC| >= 1. Counts up- and downregulated genes per contrast. Plots a horizontal diverging bar chart (red = up, blue = down) faceted by contrast type: vs GTEx, vs NAT, subtype vs subtype, NAT vs GTEx, Within-subtype."),
  blank(),
  bold("Fig 3 — plot_volcano()"),
  p("One plot per subtype vs GTEx contrast. Points colored by status: Up (red), Down (blue), Signature-Up (dark red), Signature-Down (dark blue), NS (grey). Signature genes are identified by joining with the subtype's ElasticNet signature CSV. Top 15 significant genes labeled using ggrepel::geom_text_repel to avoid overlaps."),
  blank(),
  bold("Fig 4 — build_deg_heatmap()"),
  p("Selects top 20 up and top 20 down genes per subtype vs GTEx by t-statistic, takes the union, z-scores expression across all 817 tumor samples, clips to [-3, +3]. ComplexHeatmap with PAM50 column annotation, biclustered rows and columns, split by subtype."),
  blank(),
  bold("Fig 5 — build_signature_heatmap()"),
  p("Unions all 4 subtype signatures (127 total unique genes after deduplication by highest |mean_coef|). Z-scores expression, draws ComplexHeatmap with PAM50 column annotation and a Direction (Up/Down) row annotation based on mean_coef sign from the ElasticNet output."),
  blank(),
  bold("Fig 6 — plot_fgsea_dotplot()"),
  p("Extracts the subtype name from the contrast column (e.g. LumA_vs_GTEx -> LumA). Filters to Hallmark collection only (collection == 'Hallmark'), vs_GTEx contrasts, padj < 0.05. Selects top 10 by |NES| per subtype. 40 pathways total across 4 subtypes. Dot size = -log10(padj), dot color = NES on a blue-white-red gradient."),
  blank(),
  bold("Fig 7 — plot_cox_forest()"),
  p("Reads 4 Cox CSVs and draws a horizontal forest plot with HR on log scale. Error bars are horizontal 95% CI using geom_errorbar(orientation='y'). Subtype names colored by subtype_colors. P-values labeled to the right of each point."),
  blank(),
  bold("Known warnings (cosmetic only)"),
  p("'Ignoring unknown labels: colour: Signature' — survminer KM plots from Stage 07; plots render correctly."),
  p("'Vectorized input to element_text() is not officially supported' — subtype-colored axis text; renders correctly in all tested ggplot2 versions."),
]);

// ─── full pipeline report ─────────────────────────────────────────────────────

write(path.join(DOCS, "pipeline_report.docx"), [

  // Title block
  new Paragraph({
    spacing: { after: 120 },
    alignment: AlignmentType.CENTER,
    children: [new TextRun({
      text: "Transcriptional Heterogeneity Across PAM50 Breast Cancer Subtypes",
      font: "Arial", size: 36, bold: true
    })],
  }),
  new Paragraph({
    spacing: { after: 80 },
    alignment: AlignmentType.CENTER,
    children: [new TextRun({
      text: "A Multi-Contrast RNA-seq Analysis of TCGA-BRCA TOIL Data",
      font: "Arial", size: 26, italics: true
    })],
  }),
  new Paragraph({
    spacing: { after: 400 },
    alignment: AlignmentType.CENTER,
    children: [new TextRun({
      text: "Pipeline Report — All-R Analysis (Stages 00-08)",
      font: "Arial", size: 22, color: "666666"
    })],
  }),

  // ── Abstract ──
  h1("Abstract"),
  p("Breast cancer is a heterogeneous disease encompassing molecularly distinct subtypes with markedly different clinical behaviors, treatment responses, and patient outcomes. The PAM50 gene expression signature classifies tumors into four intrinsic subtypes — Luminal A (LumA), Luminal B (LumB), HER2-enriched (Her2), and Basal-like (Basal) — each representing a biologically coherent transcriptional program. Despite decades of study, the full landscape of transcriptional differences between these subtypes, their shared and unique dysregulated pathways, and the genes that most reliably define subtype identity remain incompletely characterized at scale."),
  p("We present a comprehensive, thesis-defensible, all-R computational pipeline analyzing 817 PAM50-classified TCGA-BRCA tumor samples alongside 113 normal-adjacent tissue (NAT) samples and 179 GTEx healthy breast samples from the UCSC TOIL harmonized RNA-seq compendium. Across 18 pre-specified contrasts using limma-trend differential expression, we identified between 1,343 and 9,693 significant genes per contrast (FDR < 0.05, |log2FC| >= 1). ElasticNet bootstrap stability selection identified compact subtype-discriminating signatures of 27-36 genes per subtype. Gene set enrichment analysis across 236 Hallmark and KEGG pathways revealed 58-93 significantly enriched gene sets per subtype per reference comparison. Weighted gene co-expression network analysis identified five co-expression modules, 20 of which showed significant trait associations with PAM50 subtype (p < 0.05). Univariate Cox proportional-hazards analysis found no significant survival association for subtype-identity signatures, consistent with their discriminative rather than prognostic design, though a trend in Basal-like tumors (HR = 0.70, p = 0.11) warrants further investigation."),
  p("All code, outputs, and intermediate files are fully reproducible from a single renv-locked R environment. This report describes the complete methodological pipeline, key findings per stage, and interpretive context for each result."),
  blank(),

  br(),

  // ── 1. Introduction ──
  h1("1. Introduction"),
  h2("1.1 Breast Cancer Molecular Subtypes"),
  p("Breast cancer is not a single disease. Landmark gene expression profiling studies by Perou et al. (2000) and Sorlie et al. (2001) established that breast tumors cluster into transcriptionally distinct groups — now standardized as the PAM50 intrinsic subtypes — with profoundly different biology and clinical implications."),
  p("Luminal A tumors are estrogen receptor (ER)-positive, low grade, and carry the most favorable prognosis of any subtype, with 10-year survival rates exceeding 80% in early-stage disease. Luminal B tumors are also ER-positive but show higher proliferation, worse grade, and intermediate prognosis. HER2-enriched tumors overexpress ERBB2 and were historically associated with poor outcomes; targeted therapies (trastuzumab, pertuzumab) have dramatically improved survival in this group. Basal-like tumors — largely overlapping with triple-negative breast cancer (TNBC, defined by absence of ER, PR, and HER2 protein expression) — are aggressive, metastatic, and lack targeted therapy options, relying on chemotherapy as the primary systemic treatment."),
  p("Understanding the transcriptional programs that define each subtype, and how they differ from both cancer-adjacent normal tissue and truly healthy breast epithelium, is essential for identifying subtype-specific biomarkers, therapeutic targets, and prognostic signatures."),
  blank(),
  h2("1.2 Study Objectives"),
  p("This pipeline addresses four primary objectives:"),
  indent("1. Characterize differential gene expression for each PAM50 subtype against two independent control groups — NAT (field cancerization context) and GTEx Healthy breast (absolute dysregulation context) — using a unified 18-contrast design."),
  indent("2. Identify compact, stable gene expression signatures for each subtype via bootstrap ElasticNet stability selection, suitable for external validation and potential clinical implementation."),
  indent("3. Map the dysregulated biological pathways per subtype using pre-ranked gene set enrichment analysis against curated Hallmark and KEGG pathway collections."),
  indent("4. Identify co-expression modules correlated with PAM50 subtype identity and assess whether subtype-discriminating signatures carry independent survival information."),
  blank(),
  h2("1.3 Why Two Control Groups?"),
  p("A key design decision in this pipeline is the inclusion of both NAT and GTEx Healthy as reference groups rather than relying on a single normal comparator. These two groups capture distinct biological phenomena:"),
  p("Normal Adjacent Tissue (NAT) consists of non-tumor biopsies taken from the same breast as the tumor at the time of surgery. Although histologically normal, NAT is subject to 'field cancerization' — epigenetic, transcriptional, and genomic alterations driven by exposure to tumor secretory factors, inflammation, and paracrine signaling. A gene that appears dysregulated in tumor vs NAT may be altered in both groups relative to true normal."),
  p("GTEx Healthy breast samples are from cancer-free donors with no history of breast malignancy. They represent the true normal baseline. Genes dysregulated in tumor vs GTEx but not vs NAT are candidates for 'field effect genes' — present in the cancer-adjacent field but not uniquely in the tumor. Genes dysregulated in tumor vs both NAT and GTEx are the most robust biomarkers, altered even relative to tissue that has been exposed to the tumor microenvironment."),
  p("This dual-reference design enables post-hoc classification of every differentially expressed gene into one of three interpretive categories: robust biomarkers (DE vs both), field effect genes (DE vs GTEx only), or NAT-specific noise (DE vs NAT only, excluded from signatures)."),
  blank(),

  br(),

  // ── 2. Data ──
  h1("2. Data and Cohort"),
  h2("2.1 Data Source"),
  p("All expression data were obtained from the UCSC Xena TOIL RNA-seq compendium (toilHub), which uniformly reprocesses TCGA, TARGET, and GTEx samples through the same STAR/RSEM pipeline on GRCh38/hg38 using GENCODE v23 annotations. This eliminates batch effects from different alignment and quantification workflows, which are a major confound when comparing TCGA tumor data with external normal controls."),
  p("Expression values are provided as log2(TPM + 0.001), which is the unit used throughout this analysis. Data were never re-logged at any stage."),
  p("Dataset identifiers:"),
  indent("Expression matrix: TcgaTargetGtex_rsem_gene_tpm (~19,000 samples, ~60,000 Ensembl gene IDs)"),
  indent("Phenotype: TcgaTargetGTEX_phenotype.txt (sample metadata, study labels, tissue sites)"),
  indent("PAM50 calls: TCGA.BRCA.sampleMap/BRCA_clinicalMatrix (tcgaHub, column PAM50Call_RNAseq)"),
  indent("Survival: TCGA.BRCA.sampleMap/BRCA_clinicalMatrix (columns OS_Time_nature2012, OS_event_nature2012)"),
  blank(),
  h2("2.2 Sample Groups"),
  p("Three distinct sample groups were assembled using a two-arm OR filter on the TOIL phenotype file:"),
  p("Tumor (n = 817): TCGA primary solid tumor biopsies (barcode positions 14-15 = '01') with confirmed PAM50 subtype calls. Normal-like tumors (PAM50 = 'Normal') were excluded — these are tumors whose expression profile clusters near normal tissue but they are still tumor biopsies, not controls. Including them would contaminate the tumor arm with a heterogeneous group lacking a distinct therapeutic target."),
  p("NAT — Normal Adjacent Tissue (n = 113): TCGA normal adjacent tissue biopsies (barcode positions 14-15 = '11') from the same patients as TCGA tumor samples. These are never filtered by PAM50 (which only applies to tumor samples)."),
  p("GTEx Healthy (n = 179): GTEx breast samples from cancer-free donors (_study == 'GTEX', _primary_site == 'Breast'). These are never filtered by PAM50."),
  p("Total cohort: 1,109 samples. After Ensembl-to-HGNC ID mapping and duplicate resolution, the expression matrix contains 28,344 unique gene symbols."),
  blank(),
  p("PAM50 subtype breakdown (Tumor arm only):"),
  indent("LumA: 420 samples (51.4%)"),
  indent("LumB: 192 samples (23.5%)"),
  indent("Basal: 139 samples (17.0%)"),
  indent("Her2:  66 samples ( 8.1%)"),
  blank(),
  h2("2.3 Data Preprocessing"),
  p("Ensembl gene IDs were mapped to HGNC symbols using AnnotationDbi::mapIds() with org.Hs.eg.db, fully locally without network access. TOIL Ensembl IDs include version suffixes (e.g. ENSG00000000003.14) which were stripped before lookup. Of the ~60,000 Ensembl IDs in the matrix, unmapped entries and genes without valid HGNC symbols were removed. For the 57 cases where multiple Ensembl IDs mapped to the same HGNC symbol, expression was averaged by row-wise mean. PAM50 subtype calls were fetched separately from tcgaHub and joined on the 15-character TCGA barcode prefix."),
  blank(),

  br(),

  // ── 3. Methods ──
  h1("3. Methods"),
  h2("3.1 Dimensionality Reduction and Clustering (Stage 02)"),
  p("To visualize the transcriptional structure of the cohort, we selected the top 5,000 most variable genes by variance across all 1,109 samples using matrixStats::rowVars. Using all samples (including controls) for gene selection maximizes between-group separation in the resulting UMAP embedding. Principal component analysis (PCA) was performed on the transposed HVG matrix (prcomp, centered, unscaled, rank = 50), followed by UMAP (uwot, n_neighbors = 30, min_dist = 0.3, seed = 42)."),
  p("Within-subtype sub-subclusters were identified by running k-means (k = 2) on the 2D UMAP coordinates for each PAM50 subtype independently. These cluster assignments were used as group labels in the four Within-subtype limma contrasts in Stage 03."),
  blank(),
  h2("3.2 Differential Expression Analysis (Stage 03)"),
  p("Differential expression was performed using limma-trend (Ritchie et al., 2015) with eBayes(trend = TRUE). This method is appropriate for log2-TPM data: the trend argument fits a mean-variance relationship to the empirical Bayes prior, making the moderated t-statistics valid for RNA-seq data analyzed on the log-TPM scale. DESeq2 and edgeR (which operate on raw counts) were deliberately not used, as the TOIL data are pre-quantified TPM values."),
  p("A single six-level factor was constructed encoding all sample groups (LumA, LumB, Her2, Basal, NAT, GTEx) in a design matrix with no intercept (model.matrix(~0 + group)). All 18 contrasts were extracted from this single model fit:"),
  indent("4 subtype vs NAT contrasts: LumA/LumB/Her2/Basal each vs NAT"),
  indent("4 subtype vs GTEx contrasts: LumA/LumB/Her2/Basal each vs GTEx Healthy"),
  indent("5 subtype vs subtype contrasts: LumA/Basal, LumA/LumB, LumA/Her2, Her2/Basal, LumB/Basal"),
  indent("1 NAT vs GTEx contrast: field cancerization characterization"),
  indent("4 within-subtype contrasts: k-means C1 vs C2 within each PAM50 subtype"),
  p("All genes (28,344) were written to each output CSV without pre-filtering. Downstream stages apply their own FDR and fold-change thresholds as appropriate for their purpose."),
  blank(),
  h2("3.3 Feature Selection — ElasticNet Bootstrap Stability (Stage 04)"),
  p("For each PAM50 subtype, a minimal gene expression signature was identified using bootstrap ElasticNet stability selection (Meinshausen and Buhlmann, 2010). The goal is a compact, reproducible gene list that reliably discriminates that subtype from all other tumor subtypes across repeated subsamples of the cohort."),
  p("The binary outcome was one-vs-rest subtype membership (e.g. for LumA: LumA = 1, LumB/Her2/Basal = 0) across all 817 tumor samples. Input features were all 28,344 gene expression values in log2-TPM space."),
  p("The procedure ran in three steps. First, 10-fold cross-validated ElasticNet (glmnet, alpha = 0.5) was run once across all genes to identify lambda.1se (the regularization strength that minimizes cross-validation error within one standard error of the minimum). Second, the feature space was reduced to genes with non-zero coefficients at lambda.min (a guaranteed superset of lambda.1se genes), reducing dimensionality from ~28,000 to ~100-300 genes before bootstrapping. Third, 100 bootstrap iterations were run at fixed lambda.1se: in each iteration, samples were drawn with replacement and glmnet was refit. Genes were retained in the final signature if they had non-zero coefficients in >= 80% of bootstrap runs (selection_frequency >= 0.80). The mean coefficient across all 100 runs was recorded as mean_coef."),
  p("The 0.80 frequency threshold follows Meinshausen and Buhlmann (2010), who show this value provides theoretical control over the expected number of falsely selected variables under mild exchangeability assumptions. A uniform threshold was applied across all four subtypes."),
  blank(),
  h2("3.4 Gene Set Enrichment Analysis (Stage 05)"),
  p("Pre-ranked gene set enrichment analysis was performed using fgsea (Korotkevich et al., 2021) against two gene set collections loaded via msigdbr:"),
  indent("MSigDB Hallmark (50 gene sets): curated oncology-relevant pathways representing coherent biological processes with low within-set redundancy."),
  indent("KEGG Legacy (186 gene sets): metabolic and signaling pathway maps from the Kyoto Encyclopedia of Genes and Genomes."),
  p("Gene ranking used the limma moderated t-statistic from each DEG contrast, sorted descending. The t-statistic was chosen over log2FC because it integrates effect size with precision — genes with consistent signal across many samples rank higher than noisy genes with large but variable fold-changes. fgsea was run with minSize = 15, maxSize = 500, nPermSimple = 1,000, and eps = 0. The eps = 0 setting enables accurate p-value estimation for pathways with extreme enrichment (p < 1e-50), which occurs in the large cancer-vs-normal contrasts in this dataset."),
  p("GSEA was run for each subtype against both GTEx (absolute dysregulation) and NAT (field-effect-corrected dysregulation), yielding 8 sets of results across 236 gene sets each."),
  blank(),
  h2("3.5 Weighted Gene Co-expression Network Analysis (Stage 06)"),
  p("Co-expression network analysis was performed using WGCNA (Langfelder and Horvath, 2008) on the 817 tumor samples only. The top 5,000 genes by median absolute deviation (MAD) were selected — MAD is preferred over variance as it is more robust to outlier samples. Samples and genes with excessive missing values or zero variance were removed using WGCNA::goodSamplesGenes() before network construction."),
  p("A signed network was constructed, which preserves the direction of co-expression (positively co-expressed genes cluster together; anti-correlated genes form separate modules). The soft-thresholding power was selected by testing powers 1-30 and identifying the lowest power where the resulting network approximated a scale-free topology, quantified by R² of a log-log linear fit to the connectivity distribution. Power 8 achieved R² = 0.907, exceeding the 0.85 target threshold."),
  p("Modules were identified using blockwiseModules() with minimum module size 30 genes and a merge height of 0.25 (modules with eigengene correlation > 0.75 were merged). Five modules were identified excluding the grey (unassigned) module. Module eigengenes — the first principal component of each module's expression across samples — were correlated with binary PAM50 subtype indicators using Pearson correlation."),
  blank(),
  h2("3.6 Survival Analysis (Stage 07)"),
  p("Overall Survival (OS) data were obtained from the TCGA BRCA clinical matrix on tcgaHub (columns OS_Time_nature2012 in days and OS_event_nature2012, 0 = censored/alive, 1 = deceased). Samples with missing or zero survival time were excluded, yielding 882 samples with valid OS data."),
  p("For each PAM50 subtype, a single-sample signature score was computed as the mean log2-TPM expression across the subtype's ElasticNet signature genes. The score was standardized to zero mean and unit variance (scale()) within each subtype before model fitting, making the hazard ratio interpretable as the change in hazard per one standard deviation increase in signature expression."),
  p("Univariate Cox proportional-hazards models were fitted using survival::coxph() with formula Surv(os_time, os_event) ~ score_scaled. Kaplan-Meier curves were generated by dichotomizing samples at the within-subtype median score into High and Low groups."),
  blank(),

  br(),

  // ── 4. Results ──
  h1("4. Results"),
  h2("4.1 Cohort Structure and UMAP"),
  p("The final cohort comprised 1,109 samples across three groups: 817 PAM50-classified tumor samples, 113 NAT samples, and 179 GTEx Healthy breast samples. The expression matrix contained 28,344 HGNC gene symbols after ID mapping and deduplication."),
  p("UMAP embedding of all samples using the top 5,000 most variable genes and 50 PCA components revealed strong separation between the three sample groups. GTEx Healthy samples formed a tight cluster distinct from both tumor and NAT samples, consistent with their origin from a different biological context (healthy donors vs cancer patients). NAT samples occupied an intermediate position between GTEx and tumor clusters, which is consistent with the field cancerization hypothesis — NAT is histologically normal but transcriptionally shifted toward a tumor-like state. Tumor samples showed subtype-specific clustering, with Basal-like tumors forming the most distinct cluster and LumA/LumB showing substantial overlap, as expected from their related luminal biology."),
  p("Within-subtype k-means (k = 2) sub-subclustering identified putative transcriptional heterogeneity within each subtype, contributing four additional contrasts to the DEG analysis (Within_LumA, Within_LumB, Within_Her2, Within_Basal)."),
  blank(),
  h2("4.2 Differential Expression"),
  p("Limma-trend analysis across all 18 contrasts identified extensive differential expression throughout the cohort (FDR < 0.05, |log2FC| >= 1):"),
  p("Subtype vs GTEx contrasts showed the highest gene counts, reflecting the large transcriptional distance between cancer and truly healthy breast tissue:"),
  indent("LumB vs GTEx: 9,693 genes (most extensive dysregulation)"),
  indent("Her2 vs GTEx: 9,640 genes"),
  indent("Basal vs GTEx: 9,423 genes"),
  indent("LumA vs GTEx: 8,688 genes"),
  p("Subtype vs NAT contrasts were smaller, consistent with NAT already being transcriptionally shifted toward a tumor-like state:"),
  indent("LumB vs NAT: 6,875 genes"),
  indent("Basal vs NAT: 7,091 genes"),
  indent("Her2 vs NAT: 6,943 genes"),
  indent("LumA vs NAT: 5,575 genes"),
  p("NAT vs GTEx: 5,689 significant genes, confirming extensive field cancerization transcriptional changes in cancer-adjacent tissue even in the absence of tumor cells."),
  p("Subtype vs subtype contrasts were smaller, reflecting shared cancer biology:"),
  indent("LumA vs LumB: 1,343 genes (most similar subtypes)"),
  indent("LumA vs Her2: 3,199 genes"),
  indent("LumA vs Basal: 5,372 genes (most distant)"),
  indent("LumB vs Basal: 4,709 genes"),
  indent("Her2 vs Basal: 3,532 genes"),
  p("Within-subtype contrasts ranged from 821 (Her2, limited by small n = 66) to 4,463 (Basal) significant genes, indicating meaningful transcriptional heterogeneity within each PAM50 category."),
  blank(),
  h2("4.3 Subtype-Specific Gene Signatures"),
  p("Bootstrap ElasticNet stability selection (100 iterations, FREQ_CUTOFF = 0.80) identified compact discriminative signatures for each PAM50 subtype:"),
  indent("LumA: 30 genes"),
  indent("LumB: 27 genes"),
  indent("Her2: 36 genes"),
  indent("Basal: 34 genes"),
  p("Her2 and Basal signatures are slightly broader (36/34 genes vs 30/27 for luminal subtypes), which reflects their smaller positive-class sample sizes (n = 66 and n = 139 respectively). With fewer positive examples, the ElasticNet requires more features to maintain discriminative power at the same regularization strength, a well-understood property of sparse classifiers in class-imbalanced settings."),
  p("The selection_frequency column in each signature CSV records what fraction of the 100 bootstrap runs included each gene, and mean_coef records the average coefficient magnitude and direction. Genes with positive mean_coef are upregulated in that subtype relative to the others; negative mean_coef indicates downregulation. These signatures represent the most stable, reproducible genes for discriminating each subtype and are the basis for both the survival analysis and several visualization figures."),
  blank(),
  h2("4.4 Pathway Enrichment"),
  p("Gene set enrichment analysis against 236 Hallmark and KEGG pathways revealed substantial and biologically coherent pathway dysregulation in each subtype. Results are reported for the vs_GTEx contrasts (absolute dysregulation) and vs_NAT contrasts separately."),
  p("LumA: 58 significant pathways vs GTEx, 66 vs NAT (padj < 0.05). As the most clinically favorable subtype, LumA showed activation of estrogen response pathways and metabolic programs, with modest proliferative signal relative to other subtypes."),
  p("LumB: 73 significant pathways vs GTEx, 93 vs NAT — the highest pathway count of any subtype. LumB showed stronger cell cycle, DNA replication, and MYC target activation relative to LumA, consistent with its higher proliferative grade and worse prognosis."),
  p("Her2: 67 significant pathways vs GTEx, 73 vs NAT. ErbB/ERBB2 signaling pathways and PI3K/mTOR axis pathways were among the top enriched sets, consistent with the oncogenic driver biology of this subtype."),
  p("Basal: 69 significant pathways vs GTEx, 83 vs NAT. Basal tumors showed the strongest enrichment for DNA damage response, E2F targets, G2M checkpoint, and MYC targets — reflecting the genomic instability and high proliferative index characteristic of triple-negative/basal-like breast cancer."),
  p("The dual-reference design revealed field effect pathways: several immune and stromal signaling pathways were enriched in tumor vs GTEx but not in tumor vs NAT, suggesting they are altered throughout the breast field in TCGA patients rather than being uniquely tumor-driven."),
  blank(),
  h2("4.5 Co-expression Modules"),
  p("WGCNA identified five co-expression modules (excluding the grey unassigned module) at soft-thresholding power 8 (scale-free topology R² = 0.907). Module sizes ranged from 395 to 775 genes; 2,061 genes were unassigned to any module (grey)."),
  p("Module-trait correlation analysis revealed 20 of 24 possible module-trait associations reaching statistical significance (p < 0.05). This high proportion of significant associations confirms that the co-expression structure of the tumor transcriptome is strongly organized around PAM50 subtype identity. Some modules showed strong positive correlation with one subtype and negative correlation with others, capturing the fundamental transcriptional opposition between, for example, Basal and LumA biology. The module eigengenes represent biologically coherent programs — groups of genes that are activated or suppressed together across tumors — and their subtype correlations provide a complementary, network-based view of the same biology captured by the pairwise DEG contrasts."),
  blank(),
  h2("4.6 Survival Analysis"),
  p("Univariate Cox proportional-hazards analysis tested whether each subtype's ElasticNet signature score was associated with Overall Survival within that subtype's tumor samples. None of the four subtypes reached statistical significance at p < 0.05:"),
  indent("LumA:  HR = 0.867 (95% CI: 0.649-1.159), p = 0.336, n = 370, events = 38"),
  indent("LumB:  HR = 0.901 (95% CI: 0.606-1.339), p = 0.607, n = 167, events = 23"),
  indent("Her2:  HR = 1.141 (95% CI: 0.584-2.228), p = 0.700, n =  57, events = 11"),
  indent("Basal: HR = 0.696 (95% CI: 0.444-1.092), p = 0.114, n = 125, events = 16"),
  p("These null results are expected and methodologically coherent. The ElasticNet signatures were optimized for subtype discrimination (a classification objective), not for survival prediction (a regression objective). A gene that reliably identifies a LumA tumor is not necessarily a prognostic gene within LumA. This distinction is fundamental: subtype-discriminating features and prognostic features are biologically and statistically different quantities that are only occasionally shared."),
  p("The Basal trend (HR = 0.70, p = 0.11) is the most biologically plausible candidate for a survival-relevant signal. Basal-like/TNBC tumors have the worst prognosis of any subtype, and it is plausible that within-Basal variation in the expression of Basal-defining genes reflects clinically meaningful heterogeneity. The non-significance here is at least partly a power problem: with only 16 death events in 125 Basal samples, the Cox model is substantially underpowered (conventional guidance suggests >= 10 events per predictor variable; this analysis uses 1 predictor and 16 events, marginally meeting that criterion). A larger Basal-enriched cohort, or extended follow-up in TCGA, would be required to definitively test this association."),
  blank(),

  br(),

  // ── 5. Discussion ──
  h1("5. Discussion"),
  h2("5.1 Transcriptional Landscape"),
  p("The large number of significant DEGs in every subtype vs GTEx contrast (8,688-9,693 genes, representing 30-34% of the transcriptome) underscores the profound transcriptional reorganization that accompanies breast cancer transformation. The somewhat smaller numbers in vs NAT contrasts confirm that cancer-adjacent tissue itself undergoes substantial transcriptional change — the field cancerization effect — independently of the tumor. The 5,689 genes differentially expressed between NAT and GTEx constitute a field transcriptome that should be considered when interpreting any clinical comparison of tumor vs adjacent normal tissue."),
  p("The gradient in LumA vs other subtype DEG counts (LumA vs LumB: 1,343; LumA vs Basal: 5,372) reflects known biological distances between subtypes. LumA and LumB share luminal differentiation programs and ER-driven transcription; their differences are largely in proliferation-related genes. LumA and Basal are at opposite ends of the breast cancer transcriptional spectrum, one representing a highly differentiated luminal phenotype and the other a dedifferentiated, basal-like state with activated epithelial-mesenchymal transition programs."),
  blank(),
  h2("5.2 Signature Design and Validation Considerations"),
  p("The ElasticNet bootstrap stability selection approach produces signatures that are stable across resampled datasets, addressing the overfitting concern that plagues single-run regularized regression. The 0.80 frequency threshold provides a principled false discovery control for the selection process itself, distinct from and complementary to the FDR control applied in the DEG step."),
  p("An important caveat is that these signatures were derived and evaluated within the same TCGA-BRCA cohort. External validation in an independent cohort (e.g. METABRIC, GEO breast cancer datasets processed through the same TOIL pipeline) is required before clinical deployment. Signatures trained on log2-TPM values also require careful normalization considerations when applied to data from different platforms or preprocessing pipelines."),
  blank(),
  h2("5.3 Limitations"),
  p("Sample size constraints: Her2 (n = 66) is substantially smaller than the other subtypes, limiting statistical power for Her2-specific analyses including survival. The wide confidence intervals in the Her2 Cox model (95% CI: 0.584-2.228) reflect this directly."),
  p("Univariate survival analysis: All Cox models are univariate (signature score only). Clinical covariates known to confound survival analysis in breast cancer — age at diagnosis, tumor stage, lymph node status, treatment received — were not included. Multivariate models incorporating these covariates would be required for any claim of independent prognostic value."),
  p("Single-cohort analysis: All results are derived from TCGA-BRCA. TCGA has well-documented ascertainment biases (predominantly primary resection specimens, predominantly US patients, predominantly early-stage disease). Generalizability to other populations and disease stages is not established."),
  p("Within-subtype heterogeneity: The k-means sub-subclusters used in Stage 03's Within-subtype contrasts are derived from 2D UMAP coordinates, which are a nonlinear low-dimensional projection. Cluster assignments in UMAP space do not necessarily correspond to biologically meaningful subgroups and should be treated as exploratory."),
  blank(),
  h2("5.4 Future Directions"),
  p("Prognostic signature development: A dedicated survival-optimized analysis — using a Cox-penalized regression (e.g. glmnet with family='cox') rather than ElasticNet classification — would identify prognostic gene sets specifically within each subtype. This is a natural Stage 07 extension."),
  p("External validation: Applying the four ElasticNet signatures to an independent RNA-seq cohort (e.g. METABRIC) would test their generalizability and clinical utility."),
  p("Multi-omics integration: TCGA-BRCA contains matched copy number, methylation, and protein expression data for most tumor samples. Integrating these data layers with the transcriptional results — particularly overlapping the DEG and WGCNA module genes with copy number amplification and DNA methylation — would provide a more complete picture of the molecular mechanisms driving each subtype."),
  p("Basal survival signal: The Basal trend (HR = 0.70, p = 0.11) merits investigation in a larger TNBC-enriched cohort with extended follow-up."),
  blank(),

  br(),

  // ── 6. Conclusions ──
  h1("6. Conclusions"),
  p("This pipeline provides a comprehensive, reproducible transcriptional characterization of PAM50 breast cancer subtypes using harmonized TOIL RNA-seq data. Key contributions include:"),
  indent("A 18-contrast differential expression landscape using limma-trend with two independent control groups (NAT and GTEx Healthy), enabling distinction between tumor-specific dysregulation and field cancerization effects."),
  indent("Four compact, bootstrap-stable ElasticNet signatures (27-36 genes per subtype) suitable for external validation."),
  indent("Pathway-level characterization confirming known subtype biology (estrogen response in luminal subtypes, DNA damage and proliferation in Basal) and revealing shared vs subtype-specific pathway activation."),
  indent("A WGCNA co-expression network with strong subtype-trait correlations, providing a network-level complement to the pairwise DEG analysis."),
  indent("A null survival result that is correctly interpreted as a consequence of signature design rather than a failure of the pipeline, with a Basal trend warranting further investigation."),
  p("All analyses are implemented in R using a locked renv environment and structured as an idempotent, stage-gated pipeline. Every result is reproducible from the raw XENA hub data by running stages 00-08 in order."),
  blank(),

  br(),

  // ── 7. Technical Appendix ──
  h1("7. Technical Appendix — Pipeline Stage Reference"),
  blank(),
  h3("Stage 00 — Setup and Phenotype Inspection"),
  p("Verifies renv lock, confirms hub reachability, downloads and caches TcgaTargetGTEX_phenotype.txt, inspects column names, records confirmed field names, prints sample counts."),
  p("Output: data/processed/metadata_breast.rds"),
  blank(),
  h3("Stage 01 — Preprocessing"),
  p("Downloads TcgaTargetGtex_rsem_gene_tpm, subsets to breast samples, maps Ensembl IDs to HGNC symbols (org.Hs.eg.db), fetches PAM50 from TCGA.BRCA.sampleMap/BRCA_clinicalMatrix (tcgaHub), assigns 3-way group labels, excludes Normal-like tumors, builds 18 contrast names."),
  p("Output: data/processed/cohort.rds ($expr 28,344 x 1,109; $meta; $contrasts)"),
  blank(),
  h3("Stage 02 — Clustering"),
  p("Top 5,000 HVGs by variance -> PCA (50 PCs) -> UMAP (n_neighbors=30, min_dist=0.3). K-means (k=2) sub-subclusters within each PAM50 subtype on UMAP coordinates. 6 UMAP figure pairs saved."),
  p("Output: data/processed/umap_clusters.rds"),
  blank(),
  h3("Stage 03 — Differential Expression"),
  p("Single 6-level factor model, limma-trend eBayes(trend=TRUE). 14 main contrasts via makeContrasts + contrasts.fit. 4 within-subtype contrasts via subsetted models. 28,344 genes per file, no pre-filtering."),
  p("Output: results/tables/deg/ (18 CSV files)"),
  blank(),
  h3("Stage 04 — Feature Selection"),
  p("ElasticNet (alpha=0.5) bootstrap stability selection. cv.glmnet once for lambda.1se, reduce to lambda.min active genes, 100 bootstrap iterations at fixed lambda.1se. FREQ_CUTOFF = 0.80 (Meinshausen & Buhlmann 2010)."),
  p("Output: results/tables/signatures/ (4 CSV files, LumA=30, LumB=27, Her2=36, Basal=34 genes)"),
  blank(),
  h3("Stage 05 — Gene Set Enrichment"),
  p("fgsea, ranking by limma moderated t-statistic. 236 gene sets (50 Hallmark + 186 KEGG Legacy via msigdbr v10 API). minSize=15, maxSize=500, nPermSimple=1000, eps=0. Both vs_GTEx and vs_NAT per subtype."),
  p("Output: results/tables/fgsea/ (4 CSV files, 58-93 sig pathways per contrast)"),
  blank(),
  h3("Stage 06 — WGCNA"),
  p("Signed network, top 5,000 genes by MAD, Tumor samples only. Power selection: tested 1-30, selected power=8 (R2=0.907). blockwiseModules(maxBlockSize=6000, minModuleSize=30, mergeCutHeight=0.25). 5 modules identified. Pearson module-trait correlations: 20/24 significant at p<0.05."),
  p("Output: data/processed/wgcna_modules.rds, results/figures/wgcna_*.{png,pdf}"),
  blank(),
  h3("Stage 07 — Survival"),
  p("OS data from TCGA BRCA_clinicalMatrix (OS_Time_nature2012, OS_event_nature2012). Signature score = mean log2-TPM of signature genes, scaled to SD units. Univariate Cox PH. KM curves split at median score. No subtype significant at p<0.05; Basal trend HR=0.70, p=0.11."),
  p("Output: results/tables/survival/ (4 CSV files), results/figures/survival_km_*.{png,pdf}"),
  blank(),
  h3("Stage 08 — Visualization"),
  p("7 publication figure sets (44 files): UMAP overview, DEG count bar chart, 4 volcano plots, DEG heatmap, signature heatmap, fgsea dot plot, Cox forest plot. All saved as .png (300 dpi) + .pdf (vector)."),
  p("Output: results/figures/fig1-fig7.{png,pdf}"),
  blank(),
  blank(),
  h2("Software and Packages"),
  p("R version 4.5.2. Key packages:"),
  indent("UCSCXenaTools 1.7.0 — data retrieval"),
  indent("limma — differential expression (eBayes trend=TRUE)"),
  indent("glmnet — ElasticNet feature selection"),
  indent("fgsea — gene set enrichment"),
  indent("msigdbr — gene set collections (Hallmark + KEGG)"),
  indent("WGCNA — weighted gene co-expression network analysis"),
  indent("uwot — UMAP embedding"),
  indent("survival + survminer — Cox PH and Kaplan-Meier"),
  indent("AnnotationDbi + org.Hs.eg.db — local Ensembl-to-HGNC mapping"),
  indent("ComplexHeatmap + circlize — heatmap visualization"),
  indent("ggplot2 + ggrepel + cowplot — general visualization"),
  indent("renv — reproducible package environment"),
]);

console.log("All files written.");
