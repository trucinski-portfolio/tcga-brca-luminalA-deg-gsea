# CLAUDE.md — TCGA-BRCA Multi-Subtype Pipeline (v3.0)

> Living document. After every correction: "Update CLAUDE.md so you don't make that mistake again."
> v1.0 was architecturally clean. v3.0 must preserve that discipline and expand it.

---

## Project in One Sentence

Characterize transcriptional heterogeneity across PAM50 breast cancer subtypes (Luminal A,
Luminal B, HER2-enriched, Basal-like) using TCGA-BRCA RNA-seq — identifying subtype-specific
DEGs, co-expression modules, minimal biomarker signatures, and survival associations —
to produce a publishable-quality, thesis-defensible, all-R pipeline.

---

## Language Contract

**This is an all-R project. No Python.**
- Data acquisition: `TCGAbiolinks` or UCSC Xena downloads via R
- Preprocessing: `dplyr`, `tidyr`, `SummarizedExperiment`
- DEG: `limma-trend` (keep consistency with v1.0 method)
- GSEA: `fgsea` — NOT the external Broad GSEA .jar. That manual step is eliminated in v3.0.
- Feature selection: `glmnet` (ElasticNet/LASSO)
- Co-expression: `WGCNA`
- Survival: `survival` + `survminer`
- Visualization: `ggplot2`, `ComplexHeatmap`, `ggsurvplot`
- Environment: `renv` (lockfile committed, `.Rprofile` bootstraps renv automatically)

---

## Repo Layout (v3.0)

```
project-root/
├── CLAUDE.md
├── README.md
├── renv.lock                         ← committed, locked R environment
├── .Rprofile                         ← bootstraps renv on session start
├── .gitignore                        ← data/raw/ and data/processed/ excluded
│
├── R/                                ← all analysis scripts, sourced in order
│   ├── 00_setup.R                    ← renv check, dir creation, download validation
│   ├── 01_preprocessing.R            ← PAM50 assignment, cohort construction, QC
│   ├── 02_clustering.R               ← within-subtype UMAP + sub-subcluster detection
│   ├── 03_deg.R                      ← limma-trend DEG per subtype contrast
│   ├── 04_feature_selection.R        ← ElasticNet/LASSO → minimal biomarker signature
│   ├── 05_fgsea.R                    ← pathway enrichment (Hallmark + KEGG, fully in R)
│   ├── 06_wgcna.R                    ← co-expression module analysis
│   ├── 07_survival.R                 ← Cox regression on top DEG signatures
│   └── 08_visualization.R            ← publication-ready composite figures
│
├── data/
│   ├── raw/                          ← NOT committed
│   └── processed/                    ← .rds intermediates, NOT committed
│
├── results/
│   ├── figures/                      ← .png + .pdf outputs (committed)
│   └── tables/                       ← .csv outputs (committed)
│
└── .claude/
    └── commands/                     ← slash commands
```

---

## Data Contract

| Source | Object | Key field |
|--------|--------|-----------|
| TCGA-BRCA RNA-seq (Xena TOIL TPM) | expression matrix | TCGA barcode (col) |
| TCGA-BRCA clinical matrix (Xena) | phenotype data | `sampleID` |
| PAM50 assignments | Xena `PAM50Call_RNAseq` field | `sampleID` |

**Rules:**
- PAM50 subtypes in scope: `LumA`, `LumB`, `Her2`, `Basal`. Normal-like excluded.
- Expression values: log2(TPM + 0.001) — consistent with TOIL pipeline.
- Join always on TCGA barcode, trimmed to 15 characters (sample-level, not aliquot).
- Normal adjacent tissue samples kept as reference group (barcode position 14 = `11`).
- Do not mix GDC raw counts with Xena TOIL TPM — one source only, stay consistent.

---

## Pipeline Stages (in order, stage-gated)

```
00_setup.R             →  dirs created, raw data downloaded + validated
01_preprocessing.R     →  brca_se.rds  (SummarizedExperiment, PAM50 labeled, QC passed)
02_clustering.R        →  umap_clusters.rds  (within-subtype sub-subcluster assignments)
03_deg.R               →  results/tables/deg/  (one .csv per contrast)
04_feature_selection.R →  results/tables/signatures/  (ElasticNet gene sets per subtype)
05_fgsea.R             →  results/tables/fgsea/  (NES tables per subtype)
06_wgcna.R             →  data/processed/wgcna_modules.rds
07_survival.R          →  results/tables/survival/  (Cox HR tables)
08_visualization.R     →  results/figures/  (all publication figures)
```

Never skip stages. Never run stage N before stage N-1 output exists.
Each script checks for its required input at the top and stops with a clear error if missing.

---

## Subtype Contrast Design

Three contrast layers — document which are active in each script header:

1. **Subtype vs Normal** — each PAM50 subtype vs adjacent normal tissue
2. **Subtype vs Subtype** — pairwise where clinically meaningful (LumA vs Basal, etc.)
3. **Within-subtype sub-clustering** — UMAP on subtype-specific expression →
   Leiden or k-means → sub-subcluster labels → DEG within subtype

Do not add contrasts without updating this file first.

---

## Documentation Standard — Roxygen2 on every function

```r
#' Short title (verb phrase)
#'
#' One paragraph: what it does, when to use it, key assumptions.
#'
#' @param se SummarizedExperiment. Must have colData columns: PAM50, sample_type.
#' @param subtype Character. One of "LumA", "LumB", "Her2", "Basal".
#' @param fdr_cutoff Numeric. FDR threshold for DEG filtering. Default 0.05.
#'
#' @return data.frame with columns: gene, logFC, AveExpr, t, P.Value, adj.P.Val, B.
#'
#' @examples
#' deg <- run_limma_contrast(se = brca_se, subtype = "LumA", fdr_cutoff = 0.05)
run_limma_contrast <- function(se, subtype, fdr_cutoff = 0.05) { ... }
```

Every exported function gets full roxygen2. No stubs. No TODO comments left in.

---

## Coding Conventions

- `here::here()` for all paths. No hardcoded strings.
- `set.seed(42)` at the top of any script using random processes.
- Scripts are **idempotent** — re-running produces same output, never appends.
- Intermediate objects saved as `.rds`. Final tables as `.csv` in `results/tables/`.
- One script = one stage. Do not combine stages to "save files."
- `message()` at the start and end of every major function block.
- No `library()` calls inside functions — all at top of script.
- No `attach()`. No `<<-`.

---

## Subtype Color Palette — Fixed Across All Figures

```r
subtype_colors <- c(
  LumA   = "#2166AC",
  LumB   = "#92C5DE",
  Her2   = "#D6604D",
  Basal  = "#1A1A1A",
  Normal = "#4DAC26"
)
```

Use this object everywhere. Never hardcode hex values inline in plot calls.

---

## Feature Selection Rules (Stage 04)

- Input: full DEG result per subtype (all genes, not pre-filtered to significant only)
- `glmnet` ElasticNet (`alpha = 0.5`), 10-fold CV, `lambda.1se` for final model
- Target signature: 10–30 genes per subtype
- Stability check: bootstrap 100 iterations, report selection frequency per gene,
  retain genes selected in ≥70% of runs
- Note in comments: selected genes reflect predictive signal, not necessarily causal —
  pathway enrichment (Stage 05) is the interpretability layer

---

## fgsea Rules (Stage 05)

- Gene sets: MSigDB Hallmark (50 sets) + KEGG Legacy via `msigdbr`
- Ranking metric: limma moderated t-statistic (not log2FC alone)
- Min gene set size: 15. Max: 500. `nperm = 1000`.
- Report NES, padj, leading edge genes per pathway per subtype
- Do NOT use the Broad GSEA .jar under any circumstances in v3.0

---

## Visualization Rules (Stage 08)

- Save every publication figure as both `.png` (300 dpi) and `.pdf` (vector)
- Every figure: title, axis labels, caption string assigned in script (not just plot title)
- Use `subtype_colors` palette — consistent across all plots
- `ggplot2` or `ComplexHeatmap` only for final figures — no base R `plot()`

---

## What NOT To Do

- Do NOT refactor working logic without a reason in the commit message
- Do NOT add a new stage without adding it to this file first
- Do NOT combine two pipeline stages in one script to "simplify"
- Do NOT switch from limma-trend to DESeq2 — limma is the v1.0 foundation
- Do NOT run GSEA manually outside R — fgsea handles everything
- Do NOT leave commented-out code blocks or scratch cells in committed scripts
- Do NOT use `<<-` or `attach()`
- Do NOT skip `set.seed(42)` in any script touching random processes
- Do NOT create `.claude/launch.json` or any dev server config — v3.0 has no web servers, Jupyter, or Quarto. When asked to detect dev servers, reply that none exist for this all-R pipeline.

---

## Verification Checklist

```
□ 00: raw data files present. renv::status() clean. All dirs created.
□ 01: brca_se.rds exists. PAM50 counts approx: LumA ~500, LumB ~190, Her2 ~80, Basal ~190.
      No NA in PAM50 colData. Normal samples present.
□ 02: UMAP plots in results/figures/. Sub-subcluster n per subtype printed to console.
□ 03: One DEG .csv per contrast in results/tables/deg/. n sig genes printed per contrast.
□ 04: Signature gene list + stability % saved per subtype.
□ 05: fgsea results per subtype. ≥5 significant pathways (padj<0.05) per subtype.
□ 06: Module-trait correlation heatmap in results/figures/. Top module per subtype named.
□ 07: Cox HR table per signature. At least one signature with significant OS association.
□ 08: All figures as .png + .pdf. subtype_colors consistent across all figures.
```

---

## Slash Commands (.claude/commands/)

- `/simplify` — simplify last code block, no logic changes
- `/doccheck` — verify every function in current file has complete roxygen2
- `/stage-gate` — confirm current stage outputs exist before proceeding
- `/contrast-check` — list active contrasts, confirm they match this file

---

## If Claude Makes a Mistake

End the correction with:
> "Update CLAUDE.md so you don't make that mistake again."

Keep this file under ~175 lines. If it grows past that, consolidate.
