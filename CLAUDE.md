# CLAUDE.md — TCGA-BRCA Multi-Subtype Pipeline (v3.2)

> Living document. After every correction: "Update CLAUDE.md so you don't make that mistake again."
> v3.2 fixes dataset names, corrects the 3-way sample filter logic, clarifies
> Normal-like vs NAT vs GTEx, and adds live plot display on top of file saves.

---

## Project in One Sentence

Characterize transcriptional heterogeneity across PAM50 breast cancer subtypes (LumA, LumB,
HER2-enriched, Basal-like) using TCGA-BRCA TOIL RNA-seq — identifying subtype-specific DEGs
against NAT and GTEx Healthy controls, co-expression modules, minimal biomarker signatures,
and survival associations — to produce a publishable-quality, thesis-defensible, all-R pipeline.

---

## Current Stage

**Active:** Stage 08 complete — pipeline finished.
**Last confirmed output:** 44 publication figures across 7 figure sets. All stages 00–08 confirmed complete.

> When starting a session, read this field first. Run `/stage-gate` before writing any code.

---

## Language Contract

**All-R. No Python. No shell scripts. No manual downloads.**

| Task | Package |
|------|---------|
| Data retrieval | `UCSCXenaTools` only |
| Preprocessing | `dplyr`, `tidyr`, `SummarizedExperiment` |
| ID mapping | `AnnotationDbi` + `org.Hs.eg.db` (local, no network) — never `biomaRt` |
| DEG | `limma` (trend=TRUE) — never DESeq2 |
| GSEA | `fgsea` — never the Broad GSEA .jar |
| Feature selection | `glmnet` (ElasticNet) |
| Co-expression | `WGCNA` |
| Survival | `survival` + `survminer` |
| Visualization | `ggplot2`, `ComplexHeatmap` — never base R `plot()` |
| Environment | `renv` (lockfile committed) |

---

## UCSCXenaTools — Authoritative API Pattern

Claude Code frequently hallucinates UCSCXenaTools function names. Use ONLY this pipeline:
```r
library(UCSCXenaTools)
library(dplyr)

# --- STEP 1: Download phenotype metadata (ALWAYS before expression) ---
pheno_dl <- XenaGenerate(subset = XenaHostNames == "toilHub") |>
  XenaFilter(filterDatasets = "TcgaTargetGTEX_phenotype") |>
  XenaQuery() |>
  XenaDownload()

pheno <- XenaPrepare(pheno_dl)

# --- STEP 2: Inspect phenotype columns before hardcoding any field names ---
# Confirmed at Stage 00: disease = `primary disease or tissue`,
#                        study   = `_study`,  site = `_primary_site`
colnames(pheno)

# --- STEP 3: Build 3-way sample list (TCGA tumor/NAT + GTEx breast) ---
# This is an OR filter — never apply TCGA arm alone or GTEx drops silently
breast_samples <- pheno |>
  filter(
    `primary disease or tissue` == "Breast Invasive Carcinoma" |  # TCGA tumor (01) + NAT (11)
    (`_study` == "GTEX" & `_primary_site` == "Breast")           # GTEx healthy breast
  ) |>
  pull(sample)

# --- STEP 4: Fetch expression for filtered samples only ---
# Use the COMBINED matrix (TCGA + GTEx). tcga_RSEM_gene_tpm is TCGA-only — do not use it.
expr_dl <- XenaGenerate(subset = XenaHostNames == "toilHub") |>
  XenaFilter(filterDatasets = "TcgaTargetGtex_rsem_gene_tpm") |>
  XenaQuery() |>
  XenaDownload()

expr_raw <- XenaPrepare(expr_dl)

# --- STEP 5: Subset to breast samples ---
expr <- expr_raw[, colnames(expr_raw) %in% breast_samples]
```

**Never use:** `XenaHub()` directly, `getTCGA()`, `downloadTCGA()`, `TCGAbiolinks`,
`GEOquery`, or any unfiltered full-matrix download before phenotype subsetting.

---

## Data Contract (v3.2 — TOIL)

| Item | Value |
|------|-------|
| Hub | `https://toil.xenahubs.net` (alias: `"toilHub"` in `XenaGenerate`) |
| Expression dataset | `TcgaTargetGtex_rsem_gene_tpm` — combined TCGA + GTEx, ~19K samples |
| Phenotype dataset | `TcgaTargetGTEX_phenotype.txt` |
| Units | log2(TPM + 0.001) — **do not re-log** |
| Genome | GRCh38 (hg38) |
| Row IDs | Ensembl gene IDs → map to HGNC symbols in Stage 01 |

**Critical dataset naming notes:**
- `TcgaTargetGtex_rsem_gene_tpm` — correct combined matrix (TCGA + GTEx in same pipeline)
- `tcga_RSEM_gene_tpm` — TCGA-only, NO GTEx samples. Do not use for this project.
- `TcgaTargetGTEX_phenotype.txt` — correct phenotype file (toilHub)
- `tcga_target_gtex_samples` — does not exist on the hub
- `TCGA.BRCA.sampleMap/BRCA_clinicalMatrix` — PAM50 source (tcgaHub). The TOIL phenotype
  does NOT embed PAM50 subtypes in `detailed_category` for BRCA samples. PAM50 must be
  fetched separately from `tcgaHub` and joined on the 15-char TCGA barcode.
  Column: `PAM50Call_RNAseq`. Values: `"LumA"`, `"LumB"`, `"Her2"`, `"Basal"`, `"Normal"`.
  Overall Survival columns (confirmed): `OS_Time_nature2012` (days), `OS_event_nature2012` (0/1).
  NOT `OS.time` / `OS` — those column names do not exist in this dataset.
  **Never try to extract PAM50 from `detailed_category` — it only says "Breast Invasive Carcinoma".**

---

## 3-Way Sample Groups — Definitions and Filter Logic

There are THREE distinct concepts involving the word "normal." Never conflate them.

| Term | What it is | Barcode / ID | Included? |
|------|-----------|--------------|-----------|
| **Tumor** | TCGA Primary Solid Tumor | TCGA barcode pos 14–15 = `01` | ✅ Yes |
| **NAT** | TCGA Normal Adjacent Tissue — non-tumor biopsy from same cancer patient | TCGA barcode pos 14–15 = `11` | ✅ Yes |
| **GTEx Healthy** | Truly healthy breast from cancer-free donors | Non-TCGA sample ID, `study == "GTEX"` | ✅ Yes |
| **Normal-like** | A PAM50 **tumor subtype** — tumor sample (barcode `01`) whose expression clusters near normal | barcode `01`, PAM50 = "Normal" | ❌ Excluded |

**Normal-like is a tumor subtype label, not a tissue type.** Excluding it does not remove
any NAT or GTEx samples. NAT (barcode `11`) and GTEx are control groups — they are
never filtered by PAM50 label.

### Sample filter must be an OR across two arms:
```r
# CORRECT — confirmed column names from Stage 00 inspection
filter(
  `primary disease or tissue` == "Breast Invasive Carcinoma" |   # TCGA tumor 01 + NAT 11
  (`_study` == "GTEX" & `_primary_site` == "Breast")            # GTEx healthy breast
)

# WRONG — silently drops all GTEx samples
filter(`primary disease or tissue` == "breast invasive carcinoma")
```

**Verify GTEx field name at Stage 00** by inspecting `colnames(pheno)` — the column may be
`_primary_site` or `X_primary_site` depending on R's handling of the leading underscore.
Record the confirmed name here before Stage 01.

**Confirmed field names (recorded at Stage 00):**

All 7 columns in `TcgaTargetGTEX_phenotype.txt`: `sample`, `detailed_category`,
`` `primary disease or tissue` ``, `` `_primary_site` ``, `` `_sample_type` ``,
`` `_gender` ``, `` `_study` ``

- Disease field: `` `primary disease or tissue` `` — value is **`"Breast Invasive Carcinoma"`** (title case, n=1,212)
- Study field: `` `_study` `` (not `study`)
- GTEx site field: `` `_primary_site` `` (no `X_` prefix — confirmed exact)

Do NOT use lowercase `"breast invasive carcinoma"` — the field is title case and the filter will silently return 0 TCGA rows.

### Expected sample counts after Stage 01:
- Tumor (PAM50: LumA, LumB, Her2, Basal): ~1,100
- NAT: ~113
- GTEx Healthy Breast: ~200+
- Normal-like tumors: excluded before analysis

---

## PAM50 Subtype Filtering

Applied to Tumor group only. NAT and GTEx groups are never filtered by PAM50.
```r
# After group assignment, exclude Normal-like from the Tumor arm only
meta <- meta |>
  filter(
    group != "Tumor" |                         # keep all NAT and GTEx
    (group == "Tumor" & PAM50 %in% c("LumA", "LumB", "Her2", "Basal"))
  )
```

Subtypes in scope: `LumA`, `LumB`, `Her2`, `Basal`. Normal-like excluded.

---

## 3-Way Contrast Design — 18 Contrasts Total

| Layer | Contrasts | n |
|-------|-----------|---|
| Subtype vs NAT | LumA, LumB, Her2, Basal — each vs NAT | 4 |
| Subtype vs GTEx | LumA, LumB, Her2, Basal — each vs GTEx Healthy | 4 |
| Subtype vs Subtype | LumA/Basal, LumA/LumB, LumA/Her2, Her2/Basal, LumB/Basal | 5 |
| NAT vs GTEx | Field cancerization characterization | 1 |
| Within-subtype | UMAP sub-subclusters (Stage 02 populates) | 4 |

All three-group contrasts: **unpaired design.** GTEx is never paired with TCGA.
GTEx Healthy = primary reference (absolute dysregulation).
NAT = secondary reference (clinical margin biology).

**Interpretive categories (post-Stage-03 intersection logic):**
- **Robust biomarkers** — DE in Subtype vs NAT AND Subtype vs GTEx
- **Field effect genes** — DE in Subtype vs GTEx, NOT Subtype vs NAT
- **NAT-specific noise** — DE in Subtype vs NAT, NOT Subtype vs GTEx (exclude from signatures)

Do not add contrasts without updating this table first.

---

## Visualization Rule — Display AND Save

Every plot must do both: render to screen during the run AND save to disk.
```r
# ggplot2 — always assign, then print, then save
p <- ggplot(data, aes(...)) + geom_*() + ...

print(p)  # renders to screen when script is run interactively or via source()

ggsave(
  filename = here::here("results", "figures", "filename.png"),
  plot = p, width = 10, height = 8, dpi = 300
)
ggsave(
  filename = here::here("results", "figures", "filename.pdf"),
  plot = p, width = 10, height = 8
)
```
```r
# ComplexHeatmap — draw to screen, then save via pdf/png device
ht <- Heatmap(matrix, ...)

draw(ht)  # renders to screen

png(here::here("results", "figures", "filename.png"), width = 2400, height = 2000, res = 300)
draw(ht)
dev.off()

pdf(here::here("results", "figures", "filename.pdf"), width = 10, height = 8)
draw(ht)
dev.off()
```

**Rules:**
- `print()` or `draw()` before every `ggsave()` / device block — never save-only
- Save as both `.png` (300 dpi) AND `.pdf` (vector) every time
- Every figure has: title, axis labels, caption string variable
- Use `subtype_colors` palette — consistent across all plots
- `ggplot2` or `ComplexHeatmap` only — never base R `plot()`

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

Never hardcode hex values inline. Use this object in every plot call.

---

## Inter-Stage Data Contracts
```
Stage 00 output → data/processed/metadata_breast.rds
  $pheno: data.frame
  Columns confirmed at runtime: sample, primary_disease, study, PAM50,
  sample_type, _primary_site (or X_primary_site — verify and record above)

Stage 01 input  → data/processed/metadata_breast.rds
Stage 01 output → data/processed/cohort.rds
  $expr:      matrix — rows = HGNC symbols, cols = sample IDs
  $meta:      data.frame — cols: sample_id, group (Tumor/NAT/GTEx), PAM50
  $contrasts: character vector — length == 18, names match contrast table

Stage 02 input  → data/processed/cohort.rds
Stage 02 output → data/processed/umap_clusters.rds
  $umap:          data.frame — cols: sample_id, UMAP1, UMAP2, subtype, subcluster
  $n_per_subtype: named integer vector

Stage 03 input  → data/processed/cohort.rds + data/processed/umap_clusters.rds
Stage 03 output → results/tables/deg/*.csv (18 files)
  Each file cols: gene, logFC, AveExpr, t, P.Value, adj.P.Val, B

Stage 04 input  → results/tables/deg/*.csv (all genes, not pre-filtered)
Stage 04 output → results/tables/signatures/*.csv (one per subtype)
  Each file cols: gene, selection_frequency, mean_coef

Stage 05 input  → results/tables/deg/*.csv
Stage 05 output → results/tables/fgsea/*.csv (one per subtype)
  Each file cols: pathway, NES, padj, leadingEdge

Stage 06 input  → data/processed/cohort.rds
Stage 06 output → data/processed/wgcna_modules.rds

Stage 07 input  → results/tables/signatures/*.csv + data/processed/cohort.rds ($meta)
Stage 07 output → results/tables/survival/*.csv (Cox HR tables)

Stage 08 input  → all of the above
Stage 08 output → results/figures/*.png + results/figures/*.pdf
```

Never run stage N before stage N-1 output exists.
Each script checks for its required input at the top and stops with a clear error if missing.

---

## Pipeline Stages
```
00_setup.R             →  renv verified, hub reachable, phenotype inspected,
                           GTEx field name confirmed, metadata_breast.rds cached
01_preprocessing.R     →  cohort.rds (expr + meta with 3-way labels + 18 contrasts)
02_clustering.R        →  umap_clusters.rds (within-subtype UMAP + sub-subclusters)
03_deg.R               →  results/tables/deg/ (18 .csv files)
04_feature_selection.R →  results/tables/signatures/ (ElasticNet gene sets per subtype)
05_fgsea.R             →  results/tables/fgsea/ (NES tables per subtype)
06_wgcna.R             →  data/processed/wgcna_modules.rds
07_survival.R          →  results/tables/survival/ (Cox HR tables)
08_visualization.R     →  results/figures/ (all publication figures)
```

---

## Repo Layout
```
project-root/
├── CLAUDE.md
├── README.md
├── renv.lock
├── .Rprofile                         ← bootstraps renv on session start
├── .gitignore                        ← data/raw/ and data/processed/ excluded
├── R/
│   ├── 00_setup.R
│   ├── 01_preprocessing.R
│   ├── 02_clustering.R
│   ├── 03_deg.R
│   ├── 04_feature_selection.R
│   ├── 05_fgsea.R
│   ├── 06_wgcna.R
│   ├── 07_survival.R
│   └── 08_visualization.R
├── data/
│   ├── raw/                          ← NOT committed
│   └── processed/                    ← NOT committed
├── results/
│   ├── figures/                      ← committed (.png + .pdf)
│   └── tables/                       ← committed (.csv)
└── .claude/
    └── commands/
```

---

## Coding Conventions

- `here::here()` for all paths — no hardcoded strings
- `set.seed(42)` at top of any script using random processes
- Scripts are **idempotent** — re-running produces same output, never appends
- Intermediate objects: `.rds` in `data/processed/`; final tables: `.csv` in `results/tables/`
- One script = one stage — never combine
- `message()` at start and end of every major function block
- All `library()` calls at top of script — never inside functions
- No `attach()`. No `<<-`.

---

## Documentation — Roxygen2 on Every Function
```r
#' Run a limma-trend contrast for one subtype against one reference group
#'
#' Fits a three-level factor limma-trend model and extracts the contrast for
#' the specified subtype vs reference. Expression data must already be in
#' log2(TPM+0.001) space — do not re-log.
#'
#' @param se SummarizedExperiment. colData must have columns: PAM50, group.
#' @param subtype Character. One of "LumA", "LumB", "Her2", "Basal".
#' @param ref Character. One of "NAT", "GTEx". Reference group for contrast.
#' @param fdr_cutoff Numeric. BH-adjusted FDR threshold. Default 0.05.
#'
#' @return data.frame with cols: gene, logFC, AveExpr, t, P.Value, adj.P.Val, B.
#'
#' @examples
#' deg <- run_limma_contrast(se = brca_se, subtype = "LumA", ref = "GTEx")
run_limma_contrast <- function(se, subtype, ref, fdr_cutoff = 0.05) { ... }
```

No stubs. No TODO comments left in committed code.

---

## Stage-Specific Rules

**Stage 00:**
- Print `colnames(pheno)` and record exact GTEx tissue field name in this file
- Print sample counts: Tumor / NAT / GTEx Healthy split — confirm expected ranges
- Hub reachability check must pass before any download

**Stage 03 (DEG):**
- `limma-trend` only (`trend = TRUE` in `eBayes()`). Never DESeq2, never edgeR.
- `model.matrix` uses a three-level factor. Run two separate model fits:
  GTEx as reference level, then NAT as reference level.
- Output one .csv per contrast. 18 files total.
- Filename convention: `{SubtypeA}_vs_{SubtypeB}.csv`

**Stage 04 (Feature Selection):**
- Input: all genes from DEG result — do NOT pre-filter to significant only
- `glmnet` ElasticNet: `alpha = 0.5`, 10-fold CV, `lambda.1se`
- Target: 10–30 genes per subtype
- Bootstrap 100 iterations; retain genes selected in ≥70% of runs

**Stage 05 (fgsea):**
- Gene sets: MSigDB Hallmark + KEGG Legacy via `msigdbr`
- Ranking metric: limma moderated t-statistic (not log2FC)
- `minSize = 15`, `maxSize = 500`, `nperm = 1000`

---

## What NOT To Do

- Do NOT use `tcga_RSEM_gene_tpm` — it is TCGA-only, GTEx samples are absent
- Do NOT use `tcga_target_gtex_samples` as a dataset name — it does not exist on the hub
- Do NOT use `XenaHub()` directly — always use `XenaGenerate() |> XenaFilter() |> XenaQuery() |> XenaDownload()`
- Do NOT filter samples with only `primary_disease == "breast invasive carcinoma"` — GTEx drops silently
- Do NOT confuse Normal-like (PAM50 tumor subtype, barcode `01`) with NAT (barcode `11`) or GTEx Healthy
- Do NOT apply PAM50 subtype filter to NAT or GTEx samples — PAM50 filtering is for the Tumor group only
- Do NOT re-log expression data — TOIL is already log2(TPM+0.001)
- Do NOT download expression data before phenotype filtering — metadata first, always
- Do NOT use HiSeqV2 or DCC RSEM counts — TOIL TPM only
- Do NOT use TCGAbiolinks or GEOquery
- Do NOT switch from limma-trend to DESeq2
- Do NOT use the Broad GSEA .jar — fgsea only
- Do NOT save a plot without also printing it to screen (`print()` or `draw()` before every save)
- Do NOT save only `.png` — every figure needs both `.png` (300 dpi) and `.pdf` (vector)
- Do NOT combine two pipeline stages in one script
- Do NOT add contrasts without updating the contrast table in this file
- Do NOT refactor working logic without a reason in the commit message
- Do NOT leave commented-out code or TODO stubs in committed scripts
- Do NOT use `attach()` or `<<-`
- Do NOT skip `set.seed(42)` in scripts touching random processes
- Do NOT create `.claude/launch.json` — this is a pure R analysis pipeline, no servers

---

## Slash Commands

| Command | Purpose |
|---------|---------|
| `/stage-gate` | Confirm prior stage outputs exist before writing new code |
| `/contrast-check` | List active contrasts, verify count == 18, confirm match to table above |
| `/doccheck` | Verify every function in current file has complete roxygen2 |
| `/simplify` | Simplify last code block — no logic changes |

---

## Verification Checklist
```
□ 00: renv::status() clean. Hub reachable. GTEx field name confirmed + recorded above.
      metadata_breast.rds cached. Sample counts printed: Tumor ~1100, NAT ~113, GTEx ~200+.

☑ 01: cohort.rds exists.
      $expr: 28,344 HGNC symbols × 1,109 sample IDs.
      $meta: Tumor 817 (LumA 420, LumB 192, Her2 66, Basal 139), NAT 113, GTEx 179.
      $contrasts: length == 18. Normal-like absent from Tumor arm.
      NAT and GTEx arms present and not filtered by PAM50.
      PAM50 sourced from TCGA.BRCA.sampleMap/BRCA_clinicalMatrix (tcgaHub).

□ 02: umap_clusters.rds exists. UMAP plots printed to screen AND saved to
      results/figures/ as .png + .pdf. Sub-subcluster n per subtype printed.

□ 03: 18 .csv files in results/tables/deg/. Sig gene counts printed per contrast.

☑ 04: Signature .csv per subtype. selection_frequency column present.
      FREQ_CUTOFF = 0.80 for all subtypes.
      LumA 30 genes, LumB 27, Her2 36, Basal 34 (Her2/Basal slightly broader due to small n).

☑ 05: fgsea .csv per subtype. ≥5 pathways at padj < 0.05 per subtype.
      Both vs_GTEx and vs_NAT contrasts per subtype. 58–93 sig pathways per contrast.
      236 gene sets: 50 Hallmark + 186 KEGG Legacy.

☑ 06: wgcna_modules.rds. 5 modules (excl. grey), power=8 (R²=0.907).
      20/24 significant module-trait associations (p<0.05).
      Module-trait heatmap printed to screen AND saved as .png + .pdf.

☑ 07: Cox HR .csv per subtype. 4 files written. No subtype reached p<0.05 (expected —
      signatures trained for subtype identity, not survival). Basal trend HR=0.70, p=0.11.
      KM figures saved as .png + .pdf for all 4 subtypes.

☑ 08: All 7 figures (44 files: .png + .pdf) saved to results/figures/.
      Fig1 UMAP overview, Fig2 DEG counts, Fig3×4 volcanos, Fig4 DEG heatmap,
      Fig5 signature heatmap, Fig6 fgsea dotplot, Fig7 Cox forest plot.
      subtype_colors consistent across all figures.
```

---

## If Claude Makes a Mistake

Correct it, then end with:
> "Update CLAUDE.md so you don't make that mistake again."