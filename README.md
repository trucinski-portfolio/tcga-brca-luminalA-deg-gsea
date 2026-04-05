# TCGA-BRCA PAM50 Multi-Subtype Transcriptomics Pipeline

Characterizes transcriptional heterogeneity across PAM50 breast cancer subtypes (LumA, LumB,
HER2-enriched, Basal-like) using TCGA-BRCA TOIL RNA-seq — identifying subtype-specific DEGs,
co-expression modules, minimal biomarker signatures, and survival associations.

```
Cohort: 1,109 samples  |  Stages: 8  |  Subtypes: 4  |  Contrasts: 18  |  Figures: 22
```

---

## Overview

This is a fully reproducible, all-R computational pipeline analyzing harmonized TCGA-BRCA
RNA-seq data from the UCSC TOIL compendium. It implements a three-group design — PAM50 tumor
samples, Normal Adjacent Tissue (NAT), and GTEx Healthy breast — to separate tumor-specific
dysregulation from field cancerization effects. All 8 stages run sequentially from a single
`renv`-locked environment.

The dual-reference design enables classification of every differentially expressed gene into:
- **Robust biomarkers** — DE vs both NAT and GTEx (tumor-specific, unaffected by field)
- **Field effect genes** — DE vs GTEx only (present throughout the cancer field)
- **NAT-specific noise** — DE vs NAT only (excluded from signatures)

---

## Data Sources

| Item | Value |
|------|-------|
| Hub | UCSC Xena `toilHub` (`https://toil.xenahubs.net`) |
| Expression | `TcgaTargetGtex_rsem_gene_tpm` — TOIL harmonized, ~19K samples |
| Phenotype | `TcgaTargetGTEX_phenotype.txt` |
| PAM50 calls | `TCGA.BRCA.sampleMap/BRCA_clinicalMatrix` (tcgaHub), column `PAM50Call_RNAseq` |
| Survival | `TCGA.BRCA.sampleMap/BRCA_clinicalMatrix`, columns `OS_Time_nature2012`, `OS_event_nature2012` |
| Units | log2(TPM + 0.001) — pre-transformed by TOIL, never re-logged |
| Genome | GRCh38 / hg38 |

---

## Sample Groups

| Group | n | Notes |
|-------|---|-------|
| **Tumor** | **817** | TCGA primary solid tumors, barcode `01`, PAM50-classified |
| — LumA | 420 | Luminal A |
| — LumB | 192 | Luminal B |
| — Basal | 139 | Basal-like |
| — Her2 | 66 | HER2-enriched |
| **NAT** | **113** | TCGA normal adjacent tissue, barcode `11` |
| **GTEx Healthy** | **179** | Cancer-free donor breast tissue |
| Normal-like | — | PAM50 tumor subtype — **excluded** (not a control group) |

---

## Pipeline Stages

| Stage | Script | Key Output | Description |
|-------|--------|------------|-------------|
| 00 | `R/00_setup.R` | `metadata_breast.rds` | renv check, phenotype download, GTEx field validation |
| 01 | `R/01_preprocessing.R` | `cohort.rds` | Expression download, Ensembl→HGNC mapping, PAM50 join, 3-way labels |
| 02 | `R/02_clustering.R` | `umap_clusters.rds` | HVG→PCA→UMAP, within-subtype k-means subclusters |
| 03 | `R/03_deg.R` | `results/tables/deg/` (18 CSVs) | limma-trend, 6-level model, 18 contrasts |
| 04 | `R/04_feature_selection.R` | `results/tables/signatures/` (4 CSVs) | ElasticNet bootstrap stability selection |
| 05 | `R/05_fgsea.R` | `results/tables/fgsea/` (4 CSVs) | fgsea, Hallmark + KEGG, ranked by t-stat |
| 06 | `R/06_wgcna.R` | `wgcna_modules.rds` | Signed WGCNA, soft-threshold, module-trait correlations |
| 07 | `R/07_survival.R` | `results/tables/survival/` (4 CSVs) | Cox PH, signature score vs Overall Survival |
| 08 | `R/08_visualization.R` | `results/figures/` (44 files) | Publication figures — UMAP, volcanos, heatmaps, dotplot, forest |

---

## Key Results

- **Expression matrix:** 28,344 HGNC gene symbols × 1,109 samples
- **Significant DEGs (FDR < 0.05, |log2FC| ≥ 1):**
  - Subtype vs GTEx: 8,688–9,693 genes per subtype
  - Subtype vs NAT: 5,575–7,091 genes per subtype
  - NAT vs GTEx: 5,689 genes (field cancerization signal)
- **ElasticNet signatures:** LumA 30, LumB 27, Her2 36, Basal 34 genes (freq ≥ 0.80)
- **Pathway enrichment:** 58–93 significant pathways per subtype vs GTEx (padj < 0.05, 236 gene sets)
- **WGCNA:** 5 modules, soft-threshold power = 8 (R² = 0.907), 20/24 significant module-trait associations
- **Survival:** No significant OS association for subtype-discriminating signatures (expected —
  signatures trained for classification, not prognosis). Basal trend: HR = 0.70, p = 0.11.

---

## How to Reproduce

```r
# 1. Clone the repository
# git clone <repo-url> && cd tcga-brca-luminalA-deg-gsea

# 2. Restore the R environment (R 4.5.2, Bioconductor 3.22)
renv::restore()

# 3. Run stages in order (each checks its required inputs before running)
source("renv/activate.R")
Rscript R/00_setup.R
Rscript R/01_preprocessing.R   # ~10 min (large download on first run)
Rscript R/02_clustering.R
Rscript R/03_deg.R
Rscript R/04_feature_selection.R  # ~20 min (bootstrap ElasticNet)
Rscript R/05_fgsea.R
Rscript R/06_wgcna.R           # ~20 min (network construction)
Rscript R/07_survival.R
Rscript R/08_visualization.R
```

**Notes:**
- Stage 01 downloads ~4 GB of expression data on first run; subsequent runs use the cache in `data/raw/`
- `data/raw/` and `data/processed/` are gitignored — they must be regenerated locally
- All intermediate `.rds` files are regenerable; only `results/` is committed

---

## Repository Layout

```
.
├── CLAUDE.md                      # Living pipeline specification and data contracts
├── README.md
├── renv.lock                      # Locked R environment (R 4.5.2)
├── .Rprofile                      # Bootstraps renv on session start
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
│   ├── raw/                       # NOT committed — regenerated by Stage 00–01
│   └── processed/                 # NOT committed — regenerated by Stages 01–07
├── results/
│   ├── figures/                   # Committed — 44 .png + .pdf publication figures
│   └── tables/
│       ├── deg/                   # 18 DEG CSVs
│       ├── signatures/            # 4 ElasticNet signature CSVs
│       ├── fgsea/                 # 4 pathway enrichment CSVs
│       └── survival/              # 4 Cox HR CSVs
└── docs/
    ├── pipeline_report.docx       # Full manuscript-style report
    ├── REFERENCES.md              # Complete citation list
    ├── pipeline_summary.pptx      # Visual summary presentation
    ├── v1-final-presentation.pptx # v1.0 (archived)
    ├── v1.0-finalreport.pdf       # v1.0 (archived)
    └── step-summaries/            # Per-stage technical notes (.docx)
```

---

## References

See [`docs/REFERENCES.md`](docs/REFERENCES.md) for the complete citation list.

Key references:

- **TOIL / UCSC Xena:** Goldman et al. (2020) *Nat Commun* — harmonized RNA-seq compendium
- **PAM50:** Parker et al. (2009) *J Clin Oncol* — intrinsic subtype classifier
- **limma:** Ritchie et al. (2015) *Nucleic Acids Res* — moderated t-statistics for RNA-seq
- **Stability selection:** Meinshausen & Bühlmann (2010) *J R Stat Soc B* — ElasticNet bootstrap
- **fgsea:** Korotkevich et al. (2021) — fast gene set enrichment analysis
- **WGCNA:** Langfelder & Horvath (2008) *BMC Bioinformatics* — weighted co-expression networks
