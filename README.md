## Luminal A Invasive Ductal Carcinoma (TCGA-BRCA): Differential Expression and Pathway Enrichment Pipeline

This repository contains a reproducible analysis pipeline, emphasizing barcode harmonization, strict cohort definition, reproducible directory structure, and statistically principled ranking for GSEA. The pipeline focuses on differential gene expression and pathway enrichment (GSEA Preranked) comparing Luminal A invasive ductal carcinoma (IDC) primary tumors vs normal breast tissue using TCGA-BRCA RNA-seq and clinical phenotype data.

The workflow is structured as three main stages:
**1. Preprocessing (Python)**: load TCGA expression + phenotype tables, harmonize TCGA barcodes, filter to the study contrast, and write a clean metadata table + matched expression matrix.
**2. Differential expression (R / limma-trend)**: run limma, export full DEG results, and generate a t-statistics–based preranked list for GSEA.
**3. Visualization (R)**: volcano plot, heatmap of top DEGs, and combined GSEA “waterfall” plots across multiple gene set libraries.

**Manual step**: GSEA Preranked is run externally (UCSD/Broad Institute's GSEA_4.4.0 Software). The pipeline expects the resulting 'gsea_report_for_na_pos*.tsv' and 'gsea_report_for_na_neg*.tsv' files to be placed into 'data/GSEA_output/<GENESET_LIBRARY>/' before running the visualization notebook.

---

## Repository structure

- scripts/
  - 01_preprocessing.ipynb - preprocessing + cohort construction
  - 02_differential_expression.ipynb - limma DE + t-stat ranked list
  - 03_visualization.ipynb - DEG-annotation + plots + GSEA summaries
  - config/ - project paths and directory setup helpers
  - setup/ - helper scripts (e.g., dataset download instructions)
- data/ - raw inputs + external tool outputs (**see data/README.md**)
- results/
  - tables/ - DEG results and GSEA summary files 
  - figures/ - volcano plots, heatmaps, GSEA waterfall plots
- env/ - environment/requirements files 
- docs/ - optional project files (including 6 Page Report Detailing Methods and Findings - **The report used slightly outdated methods, and the figures are inaccurate for the current methods used**)

---

## Quickstart

### 0) Clone and move to repo root
Run everything from the repository root (the folder containing README.md, scripts/, and data/).

### 1) Preprocessing (Python)
Run scripts/01_preprocessing.ipynb to produce:

- data/processed/preprocessing_outputs/metadata_LumA_IDC_Tumor_vs_AllNormals.tsv
- data/processed/preprocessing_outputs/expr_LumA_IDC_Tumor_vs_AllNormals.tsv

### 2) Differential expression (R / limma)
Run scripts/02_differential_expression.ipynb to produce:

- DEG table(s) in results/tables/deg/
- GSEA preranked input file (t-stat based) in:
  - data/processed/GSEA_input/LumA_IDC_Tumor_vs_AllNormals.rnk

### 3) Run GSEA Preranked (manual, external tool)
Run UCSD & Broad Institute's GSEA_4.4.0 Software using **GSEA Preranked** using the .rnk file above.

For each gene set library (Hallmark, KEGG Legacy, KEGG Medicus), copy the following outputs into:

- data/GSEA_output/Hallmark-All/
- data/GSEA_output/KEGG-Legacy/
- data/GSEA_output/KEGG-Medicus/

Required files per library:
- gsea_report_for_na_pos_*.tsv
- gsea_report_for_na_neg_*.tsv

(Details and expectations are documented in data/README.md)

### 4) Visualization + summary tables (R)
Run scripts/03_visualization.ipynb to generate:
- results/tables/gsea/GSEA_ALL_summary_*.csv
- results/figures/gsea/GSEA_Waterfall_*.png
- heatmaps and DEG visualizations in results/figures/

---

## Notes on interpretation

- Differential expression is estimated via **limma-trend**.
- The GSEA ranked list is derived from **limma t-statistics** (directionality preserved).
  - Positive enrichment (POS) corresponds to pathways enriched toward the top of the ranked list.
  - Negative enrichment (NEG) corresponds to pathways enriched toward the bottom of the ranked list.

---

## Data and licensing

Raw TCGA matrices and GSEA output artifacts are not tracked in Git (size). See data/README.md for exact expectations and directory layout.

---

## Author

Thomas Rucinski
M.S. Candidate | Biomedical Engineering | University of Nevada, Las Vegas 
Email: thomasrucinski13@gmail.com
