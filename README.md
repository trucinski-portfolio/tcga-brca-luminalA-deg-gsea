# Luminal A IDC (TCGA-BRCA) — Differential Expression + GSEA Pipeline

This repository implements a **reproducible analysis pipeline** comparing  
**Luminal A invasive ductal carcinoma (IDC) primary tumors vs normal breast tissue**  
using TCGA-BRCA RNA-seq expression and phenotype metadata from **UCSC Xena**.

The pipeline emphasizes:

- **Deterministic paths**
- **Automated input acquisition + validation**
- **Reproducible Python and R environments**
- **Transparent handoff to external GSEA Preranked**
- **Regenerable figures and summary tables**

---

## What this pipeline does

**Stage 0 — Setup & validation**
- Creates the expected directory layout
- Downloads TCGA inputs from UCSC Xena
- Validates required files exist before notebook execution
- Bootstraps R libraries via `renv` 

**Stage 1 — Preprocessing (Python)**
- Reads Xena expression + phenotype tables
- Harmonizes TCGA barcodes
- Constructs the exact cohort contrast
- Writes matched `expr` matrix + `metadata`

**Stage 2 — Differential expression + ranking (R / limma)**
- Runs limma DE
- Exports full DEG results
- Produces a preranked list for GSEA

**Stage 3 — Visualization + GSEA aggregation (R)**
- Volcano plot and heatmap of top DEGs
- Reads the *two* GSEA report TSVs (pos/neg) per library
- Builds a combined enrichment summary + plots

---

## Repo structure

```text
scripts/
├── setup/
│   ├── setup_python.sh          # creates .venv, installs Python deps, registers kernel
│   ├── download_xena_inputs.sh  # downloads TCGA inputs from UCSC Xena
│   ├── prepare_xena_inputs.sh   # validates required input files
│   ├── setup_r.R                # installs R dependencies (renv-aware)
│   └── r_bootstrap.R            # notebook-safe R bootstrap helper
│
├── config/
│   └── project_paths.py         # centralized path conventions (Python)
│
├── notebooks/
│   ├── 01_preprocessing.ipynb
│   ├── 02_deg_analysis_gsea_ranking.ipynb
│   └── 03_deg_gsea_visualization.ipynb
│
data/
├── raw/preprocessing_inputs/    # downloaded Xena inputs (not tracked)
├── processed/preprocessing_outputs/
└── GSEA_output/                 # GSEA pos/neg TSVs (per gene set library)
│
results/
├── tables/                      # DEG + GSEA summary tables
└── figures/                     # volcano, heatmap, GSEA plots
│
env/
├── requirements_python.txt      # Python dependencies
│
renv/                             # renv-managed R library
renv.lock                         # locked R dependencies
README.md
```
---

## Quickstart

### 0) Repository Setup
```bash
git clone https://github.com/yourname/tcga-brca-luminalA-deg-gsea.git
cd tcga-brca-luminalA-deg-gsea

This repo uses:
- **Python** via a local virtual environment at `./.venv` (not committed)
- **R** via `renv` with a committed lockfile (`renv.lock`)

### 1) Python Environment Setup
From the repo root:

```bash
bash scripts/setup/setup_python.sh

This will:
	•	Create .venv/
	•	Install Python dependencies
	•	Register a Jupyter kernel for the project

### 2) R Environment Setup

```bash
R
install.packages("renv")
renv::restore()

This restores the exact R environment defined in renv.lock.

### 3) Download TCGA Xena Inputs

```bash
bash scripts/setup/download_xena_inputs.sh
bash scripts/setup/prepare_xena_inputs.sh

### 4) Run the analysis pipeline

Execute notebooks in order:
	1.	scripts/notebooks/01_preprocessing.ipynb (Python)
	2.	scripts/notebooks/02_differential_expression.ipynb (R)
	3.	scripts/notebooks/03_visualization.ipynb (R)

  ---

## Outputs and directory layout

After running the full pipeline, the following key outputs are produced:

### Preprocessing outputs (Script 01)
`data/processed/preprocessing_outputs/`
`├── expr_LumA_IDC_Tumor_vs_AllNormals.tsv`
`└── metadata_LumA_IDC_Tumor_vs_AllNormals.tsv`

These files define the **final cohort and expression matrix** used in all downstream analyses.

---

### Differential expression results (Script 02)
`results/tables/deg/`
`├── DEG_Results_LumA_IDC_Tumor_vs_AllNormals_limma.csv`
`└── DEG_Results_LumA_IDC_Tumor_vs_AllNormals_limma_annotated.csv`

`results/figures/deg/`
`├── Top50_DEGs_Heatmap.pdf`
`├── Top50_DEGs_Heatmap.png`
`└── Volcano_LumA_IDC_Tumor_vs_AllNormals.png`

Differential expression is performed using **limma-trend**, and the ranked list for GSEA is derived from **limma t-statistics**.

---

### Pathway enrichment summaries and figures (Script 03)
`results/tables/gsea/`
`└── GSEA_ALL_summary_LumA_IDC_Tumor_vs_AllNormals.csv`

This table merges positive and negative enrichment results across all gene set libraries and reports:
	•	Normalized Enrichment Score (NES)
	•	Enrichment direction (Up_in_Tumor / Down_in_Tumor)
	•	Gene set database source

Results are sorted by absolute enrichment strength.

`results/figures/gsea/`
`├── GSEA_Waterfall_Hallmark.png`
`├── GSEA_Waterfall_KEGG_Legacy.png`
`└── GSEA_Waterfall_KEGG_Medicus.png`

Each waterfall plot visualizes the top enriched pathways per library, ranked by NES, and colored by direction of enrichment in tumor tissue.

These outputs summarize differential gene expression and pathway enrichment across multiple gene set libraries
(Hallmark, KEGG Legacy, KEGG Medicus) and provide a comparative view of proliferation, DNA repair, proteostasis, and signaling programs active in Luminal A IDC.

---

## GSEA Preranked (external step)

GSEA Preranked is run **outside** this repository using the Broad Institute’s
`GSEA_4.x` software.

Required inputs:
- Ranked list produced by Script 02 (`.rnk` file)

Expected outputs (per gene set library):
- `gsea_report_for_na_pos_*.tsv`
- `gsea_report_for_na_neg_*.tsv`

These files must be placed in:
- `data/GSEA_output/Hallmark-All/`
- `data/GSEA_output/KEGG-Legacy/`
- `data/GSEA_output/KEGG-Medicus/`

Script 03 automatically detects **the most recent positive and negative report per
library based on file modification time** and merges them into a single enrichment
summary.

---

## Reproducibility notes

- Python dependencies are defined in `env/requirements_python.txt` and installed
  into a local `.venv`
- R dependencies are fully locked and restorable via `renv.lock`
- The `.venv/` directory is intentionally **not committed**
- Large raw TCGA matrices are not tracked in Git

This design allows the analysis to be **fully reproducible** while keeping the
repository lightweight.

---

### Directory summary

- `scripts/`   — analysis notebooks, setup scripts, and path/config helpers  
- `data/`      — raw inputs, processed data, and external GSEA outputs  
- `results/`   — tables and figures generated by the pipeline  
- `env/`       — environment specifications (Python and R)  
- `docs/`      — reports, slides, and supplementary materials  

---

## Authorship

**Thomas Rucinski**
Biomedical Engineering Masters Student at University of Nevada Las Vegas