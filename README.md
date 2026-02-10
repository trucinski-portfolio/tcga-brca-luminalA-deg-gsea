# Luminal A IDC (TCGA-BRCA) — Differential Expression + GSEA Pipeline

This repository implements a **reproducible analysis pipeline** comparing
**Luminal A invasive ductal carcinoma (IDC) primary tumors vs normal breast tissue**
using TCGA-BRCA RNA-seq expression and phenotype metadata from **UCSC Xena**.

The pipeline emphasizes:

- **Deterministic paths**
- **Automated input acquisition + validation**
- **Reproducible Python and R environments**
- **In-notebook fGSEA pathway analysis** (no external GSEA software required)
- **QC/batch-effect checks, clinical heterogeneity analysis**
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

**Stage 2 — Differential expression + fGSEA (R / limma)**
- Runs limma DE with batch-aware design
- Exports full DEG results (annotated + raw)
- Runs fGSEA (Hallmark, KEGG) in-notebook using limma t-statistics
- Produces volcano plot, heatmaps (ordered by UP/DOWN regulation), and fGSEA dotplots
- Exports ranked gene list (`.rnk`)

**Stage 3 — Clinical visualization + tumor heterogeneity (R)**
- Unsupervised clustering of tumors by DEG expression
- Silhouette-based cluster optimization
- Age and stage associations across clusters
- Tumor heterogeneity scoring (concordance of up/down-regulated genes)

---

## Repo structure

```text
scripts/
├── setup/
│   ├── setup_python.sh          # creates .venv, installs Python deps, registers kernel
│   ├── download_xena_inputs.sh  # downloads TCGA inputs from UCSC Xena
│   ├── prepare_xena_inputs.sh   # validates required input files
│   ├── r_bootstrap.R            # notebook-safe R bootstrap helper
│   ├── check_r.R                # R environment check
│   └── render_reports.sh        # renders notebooks via Quarto
│
├── config/
│   └── project_paths.py         # centralized path conventions (Python)
│
├── analysis/
│   └── batch_effect_check.R     # standalone batch-effect PCA analysis
│
├── notebooks/
│   ├── 01_preprocessing.ipynb
│   ├── 02_deg_fgsea_analysis.ipynb
│   └── 03_visualization_clinical.ipynb

data/
├── raw/preprocessing_inputs/    # downloaded Xena inputs (not tracked)
├── processed/preprocessing_outputs/
└── GSEA_output/                 # external GSEA pos/neg TSVs (legacy, optional)

results/
├── figures/
│   ├── deg/                     # volcano, heatmaps (ordered by UP/DOWN)
│   ├── gsea/                    # fGSEA dotplots, GSEA waterfall plots
│   ├── qc/                      # PCA plots by group and TSS
│   └── clinical/                # clustering, heterogeneity, age/stage plots
├── tables/
│   ├── deg/                     # DEG results, cluster assignments, heterogeneity scores
│   └── gsea/                    # fGSEA results, GSEA summaries, ranked list

env/
├── requirements_python.txt      # Python dependencies
└── requirements_python_runtime.txt

docs/
├── final-report.pdf
└── slides.pptx

renv/                            # renv-managed R library
renv.lock                        # locked R dependencies
_quarto.yml                      # Quarto rendering config
README.md
```

---

## Quickstart

### 0) Clone the repository
```bash
git clone https://github.com/trucinski-portfolio/tcga-brca-luminalA-deg-gsea.git
cd tcga-brca-luminalA-deg-gsea
```

### 1) Python environment setup
```bash
bash scripts/setup/setup_python.sh
```
This will create `.venv/`, install Python dependencies, and register a Jupyter kernel.

### 2) R environment setup
```r
install.packages("renv")
renv::restore()
```
This restores the exact R environment defined in `renv.lock`.

### 3) Download TCGA Xena inputs
```bash
bash scripts/setup/download_xena_inputs.sh
bash scripts/setup/prepare_xena_inputs.sh
```

### 4) Run the analysis pipeline

Execute notebooks in order:
1. `scripts/notebooks/01_preprocessing.ipynb` (Python)
2. `scripts/notebooks/02_deg_fgsea_analysis.ipynb` (R)
3. `scripts/notebooks/03_visualization_clinical.ipynb` (R)

Or render all as HTML reports via Quarto:
```bash
bash scripts/setup/render_reports.sh
```

---

## Outputs

### Preprocessing (Notebook 01)
```text
data/processed/preprocessing_outputs/
├── expr_LumA_IDC_Tumor_vs_AllNormals.tsv
└── metadata_LumA_IDC_Tumor_vs_AllNormals.tsv
```

### DEG + fGSEA analysis (Notebook 02)
```text
results/figures/deg/
├── Volcano_LumA_IDC_Tumor_vs_AllNormals.png
├── Top50_DEGs_Heatmap.png
├── Top50_DEGs_Heatmap_ordered.png
├── Top50_DEGs_Heatmap_ordered_by_UP.png
└── Top50_DEGs_Heatmap_ordered_by_DOWN.png

results/figures/gsea/
├── fgsea_Hallmark_dotplot.png
└── fgsea_KEGG_dotplot.png

results/tables/deg/
├── DEG_Results_LumA_IDC_Tumor_vs_AllNormals_limma.csv
└── DEG_Results_LumA_IDC_Tumor_vs_AllNormals_limma_annotated.csv

results/tables/gsea/
├── fgsea_results_all_databases.csv
└── LumA_IDC_Tumor_vs_AllNormals.rnk
```

### QC / batch-effect analysis
```text
results/figures/qc/
├── PCA_by_Group.png
├── PCA_by_TSS.png
├── PCA_by_TSS_faceted.png
└── PCA_PC3_PC4_by_TSS.png
```

### Clinical visualization + tumor heterogeneity (Notebook 03)
```text
results/figures/clinical/
├── silhouette_optimization.png
├── silhouette_optimization_tumors.png
├── 2D_space_by_stage.png
├── age_by_cluster.png
├── age_by_tumor_cluster.png
├── tumor_heterogeneity_discordance.png
├── tumor_heterogeneity_down_genes_dist.png
└── tumor_heterogeneity_down_vs_up.png

results/tables/deg/
├── sample_cluster_assignments.csv
└── tumor_heterogeneity_scores.csv
```

---

## Reproducibility notes

- Python dependencies are defined in `env/requirements_python.txt` and installed into a local `.venv`
- R dependencies are fully locked and restorable via `renv.lock`
- The `.venv/` directory is intentionally **not committed**
- Large raw TCGA matrices are not tracked in Git
- Quarto config (`_quarto.yml`) enables reproducible HTML report rendering

---

## Authorship

**Thomas Rucinski**
Biomedical Engineering Masters Student at University of Nevada Las Vegas
