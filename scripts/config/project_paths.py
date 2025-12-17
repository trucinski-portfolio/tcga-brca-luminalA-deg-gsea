# This script sets the repo root and creates data/ results/ docs/ env/
from pathlib import Path
import os

# -----------------------------
# Repo root
# -----------------------------

def get_repo_root():
    # 1) explicit override (HPC-friendly)
    env_root = os.getenv("BIOL616_REPO_ROOT")
    if env_root:
        root = Path(env_root).resolve()
    else:
        # 2) infer from this file's location:
        # .../scripts/config/project_paths.py -> repo root is 3 levels up
        root = Path(__file__).resolve().parents[2]

    required = ["README.md", "scripts", "data"]
    missing = [x for x in required if not (root / x).exists()]
    if missing:
        raise RuntimeError(
            f"Not repo root: {root}\nMissing: {missing}\n"
            "Fix by setting BIOL616_REPO_ROOT to the repo root."
        )
    return root

REPO_ROOT = get_repo_root()

# -----------------------------
# Base directories
# -----------------------------

RAW_INPUTS = REPO_ROOT / "data" / "raw" / "preprocessing_inputs"
PROCESSED = REPO_ROOT / "data" / "processed"
FIGURES = REPO_ROOT / "results" / "figures"
TABLES = REPO_ROOT / "results" / "tables"
DOCS = REPO_ROOT / "docs"
ENV = REPO_ROOT / "env"

PREPROC_OUTPUT_DIR = PROCESSED / "preprocessing_outputs"

# -----------------------------
# Directory creation
# -----------------------------

ALL_DIRS = [RAW_INPUTS, PROCESSED, FIGURES, TABLES, DOCS, ENV]

# DEG / GSEA outputs
DEG_DIR = TABLES / "deg"
GSEA_DIR = TABLES / "gsea"
PLOTS_DIR = FIGURES / "deg_gsea"

ALL_DIRS += [DEG_DIR, GSEA_DIR, PLOTS_DIR]

def ensure_dirs():
    for d in ALL_DIRS:
        d.mkdir(parents=True, exist_ok=True)

# -----------------------------
# File-level paths (defined AFTER ensure_dirs)
# -----------------------------

# Raw inputs (local-only)
EXPR_RAW = RAW_INPUTS / "HiSeqV2_PANCAN.gz"
CLIN_RAW = RAW_INPUTS / "BRCA_clinicalMatrix.tsv"

# Preprocessing outputs (inputs to downstream analysis)
META_PREPROC = PREPROC_OUTPUT_DIR / "metadata_LumA_IDC_Tumor_vs_AllNormals.tsv"
EXPR_PREPROC = PREPROC_OUTPUT_DIR / "expr_LumA_IDC_Tumor_vs_AllNormals.tsv"
