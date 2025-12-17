#!/usr/bin/env bash

TARGET_DIR="data/raw/preprocessing_inputs"
mkdir -p "$TARGET_DIR"

echo "=== UCSC Xena TCGA-BRCA Data Download Instructions ==="
echo
echo "1. Visit:"
echo "   https://xenabrowser.net/datapages/"
echo
echo "2. Download the following datasets: (see data/README.md for direct links)"
echo "  - Gene expression matrix:"
echo "    TCGA Breast Cancer (BRCA) - gene expression RNAseq - IlluminaHiSeq pancan normalized"
echo "  - Clinical metadata:"
echo "    TCGA Breast Cancer (BRCA) - dataset: phenotype - Phenotypes"
echo
echo "3. Place the downloaded files into:"
echo "   $TARGET_DIR"
echo
echo "4. Expected filenames (or adjust paths in notebooks accordingly):"
echo "   - HiSeqV2_PANCAN.gz"
echo "   - BRCA_clinicalMatrix.tsv"
echo
echo "==============================================="
