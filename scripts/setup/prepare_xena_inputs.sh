#!/usr/bin/env bash

TARGET_DIR="data/raw/preprocessing_inputs"
mkdir -p "$TARGET_DIR"

EXPECTED_FILES=(
  "HiSeqV2_PANCAN.gz"
  "BRCA_clinicalMatrix.tsv"
)

echo "Checking for required Xena input files..."

missing=0
for f in "${EXPECTED_FILES[@]}"; do
  if [[ ! -f "$TARGET_DIR/$f" ]]; then
    echo "❌ Missing: $TARGET_DIR/$f"
    missing=1
  else
    echo "✅ Found: $f"
  fi
done

if [[ $missing -eq 1 ]]; then
  echo
  echo "One or more required files are missing."
  echo "See data/README.md for download instructions."
  exit 1
fi

echo
echo "All required Xena input files are present."
