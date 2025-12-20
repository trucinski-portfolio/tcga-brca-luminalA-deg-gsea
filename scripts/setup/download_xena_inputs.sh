#!/usr/bin/env bash
set -euo pipefail

TARGET_DIR="data/raw/preprocessing_inputs"
mkdir -p "$TARGET_DIR"

# Direct Xena-hosted files (links will need update if Xena changes their paths)
EXPR_URL="https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHiSeqV2_PANCAN.gz"
CLIN_URL="https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix"

EXPR_OUT="$TARGET_DIR/HiSeqV2_PANCAN.gz"
CLIN_OUT="$TARGET_DIR/BRCA_clinicalMatrix.tsv"

echo "Target directory: $TARGET_DIR"
echo

# Pick a downloader (curl preferred, wget ok)
download() {
  local url="$1"
  local out="$2"

  if command -v curl >/dev/null 2>&1; then
    echo "Downloading: $url"
    curl -L --fail --retry 3 --retry-delay 2 -o "$out" "$url"
  elif command -v wget >/dev/null 2>&1; then
    echo "Downloading: $url"
    wget -O "$out" "$url"
  else
    echo "ERROR: Need curl or wget installed to download."
    exit 1
  fi
}

# Don’t re-download if files already exist
if [[ -f "$EXPR_OUT" ]]; then
  echo "✅ Already exists: $EXPR_OUT"
else
  download "$EXPR_URL" "$EXPR_OUT"
fi

if [[ -f "$CLIN_OUT" ]]; then
  echo "✅ Already exists: $CLIN_OUT"
else
  download "$CLIN_URL" "$CLIN_OUT"
fi

echo
echo "Verifying outputs..."
[[ -s "$EXPR_OUT" ]] || { echo "❌ Download failed or empty file: $EXPR_OUT"; exit 1; }
[[ -s "$CLIN_OUT" ]] || { echo "❌ Download failed or empty file: $CLIN_OUT"; exit 1; }

echo "✅ Done. Files present:"
ls -lh "$EXPR_OUT" "$CLIN_OUT"