#!/bin/bash
# Render Jupyter notebooks to HTML reports
#
# Usage:
#   ./scripts/setup/render_reports.sh
#
# Prerequisites:
#   Option 1: Quarto (recommended)
#     brew install quarto
#
#   Option 2: jupyter nbconvert (fallback)
#     pip install nbconvert
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

cd "$REPO_ROOT"

# Create output directory
mkdir -p docs/reports

echo "=== Rendering Jupyter Notebooks to HTML ==="
echo ""

# Check which tool is available
if command -v quarto &> /dev/null; then
    echo "Using Quarto for rendering..."
    echo ""

    # Render each notebook
    for nb in scripts/notebooks/0*.ipynb; do
        echo "Rendering: $nb"
        quarto render "$nb" --to html --output-dir docs/reports
    done

elif command -v jupyter &> /dev/null && jupyter nbconvert --version &> /dev/null; then
    echo "Quarto not found, using jupyter nbconvert..."
    echo ""

    # Render each notebook with nbconvert
    for nb in scripts/notebooks/0*.ipynb; do
        echo "Rendering: $nb"
        jupyter nbconvert --to html --output-dir "$REPO_ROOT/docs/reports" "$nb"
    done

else
    echo "ERROR: Neither Quarto nor jupyter nbconvert found."
    echo ""
    echo "Install one of the following:"
    echo "  Quarto (recommended): brew install quarto"
    echo "  jupyter nbconvert:    pip install nbconvert"
    exit 1
fi

echo ""
echo "=== Rendering Complete ==="
echo "HTML reports saved to: docs/reports/"
ls -la docs/reports/*.html 2>/dev/null || echo "No HTML files generated"
