#!/usr/bin/env bash
set -euo pipefail

# =========================================
# TCGA pipeline - Python environment setup
# - Creates .venv
# - Installs pinned deps (env/requirements_python.txt)
# - Registers a Jupyter kernel
# - Verifies key imports
# - Optionally writes an exact lockfile (env/requirements_python.lock)
# =========================================

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

VENV_DIR="${VENV_DIR:-.venv}"
REQ_FILE="${REQ_FILE:-env/requirements_python.txt}"
LOCK_FILE="${LOCK_FILE:-env/requirements_python.lock}"
KERNEL_NAME="${KERNEL_NAME:-tcga-brca-luma}"
KERNEL_DISPLAY="${KERNEL_DISPLAY:-Python (.venv) - tcga-brca-luminalA-deg-gsea}"
WRITE_LOCK="${WRITE_LOCK:-1}"   # set to 0 to skip lockfile write

echo "Repo root : $REPO_ROOT"
echo "Venv dir  : $VENV_DIR"
echo "Req file  : $REQ_FILE"
echo "Lock file : $LOCK_FILE"
echo

# ---- 0) Choose interpreter ----
if command -v python3 >/dev/null 2>&1; then
  PY=python3
elif command -v python >/dev/null 2>&1; then
  PY=python
else
  echo "❌ No python found on PATH."
  exit 1
fi

# Fail early if venv module isn't available
if ! $PY -c "import venv" >/dev/null 2>&1; then
  echo "❌ Selected Python cannot create venvs (missing 'venv' module)."
  echo "   Python used: $($PY -c 'import sys; print(sys.executable)')"
  echo "   Fix: install a proper Python 3 build, then re-run."
  exit 1
fi

echo "Using interpreter: $($PY -c 'import sys; print(sys.executable)')"
echo "Python version   : $($PY --version)"
echo

# ---- 1) Ensure requirements file exists ----
if [[ ! -f "$REQ_FILE" ]]; then
  echo "❌ Requirements file not found: $REQ_FILE"
  echo "   This repo expects pinned dependencies there."
  echo "   Create it (or commit it) and re-run."
  exit 1
fi

# ---- 2) Create venv if missing ----
if [[ ! -d "$VENV_DIR" ]]; then
  echo "Creating venv..."
  $PY -m venv "$VENV_DIR"
else
  echo "Venv already exists."
fi

# ---- 3) Activate venv ----
# shellcheck disable=SC1090
source "$VENV_DIR/bin/activate"

echo
echo "✅ Activated venv: $(python -c 'import sys; print(sys.executable)')"
echo

# ---- 4) Upgrade pip tooling ----
python -m pip install --upgrade pip setuptools wheel

# ---- 5) Install dependencies ----
echo
echo "Installing requirements from $REQ_FILE ..."
python -m pip install -r "$REQ_FILE"

# ---- 6) Ensure Jupyter kernel support & register kernel ----
python -m pip install --upgrade ipykernel
python -m ipykernel install --user --name "$KERNEL_NAME" --display-name "$KERNEL_DISPLAY"

# ---- 7) Verification ----
echo
echo "Running verification imports..."
python - <<'PY'
import sys
print("Python:", sys.executable)

mods = ["numpy", "pandas", "ipykernel"]
for m in mods:
    __import__(m)
    mod = sys.modules[m]
    ver = getattr(mod, "__version__", "unknown")
    print(f"{m}: {ver}")

print("✅ Verification OK")
PY

# ---- 8) Optional: write exact lockfile ----
if [[ "$WRITE_LOCK" == "1" ]]; then
  echo
  echo "Writing lockfile to $LOCK_FILE ..."
  python -m pip freeze --all | sort > "$LOCK_FILE"
  echo "✅ Wrote: $LOCK_FILE"
fi

echo
echo "========================================="
echo "✅ Python environment setup complete."
echo "Next:"
echo "  - VS Code: Command Palette -> 'Python: Select Interpreter' -> choose .venv"
echo "  - Notebook kernel: '$KERNEL_DISPLAY'"
echo "========================================="