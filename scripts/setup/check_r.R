# scripts/setup/setup_r.R
# Purpose: renv-based reproducible setup (no global installs)

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Ensure renv is available (minimal side effect) 
if (!requireNamespace("renv", quietly = TRUE)) {
  stop(
    "❌ Package 'renv' is not installed.\n",
    "Install it once, then rerun:\n\n",
    "  install.packages('renv', repos='https://cloud.r-project.org')\n",
    "  renv::restore()\n"
  )
}

# Ensure this is an renv project 
if (!file.exists("renv.lock")) {
  stop(
    "❌ renv.lock not found in repo root.\n",
    "This project uses renv. If you haven't created the lockfile yet, run:\n\n",
    "  renv::snapshot()\n"
  )
}

# Restore project library (this installs into ./renv/library/...) 
# Keep it non-interactive and deterministic.
renv::restore(prompt = FALSE)

# Required packages should now be available 
required_pkgs <- c(
  # CRAN
  "cowplot", "data.table", "dplyr", "ggplot2", "knitr", "readr", "stringr",
  # Bioconductor
  "limma", "ComplexHeatmap", "circlize"
)

missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing)) {
  stop(
    "❌ After renv::restore(), some packages are still missing:\n  - ",
    paste(missing, collapse = "\n  - "),
    "\n\nTry:\n  renv::restore(rebuild = TRUE, prompt = FALSE)\n"
  )
}

# Print a short environment summary (useful in notebook output)
message("✅ renv restore complete.")
message("R version: ", getRversion())
message("Library: ", .libPaths()[1])
message("Lockfile: ", normalizePath("renv.lock"))