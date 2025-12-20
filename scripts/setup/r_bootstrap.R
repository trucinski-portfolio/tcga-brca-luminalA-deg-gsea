# scripts/setup/r_bootstrap.R
# renv-first bootstrap: repo root discovery + (optional) renv restore

find_repo_root <- function(start = getwd()) {
  p <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in 1:20) {
    if (file.exists(file.path(p, "README.md")) &&
        dir.exists(file.path(p, "scripts")) &&
        dir.exists(file.path(p, "data"))) {
      return(p)
    }
    parent <- dirname(p)
    if (identical(parent, p)) break
    p <- parent
  }
  stop("Could not find repo root (expected README.md + scripts/ + data/).")
}

REPO_ROOT <- find_repo_root()
setwd(REPO_ROOT)

cat("Repo root:", REPO_ROOT, "\n")

# --- renv handling ---
if (file.exists(file.path(REPO_ROOT, "renv.lock"))) {

  # If renv isn't available in the current R, install it
  if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv", repos = "https://cloud.r-project.org")
  }

  # Activate renv for this project 
  renv::activate(project = REPO_ROOT)

  # Restore packages if the library is missing/out-of-sync
  status <- renv::status(project = REPO_ROOT)
  if (!status$synchronized) {
    cat("renv is out-of-sync. Restoring from renv.lock...\n")
    renv::restore(project = REPO_ROOT, prompt = FALSE)
  } else {
    cat("renv is synchronized.\n")
  }

} else {
  warning("No renv.lock found. Run renv::init() + renv::snapshot() to lock dependencies.")
}

cat("✅ R bootstrap complete.\n")