# R/04_feature_selection.R
# Stage 04: ElasticNet bootstrap feature selection — one minimal gene signature
#           per PAM50 subtype.
#
# Inputs:  data/processed/cohort.rds
#          results/tables/deg/*.csv   (18 files — all genes, no pre-filter)
# Outputs: results/tables/signatures/{Subtype}_signature.csv   (4 files)
#            Columns: gene, selection_frequency, mean_coef
#
# Method:
#   For each PAM50 subtype (1-vs-rest among Tumor samples):
#     1. Run cv.glmnet (alpha=0.5, 10-fold) once on the full data to find lambda.1se
#     2. Bootstrap 100 iterations: resample with replacement, fit glmnet at lambda.1se
#     3. Retain genes selected in >= 70% of bootstrap runs
#   Target: 10-30 genes per subtype. Prints a warning if outside that range.

library(here)
library(dplyr)
library(glmnet)

set.seed(42)

PAM50_LEVELS  <- c("LumA", "LumB", "Her2", "Basal")
ALPHA         <- 0.5      # ElasticNet mixing: 0 = Ridge, 1 = Lasso
N_FOLDS       <- 10L
N_BOOT        <- 100L
FREQ_CUTOFF   <- 0.80     # minimum selection frequency to include a gene

# ---------------------------------------------------------------------------
# Stage gate
# ---------------------------------------------------------------------------

#' Stop with a clear error if Stage 03 outputs are missing or incomplete
#'
#' Checks for cohort.rds and all 18 DEG .csv files.
#'
#' @return Invisible NULL.
#' @examples
#' assert_stage03_complete()
assert_stage03_complete <- function() {
  message("── Stage 04 | stage gate ────────────────────────────────────────────────")

  cohort_path <- here::here("data", "processed", "cohort.rds")
  if (!file.exists(cohort_path))
    stop("cohort.rds missing — run R/01_preprocessing.R first.")

  deg_dir <- here::here("results", "tables", "deg")
  expected <- paste0(c(
    paste0(PAM50_LEVELS, "_vs_NAT"),
    paste0(PAM50_LEVELS, "_vs_GTEx"),
    "LumA_vs_Basal", "LumA_vs_LumB", "LumA_vs_Her2",
    "Her2_vs_Basal", "LumB_vs_Basal", "NAT_vs_GTEx",
    paste0("Within_", PAM50_LEVELS)
  ), ".csv")
  missing <- expected[!file.exists(file.path(deg_dir, expected))]
  if (length(missing) > 0)
    stop("Missing DEG files: ", paste(missing, collapse = ", "),
         "\nRun R/03_deg.R first.")

  message("  cohort.rds found.")
  message("  18 DEG files found.")
  message("── Stage 04 | gate passed ───────────────────────────────────────────────")
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

#' Load cohort.rds and return expr and meta
#'
#' @return Named list: \code{$expr} (matrix), \code{$meta} (data.frame).
#' @examples
#' inputs <- load_cohort()
load_cohort <- function() {
  message("── Stage 04 | loading cohort ────────────────────────────────────────────")
  cohort <- readRDS(here::here("data", "processed", "cohort.rds"))
  message("  $expr : ", nrow(cohort$expr), " genes × ", ncol(cohort$expr), " samples")
  message("  $meta : ", nrow(cohort$meta), " rows")
  message("── Stage 04 | cohort loaded ─────────────────────────────────────────────")
  list(expr = cohort$expr, meta = cohort$meta)
}

# ---------------------------------------------------------------------------
# ElasticNet bootstrap
# ---------------------------------------------------------------------------

#' Run ElasticNet bootstrap stability selection for one PAM50 subtype
#'
#' Subsets expression to Tumor samples only and builds a binary outcome
#' (1 = this subtype, 0 = all other tumor subtypes). Runs \code{cv.glmnet}
#' once to determine \code{lambda.1se}, then bootstraps \code{n_boot} times:
#' each iteration resamples observations with replacement and fits
#' \code{glmnet} at the fixed \code{lambda.1se}. Genes are retained if their
#' selection frequency across bootstrap runs meets \code{freq_cutoff}.
#'
#' Using all 28,344 genes from the expression matrix as features — no
#' DEG significance pre-filter — lets ElasticNet perform its own selection
#' from the full gene space.
#'
#' @param expr        Numeric matrix. HGNC symbols × sample IDs.
#' @param meta        data.frame. Columns: sample_id, group, PAM50.
#' @param subtype     Character. One of PAM50_LEVELS.
#' @param n_boot      Integer. Bootstrap iterations. Default 100.
#' @param freq_cutoff Numeric. Minimum selection frequency [0,1]. Default 0.70.
#' @return data.frame: gene, selection_frequency, mean_coef. Sorted by
#'   descending selection_frequency.
#' @examples
#' sig <- run_elasticnet_subtype(expr, meta, "LumA")
run_elasticnet_subtype <- function(expr, meta,
                                   subtype,
                                   n_boot      = N_BOOT,
                                   freq_cutoff = FREQ_CUTOFF) {
  message("── Stage 04 | ", subtype, " ───────────────────────────────────────────")

  # Tumor samples only — feature selection is within the tumour compartment
  tumor_meta <- meta[meta$group == "Tumor" & !is.na(meta$PAM50), ]
  tumor_expr <- expr[, tumor_meta$sample_id, drop = FALSE]

  y <- as.integer(tumor_meta$PAM50 == subtype)
  X <- t(tumor_expr)   # samples × genes — glmnet convention

  message("  samples: ", nrow(X), " (", sum(y == 1L), " ", subtype,
          " + ", sum(y == 0L), " others)")
  message("  features: ", ncol(X), " genes")

  # Step 1: cv.glmnet once on all genes to find lambda.1se and lambda.min
  message("  running cv.glmnet (", N_FOLDS, "-fold) to determine lambda.1se ...")
  set.seed(42)
  cv_fit     <- glmnet::cv.glmnet(X, y, alpha = ALPHA, family = "binomial",
                                   nfolds = N_FOLDS, type.measure = "deviance")
  lambda_1se <- cv_fit$lambda.1se
  n_nonzero  <- sum(as.numeric(as.matrix(coef(cv_fit, s = "lambda.1se"))) != 0) - 1L
  message("  lambda.1se = ", signif(lambda_1se, 4),
          " | non-zero coefs at lambda.1se: ", n_nonzero)

  # Reduce bootstrap feature space to genes active at lambda.min.
  # Access beta matrix directly by index — avoids coef() interpolation issues
  # on sparse dgCMatrix objects for imbalanced binary outcomes.
  lmin_idx     <- which.min(abs(cv_fit$lambda - cv_fit$lambda.min))
  beta_min_col <- cv_fit$glmnet.fit$beta[, lmin_idx]
  active_genes <- rownames(cv_fit$glmnet.fit$beta)[beta_min_col != 0]

  # Fallback: if lambda.min path is degenerate, use lambda.1se genes (still fast)
  if (length(active_genes) < 2L) {
    message("  lambda.min path degenerate — falling back to lambda.1se gene set")
    l1se_idx     <- which.min(abs(cv_fit$lambda - cv_fit$lambda.1se))
    beta_1se_col <- cv_fit$glmnet.fit$beta[, l1se_idx]
    active_genes <- rownames(cv_fit$glmnet.fit$beta)[beta_1se_col != 0]
  }

  X_boot <- X[, active_genes, drop = FALSE]
  message("  bootstrap feature space: ", ncol(X_boot), " genes")

  # Step 2: bootstrap stability selection at fixed lambda.1se on reduced feature space
  message("  bootstrapping (n = ", n_boot, ") ...")
  n_samples  <- nrow(X_boot)
  gene_names <- colnames(X_boot)

  sel_matrix  <- matrix(0L,  nrow = length(gene_names), ncol = n_boot,
                         dimnames = list(gene_names, NULL))
  coef_matrix <- matrix(0.0, nrow = length(gene_names), ncol = n_boot,
                         dimnames = list(gene_names, NULL))

  for (b in seq_len(n_boot)) {
    idx   <- sample(n_samples, n_samples, replace = TRUE)
    fit_b <- glmnet::glmnet(X_boot[idx, ], y[idx], alpha = ALPHA,
                             family = "binomial", lambda = lambda_1se)
    coef_b <- as.numeric(coef(fit_b))[-1L]   # drop intercept
    sel_matrix[, b]  <- as.integer(coef_b != 0)
    coef_matrix[, b] <- coef_b
  }

  selection_freq <- rowMeans(sel_matrix)
  mean_coef      <- rowMeans(coef_matrix)

  # Step 3: retain genes meeting frequency cutoff
  keep <- selection_freq >= freq_cutoff
  result <- data.frame(
    gene                = gene_names[keep],
    selection_frequency = selection_freq[keep],
    mean_coef           = mean_coef[keep],
    stringsAsFactors    = FALSE,
    row.names           = NULL
  )
  result <- result[order(-result$selection_frequency), ]

  n_sig <- nrow(result)
  message("  genes retained (freq >= ", freq_cutoff, "): ", n_sig)
  if (n_sig < 10L)
    warning(subtype, ": only ", n_sig, " genes selected — consider lowering FREQ_CUTOFF.")
  if (n_sig > 30L)
    message("  note: ", n_sig, " genes exceeds target of 30 — signature is broader than expected.")

  message("── Stage 04 | ", subtype, " complete ──────────────────────────────────")
  result
}

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

#' Write one subtype signature to results/tables/signatures/{Subtype}_signature.csv
#'
#' @param sig     data.frame. Output of \code{run_elasticnet_subtype()}.
#' @param subtype Character. Used as filename stem.
#' @return Invisible character path.
#' @examples
#' save_signature(sig, "LumA")
save_signature <- function(sig, subtype) {
  out_dir <- here::here("results", "tables", "signatures")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, paste0(subtype, "_signature.csv"))
  write.csv(sig, file = path, row.names = FALSE, quote = FALSE)
  message("  saved: ", path, " (", nrow(sig), " genes)")
  invisible(path)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

assert_stage03_complete()

inputs <- load_cohort()
expr   <- inputs$expr
meta   <- inputs$meta

for (st in PAM50_LEVELS) {
  out_path <- here::here("results", "tables", "signatures",
                          paste0(st, "_signature.csv"))
  if (file.exists(out_path)) {
    message("── Stage 04 | ", st, " — already complete, skipping ─────────────────")
    next
  }
  sig <- run_elasticnet_subtype(expr, meta, st)
  save_signature(sig, st)
}

# Verify 4 signature files written
sig_files <- list.files(here::here("results", "tables", "signatures"),
                         pattern = "_signature\\.csv$")
message("── Stage 04 | output summary ────────────────────────────────────────────")
message("  signature files written: ", length(sig_files), " / 4 expected")
stopifnot("expected 4 signature files" = length(sig_files) == 4L)

message("✅ Stage 04 complete. Proceed to R/05_fgsea.R")
