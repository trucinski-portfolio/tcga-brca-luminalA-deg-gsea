# R/07_survival.R
# Stage 07: Cox proportional-hazards survival analysis.
#           For each PAM50 subtype, scores each tumor sample on its ElasticNet
#           signature (mean log2-TPM of signature genes), then fits a Cox model
#           against Overall Survival. HR tables saved per subtype.
#
# Inputs:  data/processed/cohort.rds
#          results/tables/signatures/{Subtype}_signature.csv   (4 files)
# Outputs: results/tables/survival/{Subtype}_cox.csv           (4 files)
#            Columns: subtype, gene_score_hr, gene_score_lower95,
#                     gene_score_upper95, gene_score_pval, n_events, n_total
#          results/figures/survival_km_{Subtype}.{png,pdf}     (4 pairs)

library(here)
library(dplyr)
library(survival)
library(survminer)
library(UCSCXenaTools)
library(ggplot2)

set.seed(42)

PAM50_LEVELS <- c("LumA", "LumB", "Her2", "Basal")

subtype_colors <- c(
  LumA  = "#2166AC",
  LumB  = "#92C5DE",
  Her2  = "#D6604D",
  Basal = "#1A1A1A"
)

# ---------------------------------------------------------------------------
# Stage gate
# ---------------------------------------------------------------------------

#' Stop with a clear error if Stage 04 outputs are missing
#'
#' @return Invisible NULL.
#' @examples
#' assert_stage04_complete()
assert_stage04_complete <- function() {
  message("── Stage 07 | stage gate ────────────────────────────────────────────────")
  cohort_path <- here::here("data", "processed", "cohort.rds")
  if (!file.exists(cohort_path))
    stop("cohort.rds missing — run R/01_preprocessing.R first.")

  sig_dir  <- here::here("results", "tables", "signatures")
  required <- paste0(PAM50_LEVELS, "_signature.csv")
  missing  <- required[!file.exists(file.path(sig_dir, required))]
  if (length(missing) > 0)
    stop("Missing signature files: ", paste(missing, collapse = ", "),
         "\nRun R/04_feature_selection.R first.")

  message("  cohort.rds found.")
  message("  signature files found for all subtypes.")
  message("── Stage 07 | gate passed ───────────────────────────────────────────────")
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Fetch survival data
# ---------------------------------------------------------------------------

#' Download TCGA-BRCA overall survival data from UCSC Xena (tcgaHub)
#'
#' Fetches the TCGA BRCA clinical matrix which contains OS_time and OS
#' (overall survival event indicator). Caches to data/processed/survival.rds.
#' Uses only columns needed: sampleID, OS, OS.time.
#'
#' @return data.frame with columns: sample_id, os_time (days), os_event (0/1).
#' @examples
#' surv_data <- fetch_survival()
fetch_survival <- function() {
  message("── Stage 07 | fetching survival data ────────────────────────────────────")
  cache_path <- here::here("data", "processed", "survival.rds")

  if (file.exists(cache_path)) {
    message("  loading cached survival data ...")
    surv <- readRDS(cache_path)
    message("  ", nrow(surv), " samples with survival data")
    message("── Stage 07 | survival data loaded ──────────────────────────────────────")
    return(surv)
  }

  message("  downloading TCGA BRCA clinical matrix from tcgaHub ...")
  dl <- XenaGenerate(subset = XenaHostNames == "tcgaHub") |>
    XenaFilter(filterDatasets = "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix") |>
    XenaQuery() |>
    XenaDownload()

  clin <- XenaPrepare(dl)
  message("  clinical matrix dimensions: ", nrow(clin), " × ", ncol(clin))

  # Extract OS columns — confirmed column names in BRCA_clinicalMatrix:
  # OS_Time_nature2012 (days) and OS_event_nature2012 (0/1)
  required_cols <- c("sampleID", "OS_Time_nature2012", "OS_event_nature2012")
  if (!all(required_cols %in% colnames(clin)))
    stop("Expected columns not found. Available: ",
         paste(head(colnames(clin), 20), collapse = ", "))

  surv <- clin |>
    dplyr::select(sample_id = sampleID,
                  os_time   = OS_Time_nature2012,
                  os_event  = OS_event_nature2012) |>
    dplyr::filter(!is.na(os_time), !is.na(os_event)) |>
    dplyr::mutate(
      os_time  = as.numeric(os_time),
      os_event = as.integer(as.character(os_event))
    ) |>
    dplyr::filter(os_time > 0)

  message("  samples with valid OS data: ", nrow(surv))
  message("  events (deaths): ", sum(surv$os_event == 1))

  dir.create(here::here("data", "processed"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(surv, cache_path)
  message("  cached: survival.rds")
  message("── Stage 07 | survival data loaded ──────────────────────────────────────")
  surv
}

# ---------------------------------------------------------------------------
# Signature score
# ---------------------------------------------------------------------------

#' Compute a single-sample signature score as mean expression of signature genes
#'
#' For each tumor sample, averages the log2-TPM values across all signature
#' genes that are present in the expression matrix. Returns a named numeric
#' vector (sample → score).
#'
#' @param expr    Numeric matrix. HGNC symbols × sample IDs.
#' @param sig_csv Character. Path to signature CSV (gene, selection_frequency, mean_coef).
#' @return Named numeric vector. Names = sample IDs, values = mean log2-TPM.
#' @examples
#' scores <- compute_signature_score(expr, here::here("results","tables","signatures","LumA_signature.csv"))
compute_signature_score <- function(expr, sig_csv) {
  sig   <- read.csv(sig_csv, stringsAsFactors = FALSE)
  genes <- sig$gene
  found <- intersect(genes, rownames(expr))

  if (length(found) == 0)
    stop("None of the signature genes found in expression matrix.")
  if (length(found) < length(genes))
    message("    ", length(genes) - length(found), " signature genes absent from expr — using ",
            length(found), " of ", length(genes))

  colMeans(expr[found, , drop = FALSE])
}

# ---------------------------------------------------------------------------
# Cox model
# ---------------------------------------------------------------------------

#' Fit a univariate Cox model: signature score → Overall Survival
#'
#' Subsets to Tumor samples of the specified PAM50 subtype, joins expression
#' scores with survival data, and fits a Cox PH model. Returns a one-row
#' data.frame with HR, 95% CI, p-value, and event counts.
#'
#' @param expr      Numeric matrix. HGNC symbols × sample IDs.
#' @param meta      data.frame. Columns: sample_id, group, PAM50.
#' @param surv_data data.frame. Columns: sample_id, os_time, os_event.
#' @param subtype   Character. One of PAM50_LEVELS.
#' @return data.frame: subtype, gene_score_hr, gene_score_lower95,
#'   gene_score_upper95, gene_score_pval, n_events, n_total.
#' @examples
#' result <- run_cox(expr, meta, surv_data, "LumA")
run_cox <- function(expr, meta, surv_data, subtype) {
  message("  Cox model: ", subtype)

  sig_path <- here::here("results", "tables", "signatures",
                          paste0(subtype, "_signature.csv"))
  scores <- compute_signature_score(
    expr[, meta$sample_id[meta$group == "Tumor" & !is.na(meta$PAM50) &
                             meta$PAM50 == subtype], drop = FALSE],
    sig_path
  )

  df <- data.frame(
    sample_id = names(scores),
    score     = as.numeric(scores),
    stringsAsFactors = FALSE
  ) |>
    dplyr::inner_join(surv_data, by = "sample_id") |>
    dplyr::filter(!is.na(os_time), !is.na(os_event), os_time > 0)

  n_total  <- nrow(df)
  n_events <- sum(df$os_event == 1)

  message("    n = ", n_total, " samples, events = ", n_events)

  if (n_events < 5) {
    warning(subtype, ": only ", n_events, " OS events — Cox estimates unreliable.")
  }

  # Scale score to SD units for interpretable HR
  df$score_scaled <- scale(df$score)[, 1]

  fit <- survival::coxph(
    survival::Surv(os_time, os_event) ~ score_scaled,
    data = df
  )

  coef_summary <- summary(fit)$coefficients
  conf_int     <- summary(fit)$conf.int

  data.frame(
    subtype            = subtype,
    gene_score_hr      = conf_int[1, "exp(coef)"],
    gene_score_lower95 = conf_int[1, "lower .95"],
    gene_score_upper95 = conf_int[1, "upper .95"],
    gene_score_pval    = coef_summary[1, "Pr(>|z|)"],
    n_events           = n_events,
    n_total            = n_total,
    stringsAsFactors   = FALSE,
    row.names          = NULL
  )
}

# ---------------------------------------------------------------------------
# Kaplan-Meier plot
# ---------------------------------------------------------------------------

#' Draw and save a Kaplan-Meier survival plot split by signature score median
#'
#' Dichotomises each tumor sample into High/Low signature expression groups
#' using the within-subtype median score as threshold. Produces a KM curve
#' with risk table, log-rank p-value, and saves as .png + .pdf.
#'
#' @param expr      Numeric matrix. HGNC symbols × sample IDs.
#' @param meta      data.frame. Columns: sample_id, group, PAM50.
#' @param surv_data data.frame. Columns: sample_id, os_time, os_event.
#' @param subtype   Character. One of PAM50_LEVELS.
#' @return Invisible NULL.
#' @examples
#' plot_km(expr, meta, surv_data, "LumA")
plot_km <- function(expr, meta, surv_data, subtype) {
  message("  KM plot: ", subtype)

  sig_path <- here::here("results", "tables", "signatures",
                          paste0(subtype, "_signature.csv"))
  tumor_ids <- meta$sample_id[meta$group == "Tumor" & !is.na(meta$PAM50) &
                                 meta$PAM50 == subtype]

  scores <- compute_signature_score(expr[, tumor_ids, drop = FALSE], sig_path)

  df <- data.frame(
    sample_id = names(scores),
    score     = as.numeric(scores),
    stringsAsFactors = FALSE
  ) |>
    dplyr::inner_join(surv_data, by = "sample_id") |>
    dplyr::filter(!is.na(os_time), !is.na(os_event), os_time > 0) |>
    dplyr::mutate(
      score_group = ifelse(score >= median(score), "High", "Low"),
      score_group = factor(score_group, levels = c("Low", "High"))
    )

  fit <- survminer::surv_fit(
    survival::Surv(os_time, os_event) ~ score_group,
    data = df
  )

  p <- survminer::ggsurvplot(
    fit,
    data           = df,
    pval           = TRUE,
    pval.method    = TRUE,
    conf.int       = TRUE,
    risk.table     = TRUE,
    risk.table.height = 0.25,
    palette        = c("steelblue", subtype_colors[[subtype]]),
    title          = paste0(subtype, " — Signature Score vs Overall Survival"),
    xlab           = "Time (days)",
    ylab           = "Overall Survival Probability",
    legend.title   = "Signature",
    legend.labs    = c("Low", "High"),
    ggtheme        = theme_bw(base_size = 12)
  )

  print(p)

  out_base <- here::here("results", "figures",
                          paste0("survival_km_", subtype))
  dir.create(here::here("results", "figures"), recursive = TRUE,
             showWarnings = FALSE)

  png(paste0(out_base, ".png"), width = 2400, height = 2000, res = 300)
  print(p)
  dev.off()

  pdf(paste0(out_base, ".pdf"), width = 10, height = 9)
  print(p)
  dev.off()

  message("    saved: survival_km_", subtype, ".png + .pdf")
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

#' Write Cox HR results to results/tables/survival/{Subtype}_cox.csv
#'
#' @param result  data.frame. One-row output of \code{run_cox()}.
#' @param subtype Character. Used as filename stem.
#' @return Invisible character path.
#' @examples
#' save_cox(result, "LumA")
save_cox <- function(result, subtype) {
  out_dir <- here::here("results", "tables", "survival")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, paste0(subtype, "_cox.csv"))
  write.csv(result, file = path, row.names = FALSE, quote = FALSE)
  message("  saved: ", basename(path),
          " (HR = ", round(result$gene_score_hr, 3),
          ", p = ", signif(result$gene_score_pval, 3), ")")
  invisible(path)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

assert_stage04_complete()

message("── Stage 07 | loading cohort ────────────────────────────────────────────")
cohort <- readRDS(here::here("data", "processed", "cohort.rds"))
expr   <- cohort$expr
meta   <- cohort$meta
message("  expr: ", nrow(expr), " genes × ", ncol(expr), " samples")
message("  meta: ", nrow(meta), " rows")
message("── Stage 07 | cohort loaded ─────────────────────────────────────────────")

surv_data <- fetch_survival()

message("── Stage 07 | fitting Cox models ────────────────────────────────────────")
cox_results <- lapply(PAM50_LEVELS, function(st) {
  result <- run_cox(expr, meta, surv_data, st)
  save_cox(result, st)
  result
})
all_cox <- dplyr::bind_rows(cox_results)

message("── Stage 07 | summary ───────────────────────────────────────────────────")
print(all_cox[, c("subtype", "gene_score_hr", "gene_score_lower95",
                   "gene_score_upper95", "gene_score_pval", "n_events")])

message("── Stage 07 | KM plots ──────────────────────────────────────────────────")
for (st in PAM50_LEVELS) {
  plot_km(expr, meta, surv_data, st)
}

# Verify all 4 cox files written
cox_files <- list.files(here::here("results", "tables", "survival"),
                         pattern = "_cox\\.csv$")
message("── Stage 07 | output summary ────────────────────────────────────────────")
message("  Cox CSV files: ", length(cox_files), " / 4 expected")
stopifnot("expected 4 cox files" = length(cox_files) == 4L)

n_sig <- sum(all_cox$gene_score_pval < 0.05, na.rm = TRUE)
message("  subtypes with significant OS association (p<0.05): ", n_sig, " / 4")

message("✅ Stage 07 complete. Proceed to R/08_visualization.R")
