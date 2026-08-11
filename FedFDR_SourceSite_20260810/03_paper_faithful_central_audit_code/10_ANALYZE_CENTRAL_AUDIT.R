# ============================================================
# Analyze paper-faithful central-step audit
# ============================================================
rm(list = ls())
options(stringsAsFactors = FALSE)

INPUT_DIR <- "results_central_audit_R10000"
input_file <- file.path(INPUT_DIR, "AUDIT_R10000_ALL_RESULTS.rds")
if (!file.exists(input_file)) stop("Cannot find ", input_file)

res <- readRDS(input_file)
v <- res[res$status == "ok", , drop = FALSE]
f <- res[res$status == "failed", , drop = FALSE]
outdir <- file.path(INPUT_DIR, "analysis")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

mean_ci <- function(x) {
  x <- as.numeric(x[is.finite(x)])
  n <- length(x)
  m <- mean(x)
  se <- if (n > 1L) stats::sd(x) / sqrt(n) else NA_real_
  c(mean = m, lower = m - 1.96 * se, upper = m + 1.96 * se, n = n)
}

summ_variant <- function(prefix, label) {
  pair_fdp <- mean_ci(v[[paste0(prefix, "_pair_fdp")]])
  pair_power <- mean_ci(v[[paste0(prefix, "_pair_power")]])
  pair_size <- mean_ci(v[[paste0(prefix, "_pair_size")]])
  final_fdp <- mean_ci(v[[paste0(prefix, "_final_fdp")]])
  final_power <- mean_ci(v[[paste0(prefix, "_final_power")]])
  final_size <- mean_ci(v[[paste0(prefix, "_final_size")]])
  data.frame(
    variant = label,
    pair_mean_FDP = pair_fdp["mean"],
    pair_mean_Power = pair_power["mean"],
    pair_mean_Size = pair_size["mean"],
    final_mean_FDP = final_fdp["mean"],
    final_FDP_CI_low = final_fdp["lower"],
    final_FDP_CI_high = final_fdp["upper"],
    final_mean_Power = final_power["mean"],
    final_Power_CI_low = final_power["lower"],
    final_Power_CI_high = final_power["upper"],
    final_mean_Size = final_size["mean"],
    final_empty_rate = mean(v[[paste0(prefix, "_final_size")]] == 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

summary_table <- rbind(
  summ_variant("A_code_theta", "A: code theta input + code central logic"),
  summ_variant("B_code_rawbeta", "B: raw beta input + code central logic"),
  summ_variant("C_paper_rawbeta", "C: raw beta input + paper Algorithm 1")
)
utils::write.csv(summary_table, file.path(outdir, "01_variant_summary.csv"), row.names = FALSE)

comparison <- data.frame(
  metric = c(
    "pair exact-match rate A vs B (input-scale effect)",
    "pair exact-match rate B vs C (central-rule effect)",
    "pair mean Jaccard A vs B",
    "pair mean Jaccard B vs C",
    "final exact-match rate A vs B",
    "final exact-match rate B vs C",
    "final mean Jaccard A vs B",
    "final mean Jaccard B vs C",
    "paper final empty rate",
    "paper pair size >= 10 rate",
    "paper tau=0 boundary-infimum convention rate"
  ),
  value = c(
    mean(v$same_A_B_pair, na.rm = TRUE),
    mean(v$same_B_C_pair, na.rm = TRUE),
    mean(v$jaccard_A_vs_B_pair, na.rm = TRUE),
    mean(v$jaccard_B_vs_C_pair, na.rm = TRUE),
    mean(v$same_A_B_final, na.rm = TRUE),
    mean(v$same_B_C_final, na.rm = TRUE),
    mean(v$jaccard_A_vs_B_final, na.rm = TRUE),
    mean(v$jaccard_B_vs_C_final, na.rm = TRUE),
    mean(v$C_paper_rawbeta_final_size == 0, na.rm = TRUE),
    mean(v$C_paper_rawbeta_pair_size >= 10, na.rm = TRUE),
    mean(grepl("^tau0_boundary", v$C_paper_rawbeta_tau_convention), na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(comparison, file.path(outdir, "02_variant_comparisons.csv"), row.names = FALSE)

# K=2 degeneracy diagnostic: pair size -> final size for exact paper Algorithm 1
ktab <- aggregate(
  v$C_paper_rawbeta_final_size,
  by = list(pair_size = v$C_paper_rawbeta_pair_size),
  FUN = function(x) c(n = length(x), mean_final_size = mean(x), empty_rate = mean(x == 0))
)
ktab2 <- data.frame(
  pair_size = ktab$pair_size,
  n = vapply(ktab$x, function(z) z["n"], numeric(1)),
  mean_final_size = vapply(ktab$x, function(z) z["mean_final_size"], numeric(1)),
  empty_rate = vapply(ktab$x, function(z) z["empty_rate"], numeric(1))
)
utils::write.csv(ktab2, file.path(outdir, "03_K2_pair_to_final_degeneracy.csv"), row.names = FALSE)

if (nrow(f) > 0L) {
  fail_tab <- sort(table(f$failure_reason, useNA = "ifany"), decreasing = TRUE)
  fail <- data.frame(
    failure_reason = names(fail_tab),
    count = as.integer(fail_tab),
    stringsAsFactors = FALSE
  )
  utils::write.csv(fail, file.path(outdir, "04_failure_reasons.csv"), row.names = FALSE)
}

cat("\n============================================================\n")
cat("PAPER-FAITHFUL CENTRAL-STEP AUDIT ANALYSIS\n")
cat("Valid:", nrow(v), " Failed:", nrow(f), "\n")
cat("============================================================\n\n")
cat("VARIANT SUMMARY\n")
print(summary_table, row.names = FALSE, digits = 5)
cat("\nCOMPARISONS\n")
print(comparison, row.names = FALSE, digits = 5)
cat("\nAnalysis saved to:", normalizePath(outdir), "\n")
