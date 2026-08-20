# ============================================================
# Task A K-grid summary and figures
# ============================================================
# By default this reads PILOT results. Change INPUT_STAGE to "FORMAL"
# only after a formal run exists.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

source("01_SourceAttack_Kgeneral_functions.R")

# ================= USER SETTINGS =============================
INPUT_STAGE <- "PILOT"   # "QUICK", "PILOT", or "FORMAL"
K_GRID <- c(2L, 3L, 5L, 10L, 20L, 50L, 100L)
# ============================================================

stage <- toupper(INPUT_STAGE)
if (!stage %in% c("QUICK", "PILOT", "FORMAL")) stop("Unknown INPUT_STAGE")

INPUT_DIR <- switch(
  stage,
  QUICK = "results_quick_Kgrid",
  PILOT = "results_pilot_Kgrid",
  FORMAL = "results_formal_Kgrid"
)

OUTPUT_DIR <- file.path(INPUT_DIR, "summary")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

read_one_K <- function(K_now) {
  f <- file.path(INPUT_DIR, sprintf("%s_K%03d_ALL_RESULTS.csv", stage, K_now))
  if (stage == "QUICK") {
    f <- file.path(INPUT_DIR, sprintf("QUICK_K%03d_RESULTS.csv", K_now))
  }
  if (!file.exists(f)) stop("Missing result file: ", f)
  utils::read.csv(f, stringsAsFactors = FALSE)
}

attack_list <- list()
diag_list <- list()
failure_list <- list()

for (K_now in K_GRID) {
  df <- read_one_K(K_now)
  sm <- summarize_K_results(df, K_now)
  attack_list[[as.character(K_now)]] <- sm$attack
  diag_list[[as.character(K_now)]] <- sm$diagnostics

  failed <- df[df$status == "failed", , drop = FALSE]
  if (nrow(failed) == 0L) {
    failure_list[[as.character(K_now)]] <- data.frame(
      K = K_now, failure_class = "none", count = 0L, proportion_attempted = 0,
      stringsAsFactors = FALSE
    )
  } else {
    tt <- sort(table(failed$failure_class), decreasing = TRUE)
    failure_list[[as.character(K_now)]] <- data.frame(
      K = K_now,
      failure_class = names(tt),
      count = as.integer(tt),
      proportion_attempted = as.integer(tt) / nrow(df),
      stringsAsFactors = FALSE
    )
  }
}

attack_summary <- do.call(rbind, attack_list)
diagnostics <- do.call(rbind, diag_list)
failures <- do.call(rbind, failure_list)
row.names(attack_summary) <- NULL
row.names(diagnostics) <- NULL
row.names(failures) <- NULL

utils::write.csv(
  attack_summary,
  file.path(OUTPUT_DIR, sprintf("%s_KGRID_ATTACK_SUMMARY.csv", stage)),
  row.names = FALSE
)
utils::write.csv(
  diagnostics,
  file.path(OUTPUT_DIR, sprintf("%s_KGRID_DIAGNOSTICS.csv", stage)),
  row.names = FALSE
)
utils::write.csv(
  failures,
  file.path(OUTPUT_DIR, sprintf("%s_KGRID_FAILURES.csv", stage)),
  row.names = FALSE
)

cat("\n============================================================\n")
cat(stage, "K-GRID SUMMARY\n")
cat("============================================================\n")

show_cols <- c(
  "K", "condition", "member_accuracy", "random_baseline",
  "member_minus_random", "fresh_accuracy", "member_minus_fresh", "n_valid"
)
print(attack_summary[, show_cols], row.names = FALSE)

cat("\nKEY R2-(X,Y) BY K\n")
key <- attack_summary[attack_summary$condition == "R2_XY", show_cols, drop = FALSE]
print(key, row.names = FALSE)

cat("\nDIAGNOSTICS BY K\n")
print(diagnostics, row.names = FALSE)

# ------------------------------------------------------------
# Plot 1: member accuracy across K, all main attack conditions
# ------------------------------------------------------------

plot_conditions <- c("R1_X", "R2_X", "R2_XY", "R3_X", "R3_XY")
plot_dat <- attack_summary[attack_summary$condition %in% plot_conditions, , drop = FALSE]

png(file.path(OUTPUT_DIR, sprintf("%s_plot1_member_accuracy_vs_K.png", stage)),
    width = 1400, height = 900, res = 140)

ylim_top <- max(c(plot_dat$member_accuracy, plot_dat$random_baseline), na.rm = TRUE)
ylim_top <- min(1, max(0.2, ylim_top * 1.10))
plot(
  K_GRID, 1 / K_GRID,
  type = "b", pch = 1, lty = 2,
  ylim = c(0, ylim_top),
  xlab = "Number of sites K",
  ylab = "Source-site accuracy",
  main = paste0(stage, ": Member source-site accuracy vs K")
)

pch_vals <- seq_along(plot_conditions) + 1L
for (ii in seq_along(plot_conditions)) {
  cc <- plot_conditions[ii]
  dd <- plot_dat[plot_dat$condition == cc, , drop = FALSE]
  dd <- dd[order(dd$K), , drop = FALSE]
  lines(dd$K, dd$member_accuracy, type = "b", pch = pch_vals[ii], lty = ii + 1L)
}
legend(
  "topright",
  legend = c("Random 1/K", plot_conditions),
  pch = c(1, pch_vals),
  lty = c(2, seq_along(plot_conditions) + 1L),
  bty = "n"
)
dev.off()

# ------------------------------------------------------------
# Plot 2: member advantage over random baseline
# ------------------------------------------------------------

png(file.path(OUTPUT_DIR, sprintf("%s_plot2_member_minus_random.png", stage)),
    width = 1400, height = 900, res = 140)

adv_max <- max(abs(plot_dat$member_minus_random), na.rm = TRUE)
adv_lim <- max(0.02, adv_max * 1.15)
plot(
  K_GRID, rep(0, length(K_GRID)),
  type = "l", lty = 2,
  ylim = c(-adv_lim, adv_lim),
  xlab = "Number of sites K",
  ylab = "Member accuracy - 1/K",
  main = paste0(stage, ": Attack advantage above random guessing")
)
for (ii in seq_along(plot_conditions)) {
  cc <- plot_conditions[ii]
  dd <- plot_dat[plot_dat$condition == cc, , drop = FALSE]
  dd <- dd[order(dd$K), , drop = FALSE]
  lines(dd$K, dd$member_minus_random, type = "b", pch = pch_vals[ii], lty = ii + 1L)
}
legend(
  "topright",
  legend = plot_conditions,
  pch = pch_vals,
  lty = seq_along(plot_conditions) + 1L,
  bty = "n"
)
dev.off()

# ------------------------------------------------------------
# Plot 3: R2-(X,Y), member vs matched fresh vs random
# ------------------------------------------------------------

r2 <- attack_summary[attack_summary$condition == "R2_XY", , drop = FALSE]
r2 <- r2[order(r2$K), , drop = FALSE]

png(file.path(OUTPUT_DIR, sprintf("%s_plot3_R2XY_member_fresh_random.png", stage)),
    width = 1400, height = 900, res = 140)

r2_top <- max(c(r2$member_accuracy, r2$fresh_accuracy, r2$random_baseline), na.rm = TRUE)
plot(
  r2$K, r2$member_accuracy,
  type = "b", pch = 16,
  ylim = c(0, min(1, max(0.2, r2_top * 1.12))),
  xlab = "Number of sites K",
  ylab = "Source-site accuracy",
  main = paste0(stage, ": R2-(X,Y) member vs matched fresh")
)
lines(r2$K, r2$fresh_accuracy, type = "b", pch = 1, lty = 2)
lines(r2$K, r2$random_baseline, type = "b", pch = 2, lty = 3)
legend(
  "topright",
  legend = c("Member", "Matched fresh", "Random 1/K"),
  pch = c(16, 1, 2), lty = c(1, 2, 3), bty = "n"
)
dev.off()

# ------------------------------------------------------------
# Plot 4: external-support size and failure rate (separate figures)
# ------------------------------------------------------------

png(file.path(OUTPUT_DIR, sprintf("%s_plot4_external_support_size.png", stage)),
    width = 1400, height = 900, res = 140)
plot(
  diagnostics$K, diagnostics$mean_external_support_size,
  type = "b", pch = 16,
  ylim = c(0, max(50, diagnostics$mean_external_support_size, na.rm = TRUE)),
  xlab = "Number of sites K",
  ylab = "Mean external-support size",
  main = paste0(stage, ": External-support growth with K")
)
abline(h = 50, lty = 2)
dev.off()

png(file.path(OUTPUT_DIR, sprintf("%s_plot5_failure_rate.png", stage)),
    width = 1400, height = 900, res = 140)
plot(
  diagnostics$K, diagnostics$failure_rate,
  type = "b", pch = 16,
  ylim = c(0, max(0.01, diagnostics$failure_rate, na.rm = TRUE) * 1.15),
  xlab = "Number of sites K",
  ylab = "Failure rate",
  main = paste0(stage, ": Invalid-repetition rate vs K")
)
dev.off()

cat("\nSummary files and figures saved in:\n  ", normalizePath(OUTPUT_DIR), "\n", sep = "")
cat("\nInterpretation priority:\n")
cat("  1) R2-(X,Y) member accuracy vs 1/K\n")
cat("  2) R2-(X,Y) member - matched fresh\n")
cat("  3) external-support growth / Hessian stabilization / failures\n")
