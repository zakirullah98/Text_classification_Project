# ============================================================
# Task A FORMAL N=10,000 summary, teacher table, and figures
# ============================================================
# Reads: results_formal_Kgrid_N10000/FORMAL_Kxxx_ALL_RESULTS.csv
# Produces:
#   - teacher-facing main table for R2-(X,Y) and R3-(X,Y)
#   - all-attack table with Monte Carlo uncertainty
#   - diagnostics table
#   - formal figures for the K-sensitivity result
#
# IMPORTANT: this script refuses to summarize a K unless exactly 10,000
# attempted repetitions are present, so it cannot silently read the old N=2,000 run.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

source("01_SourceAttack_Kgeneral_functions.R")

K_GRID <- c(2L, 3L, 5L, 10L, 20L, 50L, 100L)
EXPECTED_N <- 10000L
INPUT_DIR <- "results_formal_Kgrid_N10000"
OUTPUT_DIR <- file.path(INPUT_DIR, "teacher_summary")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

read_one_K <- function(K_now) {
  f <- file.path(INPUT_DIR, sprintf("FORMAL_K%03d_ALL_RESULTS.csv", K_now))
  if (!file.exists(f)) stop("Missing formal result file: ", f)
  df <- utils::read.csv(f, stringsAsFactors = FALSE)
  if (nrow(df) != EXPECTED_N) {
    stop(sprintf("K=%d has %d attempted repetitions; expected exactly %d.",
                 K_now, nrow(df), EXPECTED_N))
  }
  df
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
      K = K_now, failure_class = "none", count = 0L,
      proportion_attempted = 0, stringsAsFactors = FALSE
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

attack <- do.call(rbind, attack_list)
diagnostics <- do.call(rbind, diag_list)
failures <- do.call(rbind, failure_list)
row.names(attack) <- row.names(diagnostics) <- row.names(failures) <- NULL

# ------------------------------------------------------------
# Main teacher-facing wide table: R2-(X,Y) and R3-(X,Y)
# ------------------------------------------------------------
r2 <- attack[attack$condition == "R2_XY", , drop = FALSE]
r3 <- attack[attack$condition == "R3_XY", , drop = FALSE]
r2 <- r2[match(K_GRID, r2$K), , drop = FALSE]
r3 <- r3[match(K_GRID, r3$K), , drop = FALSE]

teacher_main <- data.frame(
  K = K_GRID,
  Random_1_over_K = 1 / K_GRID,
  R2XY_Member = r2$member_accuracy,
  R2XY_Member_minus_Random = r2$member_minus_random,
  R2XY_Fresh = r2$fresh_accuracy,
  R2XY_Member_minus_Fresh = r2$member_minus_fresh,
  R3XY_Member = r3$member_accuracy,
  R3XY_Member_minus_Random = r3$member_minus_random,
  R3XY_Fresh = r3$fresh_accuracy,
  R3XY_Member_minus_Fresh = r3$member_minus_fresh,
  Valid = r2$n_valid,
  stringsAsFactors = FALSE
)

teacher_uncertainty <- data.frame(
  K = K_GRID,
  R2XY_Member_MCSE = r2$member_mcse,
  R2XY_Member_CI_L = r2$member_ci_lower,
  R2XY_Member_CI_U = r2$member_ci_upper,
  R2XY_MemberFresh_MCSE = r2$member_minus_fresh_mcse,
  R2XY_MemberFresh_CI_L = r2$member_minus_fresh_ci_lower,
  R2XY_MemberFresh_CI_U = r2$member_minus_fresh_ci_upper,
  R3XY_Member_MCSE = r3$member_mcse,
  R3XY_Member_CI_L = r3$member_ci_lower,
  R3XY_Member_CI_U = r3$member_ci_upper,
  R3XY_MemberFresh_MCSE = r3$member_minus_fresh_mcse,
  R3XY_MemberFresh_CI_L = r3$member_minus_fresh_ci_lower,
  R3XY_MemberFresh_CI_U = r3$member_minus_fresh_ci_upper,
  stringsAsFactors = FALSE
)

utils::write.csv(teacher_main,
                 file.path(OUTPUT_DIR, "TABLE1_TEACHER_MAIN_R2_R3_XY.csv"),
                 row.names = FALSE)
utils::write.csv(teacher_uncertainty,
                 file.path(OUTPUT_DIR, "TABLE2_R2_R3_XY_UNCERTAINTY.csv"),
                 row.names = FALSE)
utils::write.csv(attack,
                 file.path(OUTPUT_DIR, "TABLE3_ALL_ATTACKS_LONG.csv"),
                 row.names = FALSE)
utils::write.csv(diagnostics,
                 file.path(OUTPUT_DIR, "TABLE4_DIAGNOSTICS_BY_K.csv"),
                 row.names = FALSE)
utils::write.csv(failures,
                 file.path(OUTPUT_DIR, "TABLE5_FAILURES_BY_K.csv"),
                 row.names = FALSE)

# Pretty console table
fmt4 <- function(x) ifelse(is.finite(x), sprintf("%.4f", x), "NA")
pretty <- teacher_main
for (nm in setdiff(names(pretty), c("K", "Valid"))) pretty[[nm]] <- fmt4(pretty[[nm]])

cat("\n=======================================================================\n")
cat("TASK A FORMAL N=10,000 PER K -- TEACHER MAIN TABLE\n")
cat("=======================================================================\n")
print(pretty, row.names = FALSE)
cat("\nInterpretation columns:\n")
cat("  Member_minus_Random = member accuracy - 1/K\n")
cat("  Member_minus_Fresh  = member accuracy - matched-fresh accuracy\n")

# ------------------------------------------------------------
# Plot helpers: equally spaced K positions, actual K shown as labels
# ------------------------------------------------------------
x <- seq_along(K_GRID)
xlab <- as.character(K_GRID)
errbar <- function(x, y, lower, upper, cap = 0.06) {
  ok <- is.finite(x) & is.finite(y) & is.finite(lower) & is.finite(upper)
  if (!any(ok)) return(invisible(NULL))
  arrows(x[ok], lower[ok], x[ok], upper[ok], angle = 90, code = 3,
         length = cap)
}

# Plot 1: main result -- member accuracy vs random baseline
png(file.path(OUTPUT_DIR, "FIG1_MEMBER_ACCURACY_R2_R3_VS_RANDOM.png"),
    width = 1500, height = 950, res = 150)
ytop <- min(1, max(c(r2$member_ci_upper, r3$member_ci_upper, 1 / K_GRID), na.rm = TRUE) * 1.10)
plot(x, r2$member_accuracy, type = "b", pch = 16, lty = 1,
     xaxt = "n", ylim = c(0, ytop), xlab = "Number of sites K",
     ylab = "Member source-site accuracy",
     main = "Task A: Member source-site accuracy as K increases")
axis(1, at = x, labels = xlab)
errbar(x, r2$member_accuracy, r2$member_ci_lower, r2$member_ci_upper)
lines(x, r3$member_accuracy, type = "b", pch = 17, lty = 2)
errbar(x, r3$member_accuracy, r3$member_ci_lower, r3$member_ci_upper)
lines(x, 1 / K_GRID, type = "b", pch = 1, lty = 3)
legend("topright", legend = c("R2-(X,Y) member", "R3-(X,Y) member", "Random 1/K"),
       pch = c(16, 17, 1), lty = c(1, 2, 3), bty = "n")
dev.off()

# Plot 2: member advantage over random, with shifted member CI
png(file.path(OUTPUT_DIR, "FIG2_MEMBER_MINUS_RANDOM_R2_R3.png"),
    width = 1500, height = 950, res = 150)
r2_adv_l <- r2$member_ci_lower - 1 / K_GRID
r2_adv_u <- r2$member_ci_upper - 1 / K_GRID
r3_adv_l <- r3$member_ci_lower - 1 / K_GRID
r3_adv_u <- r3$member_ci_upper - 1 / K_GRID
alim <- max(abs(c(r2_adv_l, r2_adv_u, r3_adv_l, r3_adv_u)), na.rm = TRUE) * 1.15
plot(x, r2$member_minus_random, type = "b", pch = 16, lty = 1,
     xaxt = "n", ylim = c(-alim, alim), xlab = "Number of sites K",
     ylab = "Member accuracy - 1/K",
     main = "Task A: Attack advantage above random guessing")
axis(1, at = x, labels = xlab)
abline(h = 0, lty = 3)
errbar(x, r2$member_minus_random, r2_adv_l, r2_adv_u)
lines(x, r3$member_minus_random, type = "b", pch = 17, lty = 2)
errbar(x, r3$member_minus_random, r3_adv_l, r3_adv_u)
legend("topright", legend = c("R2-(X,Y)", "R3-(X,Y)"),
       pch = c(16,17), lty = c(1,2), bty = "n")
dev.off()

# Plot 3: paired member-fresh difference, with paired 95% CI
png(file.path(OUTPUT_DIR, "FIG3_MEMBER_MINUS_FRESH_R2_R3.png"),
    width = 1500, height = 950, res = 150)
dlim <- max(abs(c(r2$member_minus_fresh_ci_lower, r2$member_minus_fresh_ci_upper,
                  r3$member_minus_fresh_ci_lower, r3$member_minus_fresh_ci_upper)),
            na.rm = TRUE) * 1.15
plot(x, r2$member_minus_fresh, type = "b", pch = 16, lty = 1,
     xaxt = "n", ylim = c(-dlim, dlim), xlab = "Number of sites K",
     ylab = "Member accuracy - matched-fresh accuracy",
     main = "Task A: Training-related source-site signal")
axis(1, at = x, labels = xlab)
abline(h = 0, lty = 3)
errbar(x, r2$member_minus_fresh,
       r2$member_minus_fresh_ci_lower, r2$member_minus_fresh_ci_upper)
lines(x, r3$member_minus_fresh, type = "b", pch = 17, lty = 2)
errbar(x, r3$member_minus_fresh,
       r3$member_minus_fresh_ci_lower, r3$member_minus_fresh_ci_upper)
legend("topright", legend = c("R2-(X,Y)", "R3-(X,Y)"),
       pch = c(16,17), lty = c(1,2), bty = "n")
dev.off()

# Plot 4: R2-(X,Y) member vs fresh vs random
png(file.path(OUTPUT_DIR, "FIG4_R2XY_MEMBER_FRESH_RANDOM.png"),
    width = 1500, height = 950, res = 150)
ytop2 <- min(1, max(c(r2$member_ci_upper, r2$fresh_ci_upper, 1 / K_GRID), na.rm = TRUE) * 1.10)
plot(x, r2$member_accuracy, type = "b", pch = 16, lty = 1,
     xaxt = "n", ylim = c(0, ytop2), xlab = "Number of sites K",
     ylab = "Source-site accuracy",
     main = "Task A: R2-(X,Y), member vs matched fresh")
axis(1, at = x, labels = xlab)
lines(x, r2$fresh_accuracy, type = "b", pch = 1, lty = 2)
lines(x, 1 / K_GRID, type = "b", pch = 2, lty = 3)
legend("topright", legend = c("Member", "Matched fresh", "Random 1/K"),
       pch = c(16,1,2), lty = c(1,2,3), bty = "n")
dev.off()

# Plot 5: external-support growth
png(file.path(OUTPUT_DIR, "FIG5_EXTERNAL_SUPPORT_SIZE.png"),
    width = 1500, height = 950, res = 150)
plot(x, diagnostics$mean_external_support_size, type = "b", pch = 16,
     xaxt = "n", ylim = c(0, max(50, diagnostics$mean_external_support_size, na.rm = TRUE)),
     xlab = "Number of sites K", ylab = "Mean external-support size",
     main = "Task A diagnostic: External-support growth")
axis(1, at = x, labels = xlab)
abline(h = 50, lty = 3)
dev.off()

# Plot 6: Hessian stabilization rate
png(file.path(OUTPUT_DIR, "FIG6_HESSIAN_REGULARIZATION_RATE.png"),
    width = 1500, height = 950, res = 150)
plot(x, diagnostics$mean_hessian_regularization_rate, type = "b", pch = 16,
     xaxt = "n", ylim = c(0, 1), xlab = "Number of sites K",
     ylab = "Mean Hessian regularization rate",
     main = "Task A diagnostic: Refined-estimation stabilization")
axis(1, at = x, labels = xlab)
dev.off()

# Text summary for quick teacher display
sink(file.path(OUTPUT_DIR, "TEACHER_README_RESULT_SUMMARY.txt"))
cat("Task A: Source-site attack sensitivity to the number of sites\n")
cat("Formal design: N=10,000 attempted repetitions for each K\n")
cat("K grid: ", paste(K_GRID, collapse = ", "), "\n\n", sep = "")
cat("Main quantities to discuss:\n")
cat("1) Member accuracy compared with random baseline 1/K.\n")
cat("2) Member accuracy compared with matched fresh accuracy.\n")
cat("3) R2-(X,Y) versus R3-(X,Y) as K increases.\n")
cat("4) External-support size and Hessian stabilization as structural diagnostics.\n\n")
print(pretty, row.names = FALSE)
sink()

cat("\n=======================================================================\n")
cat("FORMAL N=10,000 SUMMARY COMPLETE\n")
cat("Teacher-facing outputs saved in:\n  ", normalizePath(OUTPUT_DIR), "\n", sep = "")
cat("Start with TABLE1_TEACHER_MAIN_R2_R3_XY.csv and FIG1--FIG3.\n")
cat("=======================================================================\n")
