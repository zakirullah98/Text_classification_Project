# ============================================================
# Task B formal summary + teacher-ready tables/figures
# Reads ONLY N=10,000 formal MAIN-population outputs for attack results.
# If Pilot DGP diagnostics are available, they are copied separately as
# diagnostic context; they are not mixed into the formal attack estimands.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)
source("01_TaskB_PopulationAttack_functions.R")

INPUT_DIR <- "results_taskB_formal_N10000"
PILOT_DIR <- "results_taskB_pilot"
OUT_DIR <- file.path(INPUT_DIR, "teacher_summary")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

all_results <- list()
attack_summaries <- list()
status_rows <- list()

for (pop in TASKB_MAIN_POPULATIONS) {
  f <- file.path(INPUT_DIR, pop, "FORMAL_N10000_ALL_RESULTS.rds")
  if (!file.exists(f)) stop("Missing formal result: ", f)
  df <- readRDS(f)
  if (nrow(df) != 10000L || !all(as.integer(df$rep_id) == 1:10000)) {
    stop("Formal N=10,000 audit failed for ", pop)
  }
  all_results[[pop]] <- df
  attack_summaries[[pop]] <- summarize_taskb_attacks(df)
  status_rows[[pop]] <- data.frame(
    population = pop,
    population_label = population_label(pop),
    attempted = 10000L,
    valid = sum(df$status == "ok", na.rm = TRUE),
    failed = sum(df$status == "failed", na.rm = TRUE),
    failure_rate = mean(df$status == "failed", na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

attack_all <- do.call(rbind, attack_summaries)
status_all <- do.call(rbind, status_rows)

utils::write.csv(status_all, file.path(OUT_DIR, "TABLE0_RUN_STATUS.csv"), row.names = FALSE)
utils::write.csv(attack_all, file.path(OUT_DIR, "TABLE_FULL_ALL_REGIMES.csv"), row.names = FALSE)

# Teacher main table: response-aware coefficient and full-information attacks.
main <- attack_all[attack_all$condition %in% c("R2_XY", "R3_XY"), , drop = FALSE]
main$random_baseline <- 0.5
main$member_minus_random <- main$member_accuracy - 0.5
main <- main[, c(
  "population", "population_label", "condition",
  "member_accuracy", "random_baseline", "member_minus_random",
  "member_ci_lower", "member_ci_upper",
  "fresh_accuracy", "fresh_ci_lower", "fresh_ci_upper",
  "delta_train", "delta_train_ci_lower", "delta_train_ci_upper",
  "n_valid"
)]
utils::write.csv(main, file.path(OUT_DIR, "TABLE1_TEACHER_MAIN_R2_R3_XY.csv"), row.names = FALSE)

# Wide version for a slide/table.
r2 <- main[main$condition == "R2_XY", ]
r3 <- main[main$condition == "R3_XY", ]
wide <- data.frame(
  population = r2$population,
  population_label = r2$population_label,
  random = 0.5,
  R2_XY_member = r2$member_accuracy,
  R2_XY_member_minus_random = r2$member_minus_random,
  R2_XY_fresh = r2$fresh_accuracy,
  R2_XY_member_minus_fresh = r2$delta_train,
  R3_XY_member = r3$member_accuracy[match(r2$population, r3$population)],
  R3_XY_member_minus_random = r3$member_minus_random[match(r2$population, r3$population)],
  R3_XY_fresh = r3$fresh_accuracy[match(r2$population, r3$population)],
  R3_XY_member_minus_fresh = r3$delta_train[match(r2$population, r3$population)],
  stringsAsFactors = FALSE
)
utils::write.csv(wide, file.path(OUT_DIR, "TABLE2_TEACHER_WIDE.csv"), row.names = FALSE)

# Full regime member/fresh table.
full_teacher <- attack_all[, c(
  "population", "population_label", "condition",
  "member_accuracy", "attack_advantage", "fresh_accuracy", "delta_train",
  "member_ci_lower", "member_ci_upper", "delta_train_ci_lower", "delta_train_ci_upper", "n_valid"
)]
utils::write.csv(full_teacher, file.path(OUT_DIR, "TABLE3_ALL_ATTACKS_TEACHER.csv"), row.names = FALSE)

# Optional Pilot DGP diagnostics kept separate.
pilot_diag_file <- file.path(PILOT_DIR, "PILOT_DGP_DIAGNOSTICS.csv")
if (file.exists(pilot_diag_file)) {
  pilot_diag <- utils::read.csv(pilot_diag_file, stringsAsFactors = FALSE)
  utils::write.csv(pilot_diag, file.path(OUT_DIR, "TABLE4_PILOT_DGP_DIAGNOSTICS_CONTEXT.csv"), row.names = FALSE)
}

# ----------------------- Figures -----------------------------
pop_order <- TASKB_MAIN_POPULATIONS
pop_short <- c("P0\nGaussian", "P1\nCopula Uniform", "P2\nBlock sparse", "P3\nNear-singular", "P4\nPoisson")

get_ordered <- function(cond) {
  z <- main[main$condition == cond, , drop = FALSE]
  z[match(pop_order, z$population), , drop = FALSE]
}

r2o <- get_ordered("R2_XY")
r3o <- get_ordered("R3_XY")

png(file.path(OUT_DIR, "FIG1_MEMBER_ACCURACY_R2_R3_VS_RANDOM.png"), width = 1500, height = 900, res = 160)
x <- seq_along(pop_order)
ylim <- range(c(r2o$member_ci_lower, r2o$member_ci_upper, r3o$member_ci_lower, r3o$member_ci_upper, 0.5), na.rm = TRUE)
ylim <- c(max(0, ylim[1] - 0.03), min(1, ylim[2] + 0.03))
plot(x, r2o$member_accuracy, type = "b", pch = 16, xaxt = "n", ylim = ylim,
     xlab = "Population", ylab = "Member source-site accuracy",
     main = "Task B: Member source-site accuracy across populations")
axis(1, at = x, labels = pop_short)
segments(x, r2o$member_ci_lower, x, r2o$member_ci_upper)
lines(x, r3o$member_accuracy, type = "b", pch = 17)
segments(x + 0.04, r3o$member_ci_lower, x + 0.04, r3o$member_ci_upper)
abline(h = 0.5, lty = 2)
legend("topright", legend = c("R2-(X,Y)", "R3-(X,Y)", "Random = 0.5"),
       pch = c(16, 17, NA), lty = c(1, 1, 2), bty = "n")
dev.off()

png(file.path(OUT_DIR, "FIG2_MEMBER_MINUS_RANDOM_R2_R3.png"), width = 1500, height = 900, res = 160)
r2_adv_lo <- r2o$member_ci_lower - 0.5
r2_adv_hi <- r2o$member_ci_upper - 0.5
r3_adv_lo <- r3o$member_ci_lower - 0.5
r3_adv_hi <- r3o$member_ci_upper - 0.5
ylim <- range(c(r2_adv_lo, r2_adv_hi, r3_adv_lo, r3_adv_hi, 0), na.rm = TRUE)
plot(x, r2o$member_minus_random, type = "b", pch = 16, xaxt = "n", ylim = ylim,
     xlab = "Population", ylab = "Member accuracy - 0.5",
     main = "Task B: Attack advantage above random guessing")
axis(1, at = x, labels = pop_short)
segments(x, r2_adv_lo, x, r2_adv_hi)
lines(x, r3o$member_minus_random, type = "b", pch = 17)
segments(x + 0.04, r3_adv_lo, x + 0.04, r3_adv_hi)
abline(h = 0, lty = 2)
legend("topright", legend = c("R2-(X,Y)", "R3-(X,Y)"), pch = c(16,17), lty = 1, bty = "n")
dev.off()

png(file.path(OUT_DIR, "FIG3_MEMBER_MINUS_FRESH_R2_R3.png"), width = 1500, height = 900, res = 160)
ylim <- range(c(r2o$delta_train_ci_lower, r2o$delta_train_ci_upper,
                r3o$delta_train_ci_lower, r3o$delta_train_ci_upper, 0), na.rm = TRUE)
plot(x, r2o$delta_train, type = "b", pch = 16, xaxt = "n", ylim = ylim,
     xlab = "Population", ylab = "Member accuracy - matched-fresh accuracy",
     main = "Task B: Training-related source-site signal")
axis(1, at = x, labels = pop_short)
segments(x, r2o$delta_train_ci_lower, x, r2o$delta_train_ci_upper)
lines(x, r3o$delta_train, type = "b", pch = 17)
segments(x + 0.04, r3o$delta_train_ci_lower, x + 0.04, r3o$delta_train_ci_upper)
abline(h = 0, lty = 2)
legend("topright", legend = c("R2-(X,Y)", "R3-(X,Y)"), pch = c(16,17), lty = 1, bty = "n")
dev.off()

png(file.path(OUT_DIR, "FIG4_R2XY_MEMBER_FRESH_RANDOM.png"), width = 1500, height = 900, res = 160)
ylim <- range(c(r2o$member_ci_lower, r2o$member_ci_upper, r2o$fresh_ci_lower, r2o$fresh_ci_upper, 0.5), na.rm = TRUE)
plot(x, r2o$member_accuracy, type = "b", pch = 16, xaxt = "n", ylim = ylim,
     xlab = "Population", ylab = "Accuracy",
     main = "R2-(X,Y): member vs matched fresh across populations")
axis(1, at = x, labels = pop_short)
segments(x, r2o$member_ci_lower, x, r2o$member_ci_upper)
lines(x, r2o$fresh_accuracy, type = "b", pch = 1)
segments(x + 0.04, r2o$fresh_ci_lower, x + 0.04, r2o$fresh_ci_upper)
abline(h = 0.5, lty = 2)
legend("topright", legend = c("Member", "Matched fresh", "Random = 0.5"),
       pch = c(16, 1, NA), lty = c(1,1,2), bty = "n")
dev.off()

cat("\n============================================================\n")
cat("TASK B SUMMARY COMPLETE\n")
cat("Teacher-ready outputs:", normalizePath(OUT_DIR), "\n")
cat("Main table: TABLE2_TEACHER_WIDE.csv\n")
cat("Main figures: FIG1--FIG4\n")
cat("============================================================\n")
