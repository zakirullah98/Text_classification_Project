# ============================================================
# Task B QUICK: all five MAIN populations
# Purpose: code-path validation + DGP sanity diagnostics.
# Attack accuracies here are NOT formal results.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)
source("01_TaskB_PopulationAttack_functions.R")
check_required_packages()

N_QUICK <- 5L
DIAGNOSTIC_N <- 1000L
BASE_SEED <- 20260815L
OUTPUT_DIR <- "results_taskB_quick"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("\n============================================================\n")
cat("TASK B QUICK -- MAIN POPULATIONS\n")
cat("K=2 | n/site=50 | p=50 | N/population=", N_QUICK,
    " | oracle diagnostic M=", DIAGNOSTIC_N, "\n", sep = "")
cat("ATTACK FIXED; ONLY X POPULATION CHANGES\n")
cat("============================================================\n")

all_results <- list()
all_dgp <- list()
all_attack <- list()

for (pop in TASKB_MAIN_POPULATIONS) {
  cfg <- make_taskb_config(pop, n_reps = N_QUICK, base_seed = BASE_SEED,
                           diagnostic_n = DIAGNOSTIC_N)
  cat("\n------------------------------------------------------------\n")
  cat(cfg$population_label, "\n")
  cat("------------------------------------------------------------\n")

  rows <- vector("list", N_QUICK)
  tpop <- Sys.time()
  for (r in seq_len(N_QUICK)) {
    t0 <- Sys.time()
    rows[[r]] <- run_one_rep_taskb(r, cfg, return_debug = FALSE)
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("rep %d/%d | %-6s | %.2f sec\n",
                r, N_QUICK, rows[[r]]$status %||% "?", dt))
  }

  df <- lists_to_data_frame(rows)
  all_results[[pop]] <- df
  dgp <- summarize_taskb_dgp(df)
  atk <- summarize_taskb_attacks(df)
  all_dgp[[pop]] <- dgp
  all_attack[[pop]] <- atk

  pdir <- file.path(OUTPUT_DIR, pop)
  dir.create(pdir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(df, file.path(pdir, "quick_results.rds"))
  utils::write.csv(df, file.path(pdir, "quick_results.csv"), row.names = FALSE)
  saveRDS(cfg, file.path(pdir, "config.rds"))
  write_taskb_p2_permutation(cfg, file.path(OUTPUT_DIR, "P2_FIXED_PERMUTATION.csv"))

  ok <- df[df$status == "ok", , drop = FALSE]
  cat("Valid:", nrow(ok), " Failed:", sum(df$status == "failed", na.rm = TRUE), "\n")
  if (nrow(dgp) > 0L) {
    cat("DGP diagnostics (means over valid reps):\n")
    print(dgp[, c(
      "train_ybar_site1", "train_ybar_site2",
      "oracle_sd_xbeta_site1", "oracle_sd_xbeta_site2",
      "oracle_extreme_pi_site1", "oracle_extreme_pi_site2",
      "oracle_auc_site1", "oracle_auc_site2",
      "cov_lambda_min_site1", "cov_kappa_site1", "cov_det_root_site1"
    )], row.names = FALSE)
  }
  if (nrow(atk) > 0L) {
    cat("Quick attack check (DO NOT INTERPRET FORMALLY):\n")
    print(atk[atk$condition %in% c("R2_XY", "R3_XY"),
              c("condition", "member_accuracy", "fresh_accuracy", "delta_train")], row.names = FALSE)
  }
  cat("Population runtime:", format_seconds(as.numeric(difftime(Sys.time(), tpop, units = "secs"))), "\n")
}

combined <- do.call(rbind, all_results)
dgp_table <- do.call(rbind, all_dgp)
attack_table <- do.call(rbind, all_attack)

saveRDS(combined, file.path(OUTPUT_DIR, "QUICK_ALL_POPULATIONS.rds"))
utils::write.csv(combined, file.path(OUTPUT_DIR, "QUICK_ALL_POPULATIONS.csv"), row.names = FALSE)
utils::write.csv(dgp_table, file.path(OUTPUT_DIR, "QUICK_DGP_DIAGNOSTICS.csv"), row.names = FALSE)
utils::write.csv(attack_table, file.path(OUTPUT_DIR, "QUICK_ATTACK_SUMMARY.csv"), row.names = FALSE)

cat("\n============================================================\n")
cat("QUICK COMPLETE\n")
cat("Inspect DGP diagnostics before Pilot.\n")
cat("Next if sane: 04_RUN_PILOT_MAIN_POPULATIONS.R\n")
cat("============================================================\n")
