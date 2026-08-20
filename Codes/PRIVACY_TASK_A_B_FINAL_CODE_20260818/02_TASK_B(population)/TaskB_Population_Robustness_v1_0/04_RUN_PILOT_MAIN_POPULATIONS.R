# ============================================================
# Task B PILOT: five MAIN populations
# Purpose: verify DGP is not degenerate and estimate runtime/failure rates.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)
source("01_TaskB_PopulationAttack_functions.R")
check_required_packages()

N_PILOT <- 100L
DIAGNOSTIC_N <- 2000L
BASE_SEED <- 20260815L
N_WORKERS <- 8L
OUTPUT_DIR <- "results_taskB_pilot"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
common_file <- normalizePath("01_TaskB_PopulationAttack_functions.R")

cat("\n============================================================\n")
cat("TASK B PILOT -- MAIN POPULATIONS\n")
cat("K=2 | n/site=50 | p=50 | N/population=", N_PILOT,
    " | oracle diagnostic M=", DIAGNOSTIC_N,
    " | workers=", N_WORKERS, "\n", sep = "")
cat("============================================================\n")

cl <- NULL
if (N_WORKERS > 1L) {
  cl <- parallel::makeCluster(N_WORKERS)
  parallel::clusterEvalQ(cl, { library(glmnet); library(mvtnorm); NULL })
  parallel::clusterCall(cl, function(path) { source(path, local = .GlobalEnv); NULL }, common_file)
}

all_results <- list()
all_dgp <- list()
all_attack <- list()
runtime_rows <- list()

for (pop in TASKB_MAIN_POPULATIONS) {
  cfg <- make_taskb_config(pop, n_reps = N_PILOT, base_seed = BASE_SEED,
                           diagnostic_n = DIAGNOSTIC_N)
  cat("\n------------------------------------------------------------\n")
  cat(cfg$population_label, "\n")
  cat("------------------------------------------------------------\n")
  t0 <- Sys.time()

  ids <- seq_len(N_PILOT)
  if (is.null(cl)) {
    rows <- lapply(ids, function(i) run_one_rep_taskb(i, cfg, FALSE))
  } else {
    rows <- parallel::parLapplyLB(
      cl, ids,
      function(i, cfg_local) run_one_rep_taskb(i, cfg_local, FALSE),
      cfg_local = cfg
    )
  }

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  df <- lists_to_data_frame(rows)
  all_results[[pop]] <- df
  dgp <- summarize_taskb_dgp(df)
  atk <- summarize_taskb_attacks(df)
  all_dgp[[pop]] <- dgp
  all_attack[[pop]] <- atk

  valid_n <- sum(df$status == "ok", na.rm = TRUE)
  fail_n <- sum(df$status == "failed", na.rm = TRUE)
  runtime_rows[[pop]] <- data.frame(
    population = pop,
    population_label = cfg$population_label,
    attempted = N_PILOT,
    valid = valid_n,
    failed = fail_n,
    failure_rate = fail_n / N_PILOT,
    elapsed_seconds = elapsed,
    seconds_per_attempt = elapsed / N_PILOT,
    stringsAsFactors = FALSE
  )

  pdir <- file.path(OUTPUT_DIR, pop)
  dir.create(pdir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(df, file.path(pdir, "pilot_results.rds"))
  utils::write.csv(df, file.path(pdir, "pilot_results.csv"), row.names = FALSE)
  saveRDS(cfg, file.path(pdir, "config.rds"))

  cat("Attempted:", N_PILOT, " Valid:", valid_n, " Failed:", fail_n,
      " Runtime:", format_seconds(elapsed), "\n")
  if (fail_n > 0L) {
    cat("Failure stages:\n")
    print(sort(table(df$failure_stage[df$status == "failed"]), decreasing = TRUE))
  }
  if (nrow(dgp) > 0L) {
    cat("DGP sanity diagnostics:\n")
    print(dgp[, c(
      "train_ybar_site1", "train_ybar_site2",
      "oracle_sd_xbeta_site1", "oracle_sd_xbeta_site2",
      "oracle_extreme_pi_site1", "oracle_extreme_pi_site2",
      "oracle_mean_pi_site1", "oracle_mean_pi_site2",
      "oracle_auc_site1", "oracle_auc_site2",
      "cov_lambda_min_site1", "cov_kappa_site1", "cov_det_root_site1"
    )], row.names = FALSE)
  }
  cat("Pilot attack summary (trend only):\n")
  print(atk[atk$condition %in% c("R2_XY", "R3_XY"),
            c("condition", "member_accuracy", "fresh_accuracy", "delta_train")], row.names = FALSE)
}

if (!is.null(cl)) parallel::stopCluster(cl)

combined <- do.call(rbind, all_results)
dgp_table <- do.call(rbind, all_dgp)
attack_table <- do.call(rbind, all_attack)
runtime_table <- do.call(rbind, runtime_rows)

saveRDS(combined, file.path(OUTPUT_DIR, "PILOT_ALL_POPULATIONS.rds"))
utils::write.csv(combined, file.path(OUTPUT_DIR, "PILOT_ALL_POPULATIONS.csv"), row.names = FALSE)
utils::write.csv(dgp_table, file.path(OUTPUT_DIR, "PILOT_DGP_DIAGNOSTICS.csv"), row.names = FALSE)
utils::write.csv(attack_table, file.path(OUTPUT_DIR, "PILOT_ATTACK_SUMMARY.csv"), row.names = FALSE)
utils::write.csv(runtime_table, file.path(OUTPUT_DIR, "PILOT_RUNTIME_FAILURES.csv"), row.names = FALSE)

cfg0 <- make_taskb_config(TASKB_MAIN_POPULATIONS[1], n_reps = N_PILOT)
write_taskb_p2_permutation(cfg0, file.path(OUTPUT_DIR, "P2_FIXED_PERMUTATION.csv"))

cat("\n============================================================\n")
cat("PILOT COMPLETE\n")
cat("Do NOT tune parameters based on attack accuracy.\n")
cat("Only stop/revise if DGP diagnostics reveal genuine degeneration or invalid covariance.\n")
cat("If sane, next: 05_RUN_FORMAL_MAIN_POPULATIONS_N10000.R\n")
cat("============================================================\n")
