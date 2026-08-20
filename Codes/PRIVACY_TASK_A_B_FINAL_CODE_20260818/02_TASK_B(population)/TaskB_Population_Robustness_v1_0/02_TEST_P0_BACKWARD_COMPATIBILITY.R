# ============================================================
# Task B gate 1: P0 backward-compatibility test
# The new Task-B pipeline MUST reproduce Frozen v1.1.1 exactly
# when population = P0_TOEPLITZ_GAUSSIAN.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

source(file.path("reference_old", "00_REFERENCE_FROZEN_K2_functions.R"))
old_make_config <- make_config
old_run_one_rep <- run_one_rep

source("01_TaskB_PopulationAttack_functions.R")
check_required_packages()

N_TEST <- 5L
BASE_SEED <- 20260807L
TOL <- 1e-12

old_cfg <- old_make_config(n_reps = N_TEST, base_seed = BASE_SEED)
new_cfg <- make_taskb_config(
  population = "P0_TOEPLITZ_GAUSSIAN",
  n_reps = N_TEST,
  base_seed = BASE_SEED,
  diagnostic_n = 0L
)

cat("\n============================================================\n")
cat("TASK B P0 BACKWARD-COMPATIBILITY TEST\n")
cat("Frozen old K=2 pipeline vs new Task-B P0 pipeline\n")
cat("============================================================\n\n")

all_pass <- TRUE
max_numeric_error <- 0

for (r in seq_len(N_TEST)) {
  old <- old_run_one_rep(r, old_cfg, return_debug = TRUE)
  new <- run_one_rep_taskb(r, new_cfg, return_debug = TRUE)

  pass <- TRUE
  notes <- character(0)

  if (!identical(old$summary$status, new$summary$status)) {
    pass <- FALSE
    notes <- c(notes, "status mismatch")
  }

  if (identical(old$summary$status, "ok") && identical(new$summary$status, "ok")) {
    # Exact data/candidate checks.
    for (k in 1:2) {
      if (!isTRUE(all.equal(old$debug$site_data[[k]]$X, new$debug$site_data[[k]]$X, tolerance = 0))) {
        pass <- FALSE; notes <- c(notes, paste0("site", k, " X mismatch"))
      }
      if (!identical(old$debug$site_data[[k]]$y, new$debug$site_data[[k]]$y)) {
        pass <- FALSE; notes <- c(notes, paste0("site", k, " y mismatch"))
      }
      if (!identical(old$debug$fedfit$support_local[[k]], new$debug$fedfit$support_local[[k]])) {
        pass <- FALSE; notes <- c(notes, paste0("site", k, " first-stage support mismatch"))
      }
    }

    beta_err <- max(abs(old$debug$fedfit$beta_refined - new$debug$fedfit$beta_refined))
    if (!is.finite(beta_err)) beta_err <- Inf
    max_numeric_error <- max(max_numeric_error, beta_err)
    if (beta_err > TOL) {
      pass <- FALSE; notes <- c(notes, sprintf("refined beta error %.3e", beta_err))
    }

    if (!identical(old$debug$member$y, new$debug$member$y) ||
        !isTRUE(all.equal(old$debug$member$x, new$debug$member$x, tolerance = 0))) {
      pass <- FALSE; notes <- c(notes, "member candidate mismatch")
    }
    if (!identical(old$debug$fresh$y, new$debug$fresh$y) ||
        !isTRUE(all.equal(old$debug$fresh$x, new$debug$fresh$x, tolerance = 0))) {
      pass <- FALSE; notes <- c(notes, "fresh candidate mismatch")
    }

    conds <- c("R0_X", "R0_XY", "R1_X", "R1_XY", "R2_X", "R2_XY", "R3_X", "R3_XY")
    for (cc in conds) {
      om <- old$debug$member_attack$credits[[cc]]
      nm <- new$debug$member_attack$credits[[cc]]
      of <- old$debug$fresh_attack$credits[[cc]]
      nf <- new$debug$fresh_attack$credits[[cc]]
      if (!isTRUE(all.equal(om, nm, tolerance = TOL)) || !isTRUE(all.equal(of, nf, tolerance = TOL))) {
        pass <- FALSE; notes <- c(notes, paste0(cc, " credit mismatch"))
      }
    }
  } else if (identical(old$summary$status, "failed") && identical(new$summary$status, "failed")) {
    if (!identical(old$summary$failure_stage, new$summary$failure_stage) ||
        !identical(old$summary$failure_reason, new$summary$failure_reason)) {
      pass <- FALSE; notes <- c(notes, "failure reason mismatch")
    }
  }

  all_pass <- all_pass && pass
  cat(sprintf("rep %d/%d : %s", r, N_TEST, if (pass) "PASS" else "FAIL"))
  if (length(notes) > 0L) cat(" | ", paste(unique(notes), collapse = "; "), sep = "")
  cat("\n")
}

cat("\nMaximum refined-beta error:", format(max_numeric_error, scientific = TRUE), "\n")
cat("Overall:", if (all_pass) "PASS" else "FAIL", "\n")

if (!all_pass) {
  stop("Task-B P0 compatibility test FAILED. Do not run Quick/Pilot until fixed.")
}

cat("\nPASS: Task-B P0 is exactly backward-compatible with Frozen v1.1.1.\n")
cat("Next: run 03_RUN_QUICK_MAIN_POPULATIONS.R\n")
