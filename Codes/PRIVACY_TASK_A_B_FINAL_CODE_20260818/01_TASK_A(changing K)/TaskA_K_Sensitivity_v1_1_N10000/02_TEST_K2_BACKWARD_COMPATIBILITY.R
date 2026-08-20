# ============================================================
# Task A -- K=2 backward-compatibility regression test
# ============================================================
# MUST PASS before running any K>2 experiment.
#
# It compares the new K-general implementation with the exact Frozen
# Protocol v1.1.1 implementation on the same K=2 repetition IDs.
# The goal is not "similar" results: on valid repetitions the attack
# decision quantities should agree to numerical tolerance.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ---------------- USER-ADJUSTABLE TEST SIZE ------------------
N_TEST <- 5L
TOL <- 1e-10
# -------------------------------------------------------------

OLD_FILE <- file.path("reference_old", "00_REFERENCE_FROZEN_K2_functions.R")
NEW_FILE <- "01_SourceAttack_Kgeneral_functions.R"

if (!file.exists(OLD_FILE)) stop("Missing old reference file: ", OLD_FILE)
if (!file.exists(NEW_FILE)) stop("Missing new function file: ", NEW_FILE)

old_env <- new.env(parent = globalenv())
new_env <- new.env(parent = globalenv())
sys.source(OLD_FILE, envir = old_env)
sys.source(NEW_FILE, envir = new_env)

old_env$check_required_packages()
new_env$check_required_packages()

old_cfg <- old_env$make_config(n_reps = N_TEST)
new_cfg <- new_env$make_config(n_reps = N_TEST, K = 2L)

conditions <- c("R1_X", "R1_XY", "R2_X", "R2_XY", "R3_X", "R3_XY")

cat("\n============================================================\n")
cat("Task A K=2 BACKWARD-COMPATIBILITY TEST\n")
cat("Repetitions:", N_TEST, "| tolerance:", format(TOL, scientific = TRUE), "\n")
cat("Old: Frozen Protocol v1.1.1\n")
cat("New: K-General Source Attack v1.0\n")
cat("============================================================\n\n")

rows <- vector("list", N_TEST)
max_score_error <- 0
max_beta_error <- 0
n_valid_compared <- 0L
n_failed_both <- 0L
all_pass <- TRUE

for (rep_id in seq_len(N_TEST)) {
  old_res <- old_env$run_one_rep(rep_id, old_cfg, return_debug = TRUE)
  new_res <- new_env$run_one_rep_K(rep_id, new_cfg, return_debug = TRUE)

  old_status <- old_res$summary$status
  new_status <- new_res$summary$status
  status_match <- identical(old_status, new_status)

  rep_pass <- status_match
  msg <- character(0)

  if (!status_match) {
    msg <- c(msg, paste0("status mismatch old=", old_status, " new=", new_status))
  } else if (identical(old_status, "failed")) {
    n_failed_both <- n_failed_both + 1L
    # Same status is enough for a failed repetition; exact error text can vary.
  } else {
    n_valid_compared <- n_valid_compared + 1L

    od <- old_res$debug
    nd <- new_res$debug

    # Candidate identity must be exactly preserved.
    candidate_match <- identical(as.integer(old_res$summary$source_site), as.integer(new_res$summary$source_site)) &&
      identical(as.integer(old_res$summary$member_index), as.integer(new_res$summary$member_index)) &&
      identical(as.integer(old_res$summary$member_y), as.integer(new_res$summary$member_y)) &&
      identical(as.integer(old_res$summary$fresh_y), as.integer(new_res$summary$fresh_y))

    if (!candidate_match) {
      rep_pass <- FALSE
      msg <- c(msg, "candidate/source RNG mismatch")
    }

    # First-stage supports and refined coefficient matrices must match.
    support_match <- all(vapply(seq_len(2L), function(k) {
      identical(od$fedfit$support_local[[k]], nd$fedfit$support_local[[k]])
    }, logical(1)))

    if (!support_match) {
      rep_pass <- FALSE
      msg <- c(msg, "first-stage support mismatch")
    }

    beta_err <- max(abs(od$fedfit$beta_refined - nd$fedfit$beta_refined), na.rm = TRUE)
    max_beta_error <- max(max_beta_error, beta_err)
    if (!is.finite(beta_err) || beta_err > TOL) {
      rep_pass <- FALSE
      msg <- c(msg, sprintf("refined beta mismatch %.3e", beta_err))
    }

    central_match <- identical(
      as.integer(od$fedfit$central$final_selection),
      as.integer(nd$fedfit$central$final_selection)
    )
    if (!central_match) {
      rep_pass <- FALSE
      msg <- c(msg, "central selection mismatch")
    }

    # Old D-score must equal Q_site1 - Q_site2 in the K-general definition.
    for (prefix in c("member", "fresh")) {
      old_attack <- od[[paste0(prefix, "_attack")]]
      new_attack <- nd[[paste0(prefix, "_attack")]]

      for (cc in conditions) {
        old_D <- old_attack$scores[[cc]]
        new_Q <- new_attack$scores[[cc]]
        new_D <- new_Q[1L] - new_Q[2L]
        err <- abs(old_D - new_D)
        max_score_error <- max(max_score_error, err, na.rm = TRUE)

        if (!is.finite(err) || err > TOL) {
          rep_pass <- FALSE
          msg <- c(msg, sprintf("%s %s D/Q mismatch %.3e", prefix, cc, err))
        }

        old_credit <- old_attack$credits[[cc]]
        new_credit <- new_attack$decisions[[cc]]$credit
        if (!isTRUE(all.equal(old_credit, new_credit, tolerance = TOL))) {
          rep_pass <- FALSE
          msg <- c(msg, sprintf("%s %s credit mismatch", prefix, cc))
        }
      }

      # R0 must reduce exactly to the old 0.5 structural credit.
      for (cc in c("R0_X", "R0_XY")) {
        old_credit <- old_attack$credits[[cc]]
        new_credit <- new_attack$decisions[[cc]]$credit
        if (!isTRUE(all.equal(old_credit, new_credit, tolerance = TOL))) {
          rep_pass <- FALSE
          msg <- c(msg, sprintf("%s %s R0 credit mismatch", prefix, cc))
        }
      }
    }
  }

  all_pass <- all_pass && rep_pass
  rows[[rep_id]] <- data.frame(
    rep_id = rep_id,
    old_status = old_status,
    new_status = new_status,
    pass = rep_pass,
    note = if (length(msg) == 0L) "" else paste(unique(msg), collapse = "; "),
    stringsAsFactors = FALSE
  )

  cat(sprintf("rep %d/%d : %s\n", rep_id, N_TEST, if (rep_pass) "PASS" else "FAIL"))
  if (length(msg) > 0L) cat("  ", paste(unique(msg), collapse = "; "), "\n")
}

report <- do.call(rbind, rows)
utils::write.csv(report, "K2_backward_compatibility_report.csv", row.names = FALSE)

cat("\n============================================================\n")
cat("K=2 REGRESSION TEST SUMMARY\n")
cat("Valid repetitions compared:", n_valid_compared, "\n")
cat("Failed in both implementations:", n_failed_both, "\n")
cat("Maximum refined-beta error:", format(max_beta_error, scientific = TRUE), "\n")
cat("Maximum old-D vs new-Q-difference error:", format(max_score_error, scientific = TRUE), "\n")
cat("Overall:", if (all_pass) "PASS" else "FAIL", "\n")
cat("Report: K2_backward_compatibility_report.csv\n")
cat("============================================================\n")

if (!all_pass) {
  stop("Backward-compatibility test FAILED. Do NOT run the K-grid until this is resolved.")
}

cat("\nPASS: the K-general attack definition reduces to Frozen v1.1.1 at K=2.\n")
cat("Next: run 03_RUN_QUICK_KGRID.R\n")
