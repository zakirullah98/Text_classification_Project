# ============================================================
# Task A QUICK K-grid run
# Purpose: code/shape/runtime validation ONLY -- not inference.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

source("01_SourceAttack_Kgeneral_functions.R")
check_required_packages()

# ============================================================
# USER-ADJUSTABLE COMPUTING SETTINGS
# ============================================================
K_GRID <- c(2L, 3L, 5L, 10L, 20L, 50L, 100L)
N_QUICK_BY_K <- c(
  `2` = 3L,
  `3` = 3L,
  `5` = 3L,
  `10` = 2L,
  `20` = 2L,
  `50` = 1L,
  `100` = 1L
)
N_WORKERS <- 8L
OUTPUT_DIR <- "results_quick_Kgrid"
# ============================================================

if (!file.exists("K2_backward_compatibility_report.csv")) {
  stop(
    "Run 02_TEST_K2_BACKWARD_COMPATIBILITY.R first. ",
    "The K-grid is intentionally blocked until the K=2 regression test exists."
  )
}

backward <- utils::read.csv("K2_backward_compatibility_report.csv", stringsAsFactors = FALSE)
if (!all(backward$pass)) {
  stop("K=2 backward-compatibility report contains failures. Do NOT run K-grid.")
}

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
common_file <- normalizePath("01_SourceAttack_Kgeneral_functions.R")

cat("\n============================================================\n")
cat("Task A QUICK K-GRID\n")
cat("K grid:", paste(K_GRID, collapse = ", "), "\n")
cat("n_k=50 | p=50 | baseline DGP unchanged\n")
cat("Purpose: code validation and rough runtime only.\n")
cat("Workers:", N_WORKERS, "\n")
cat("============================================================\n\n")

all_timing <- list()

for (K_now in K_GRID) {
  N_now <- as.integer(N_QUICK_BY_K[as.character(K_now)])
  if (!is.finite(N_now) || N_now < 1L) stop("Missing quick N for K=", K_now)

  cfg <- make_config(n_reps = N_now, K = K_now)
  ids <- seq_len(N_now)

  cat("\n------------------------------------------------------------\n")
  cat("QUICK K =", K_now, "| N =", N_now, "| random baseline =", sprintf("%.4f", 1 / K_now), "\n")
  cat("------------------------------------------------------------\n")

  t0 <- Sys.time()

  if (N_WORKERS > 1L && N_now > 1L) {
    n_use <- min(N_WORKERS, N_now)
    cl <- parallel::makeCluster(n_use)
    parallel::clusterEvalQ(cl, {
      library(glmnet)
      library(mvtnorm)
      NULL
    })
    parallel::clusterCall(cl, function(path) {
      source(path, local = .GlobalEnv)
      NULL
    }, common_file)

    rows <- parallel::parLapplyLB(
      cl,
      ids,
      function(i, cfg_local) run_one_rep_K(i, cfg_local, return_debug = FALSE),
      cfg_local = cfg
    )
    parallel::stopCluster(cl)
    cl <- NULL
  } else {
    rows <- lapply(ids, function(i) run_one_rep_K(i, cfg, return_debug = FALSE))
  }

  df <- lists_to_data_frame(rows)
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  out_csv <- file.path(OUTPUT_DIR, sprintf("QUICK_K%03d_RESULTS.csv", K_now))
  out_rds <- file.path(OUTPUT_DIR, sprintf("QUICK_K%03d_RESULTS.rds", K_now))
  utils::write.csv(df, out_csv, row.names = FALSE)
  saveRDS(df, out_rds)
  saveRDS(cfg, file.path(OUTPUT_DIR, sprintf("QUICK_K%03d_CONFIG.rds", K_now)))

  valid <- df[df$status == "ok", , drop = FALSE]
  n_valid <- nrow(valid)
  n_fail <- nrow(df) - n_valid

  cat("Runtime:", format_seconds(sec), "\n")
  cat("Valid:", n_valid, "| Failed:", n_fail, "\n")

  if (n_valid > 0L) {
    cat("Mean external support size:", round(mean(valid$mean_external_support_size, na.rm = TRUE), 2), "\n")
    cat("Mean Hessian regularization rate:", round(mean(valid$prop_hessian_regularized, na.rm = TRUE), 4), "\n")
    cat("R2-(X,Y) member accuracy [QUICK ONLY]:",
        round(mean(valid$member_credit_R2_XY, na.rm = TRUE), 4), "\n")
  }

  if (n_fail > 0L) {
    cat("Failure classes:\n")
    print(sort(table(df$failure_class[df$status == "failed"]), decreasing = TRUE))
  }

  all_timing[[length(all_timing) + 1L]] <- data.frame(
    K = K_now,
    N = N_now,
    runtime_seconds = sec,
    seconds_per_attempt = sec / N_now,
    valid = n_valid,
    failed = n_fail,
    failure_rate = n_fail / N_now,
    stringsAsFactors = FALSE
  )
}

timing <- do.call(rbind, all_timing)
utils::write.csv(timing, file.path(OUTPUT_DIR, "QUICK_RUNTIME_BY_K.csv"), row.names = FALSE)

cat("\n============================================================\n")
cat("QUICK K-GRID COMPLETE\n")
cat("Results:", normalizePath(OUTPUT_DIR), "\n")
cat("IMPORTANT: these tiny-N accuracies are NOT scientific conclusions.\n")
cat("Inspect failures/runtime, then run 04_RUN_PILOT_KGRID.R.\n")
cat("============================================================\n")
