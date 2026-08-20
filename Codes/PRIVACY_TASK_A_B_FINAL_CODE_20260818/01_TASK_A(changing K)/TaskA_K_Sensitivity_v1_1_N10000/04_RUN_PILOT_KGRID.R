# ============================================================
# Task A PILOT K-grid run
# Purpose: runtime, numerical stability, failure-rate diagnostics.
# This is still not the final formal experiment.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

source("01_SourceAttack_Kgeneral_functions.R")
check_required_packages()

# ============================================================
# USER-ADJUSTABLE COMPUTING SETTINGS
# ============================================================
K_GRID <- c(2L, 3L, 5L, 10L, 20L, 50L, 100L)
N_PILOT <- 20L
N_WORKERS <- 8L
BATCH_SIZE <- 5L
OUTPUT_DIR <- "results_pilot_Kgrid"
# ============================================================

if (!file.exists("K2_backward_compatibility_report.csv")) {
  stop("Run 02_TEST_K2_BACKWARD_COMPATIBILITY.R first.")
}

backward <- utils::read.csv("K2_backward_compatibility_report.csv", stringsAsFactors = FALSE)
if (!all(backward$pass)) stop("K=2 backward-compatibility test did not pass.")

if (!dir.exists("results_quick_Kgrid")) {
  stop("Run 03_RUN_QUICK_KGRID.R first and inspect the quick outputs.")
}

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
common_file <- normalizePath("01_SourceAttack_Kgeneral_functions.R")

quick_runtime_file <- file.path("results_quick_Kgrid", "QUICK_RUNTIME_BY_K.csv")
if (file.exists(quick_runtime_file)) {
  quick_timing <- utils::read.csv(quick_runtime_file, stringsAsFactors = FALSE)
  cat("\nRough pilot runtime projection from QUICK run (very approximate):\n")
  for (K_now in K_GRID) {
    q <- quick_timing[quick_timing$K == K_now, , drop = FALSE]
    if (nrow(q) == 1L && is.finite(q$seconds_per_attempt)) {
      # Rough sequential estimate; actual parallel runtime should be lower.
      est_seq <- q$seconds_per_attempt * N_PILOT
      cat(sprintf("  K=%3d : sequential-equivalent ~ %s before parallel speedup\n",
                  K_now, format_seconds(est_seq)))
    }
  }
}

cat("\n============================================================\n")
cat("Task A PILOT K-GRID\n")
cat("K grid:", paste(K_GRID, collapse = ", "), "\n")
cat("N pilot per K:", N_PILOT, "\n")
cat("workers:", N_WORKERS, "| batch size:", BATCH_SIZE, "\n")
cat("Purpose: runtime/failure/stability diagnostics, not final claims.\n")
cat("============================================================\n\n")

all_timing <- list()

for (K_now in K_GRID) {
  cfg <- make_config(n_reps = N_PILOT, K = K_now)
  saveRDS(cfg, file.path(OUTPUT_DIR, sprintf("PILOT_K%03d_CONFIG.rds", K_now)))

  k_dir <- file.path(OUTPUT_DIR, sprintf("K%03d", K_now))
  dir.create(k_dir, showWarnings = FALSE, recursive = TRUE)

  ids_all <- seq_len(N_PILOT)
  batches <- split(ids_all, ceiling(ids_all / BATCH_SIZE))

  cat("\n------------------------------------------------------------\n")
  cat("PILOT K =", K_now, "| random baseline =", sprintf("%.4f", 1 / K_now), "\n")
  cat("------------------------------------------------------------\n")

  cl <- NULL
  if (N_WORKERS > 1L) {
    cl <- parallel::makeCluster(min(N_WORKERS, N_PILOT))
    parallel::clusterEvalQ(cl, {
      library(glmnet)
      library(mvtnorm)
      NULL
    })
    parallel::clusterCall(cl, function(path) {
      source(path, local = .GlobalEnv)
      NULL
    }, common_file)
  }

  k_start <- Sys.time()
  newly_completed <- 0L
  already_completed <- 0L

  for (b in seq_along(batches)) {
    ids <- batches[[b]]
    batch_file <- file.path(k_dir, sprintf("batch_%03d.rds", b))

    if (file.exists(batch_file)) {
      old <- tryCatch(readRDS(batch_file), error = function(e) NULL)
      if (!is.null(old) && nrow(old) == length(ids) && all(as.integer(old$rep_id) == ids)) {
        already_completed <- already_completed + length(ids)
        cat(sprintf("K=%d batch %03d/%03d already exists; skipped.\n", K_now, b, length(batches)))
        next
      } else {
        stop("Incomplete/corrupt pilot batch: ", batch_file)
      }
    }

    b0 <- Sys.time()
    if (is.null(cl)) {
      rows <- lapply(ids, function(i) run_one_rep_K(i, cfg, return_debug = FALSE))
    } else {
      rows <- parallel::parLapplyLB(
        cl,
        ids,
        function(i, cfg_local) run_one_rep_K(i, cfg_local, return_debug = FALSE),
        cfg_local = cfg
      )
    }

    batch_df <- lists_to_data_frame(rows)
    saveRDS(batch_df, batch_file)
    newly_completed <- newly_completed + length(ids)

    bsec <- as.numeric(difftime(Sys.time(), b0, units = "secs"))
    elapsed <- as.numeric(difftime(Sys.time(), k_start, units = "secs"))
    done <- already_completed + newly_completed
    eta <- if (newly_completed > 0L) elapsed / newly_completed * (N_PILOT - done) else NA_real_

    cat(sprintf(
      "K=%3d batch %03d/%03d | reps %d-%d | %s | %.1f%% | ETA %s\n",
      K_now, b, length(batches), min(ids), max(ids), format_seconds(bsec),
      100 * done / N_PILOT, format_seconds(eta)
    ))
  }

  if (!is.null(cl)) parallel::stopCluster(cl)

  expected <- file.path(k_dir, sprintf("batch_%03d.rds", seq_along(batches)))
  if (any(!file.exists(expected))) stop("Pilot K=", K_now, " incomplete.")

  dfs <- lapply(expected, readRDS)
  results <- do.call(rbind, dfs)
  results <- results[order(as.integer(results$rep_id)), , drop = FALSE]

  out_csv <- file.path(OUTPUT_DIR, sprintf("PILOT_K%03d_ALL_RESULTS.csv", K_now))
  out_rds <- file.path(OUTPUT_DIR, sprintf("PILOT_K%03d_ALL_RESULTS.rds", K_now))
  utils::write.csv(results, out_csv, row.names = FALSE)
  saveRDS(results, out_rds)

  k_sec <- as.numeric(difftime(Sys.time(), k_start, units = "secs"))
  valid <- results[results$status == "ok", , drop = FALSE]
  n_valid <- nrow(valid)
  n_fail <- nrow(results) - n_valid

  cat("K=", K_now, " complete | runtime ", format_seconds(k_sec),
      " | valid ", n_valid, " | failed ", n_fail, "\n", sep = "")

  if (n_valid > 0L) {
    cat("  mean external support size =", round(mean(valid$mean_external_support_size, na.rm = TRUE), 2), "\n")
    cat("  mean Hessian regularization rate =", round(mean(valid$prop_hessian_regularized, na.rm = TRUE), 4), "\n")
    cat("  R2-(X,Y) member accuracy [PILOT] =", round(mean(valid$member_credit_R2_XY, na.rm = TRUE), 4), "\n")
  }

  if (n_fail > 0L) {
    cat("  failure classes:\n")
    print(sort(table(results$failure_class[results$status == "failed"]), decreasing = TRUE))
  }

  all_timing[[length(all_timing) + 1L]] <- data.frame(
    K = K_now,
    N = N_PILOT,
    runtime_seconds = k_sec,
    seconds_per_attempt_wall = k_sec / N_PILOT,
    valid = n_valid,
    failed = n_fail,
    failure_rate = n_fail / N_PILOT,
    stringsAsFactors = FALSE
  )
}

timing <- do.call(rbind, all_timing)
utils::write.csv(timing, file.path(OUTPUT_DIR, "PILOT_RUNTIME_BY_K.csv"), row.names = FALSE)

cat("\n============================================================\n")
cat("PILOT K-GRID COMPLETE\n")
cat("Output:", normalizePath(OUTPUT_DIR), "\n")
cat("Next: run 06_SUMMARIZE_KGRID.R on pilot results.\n")
cat("Then choose the formal N based on runtime and Monte Carlo precision.\n")
cat("============================================================\n")
