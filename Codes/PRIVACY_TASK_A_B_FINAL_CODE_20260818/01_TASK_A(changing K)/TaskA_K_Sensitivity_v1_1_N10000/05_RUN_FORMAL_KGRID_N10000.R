# ============================================================
# Task A FORMAL K-grid: N=10,000 per K
# ============================================================
# Do NOT run this until:
#   1) 02_TEST_K2_BACKWARD_COMPATIBILITY.R passes;
#   2) 03_RUN_QUICK_KGRID.R is inspected;
#   3) 04_RUN_PILOT_KGRID.R is inspected;
#   4) N_FORMAL is chosen from pilot runtime / desired MC precision.
#
# Formal protocol fixed for the teacher-requested N=10,000 attempted repetitions
# for each K in the grid. Checkpoints are saved batch-by-batch so the run can resume.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

source("01_SourceAttack_Kgeneral_functions.R")
check_required_packages()

# ============================================================
# FROZEN FORMAL SETTINGS
# ============================================================
K_GRID <- c(2L, 3L, 5L, 10L, 20L, 50L, 100L)
N_FORMAL <- 10000L
N_WORKERS <- 8L
BATCH_SIZE <- 50L
OUTPUT_DIR <- "results_formal_Kgrid_N10000"
# ============================================================

if (!file.exists("K2_backward_compatibility_report.csv")) {
  stop("Missing K=2 backward-compatibility report.")
}
backward <- utils::read.csv("K2_backward_compatibility_report.csv", stringsAsFactors = FALSE)
if (!all(backward$pass)) stop("K=2 backward compatibility failed.")

if (!dir.exists("results_pilot_Kgrid")) {
  stop("Missing pilot results. Run 04_RUN_PILOT_KGRID.R first.")
}

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
common_file <- normalizePath("01_SourceAttack_Kgeneral_functions.R")

cat("\n============================================================\n")
cat("Task A FORMAL K-GRID\n")
cat("K grid:", paste(K_GRID, collapse = ", "), "\n")
cat("N formal per K:", N_FORMAL, "\n")
cat("n_k=50 | p=50 | baseline DGP unchanged\n")
cat("workers:", N_WORKERS, "| batch size:", BATCH_SIZE, "\n")
cat("============================================================\n\n")

all_timing <- list()

for (K_now in K_GRID) {
  cfg <- make_config(n_reps = N_FORMAL, K = K_now)
  saveRDS(cfg, file.path(OUTPUT_DIR, sprintf("FORMAL_K%03d_CONFIG.rds", K_now)))

  k_dir <- file.path(OUTPUT_DIR, sprintf("K%03d", K_now))
  dir.create(k_dir, showWarnings = FALSE, recursive = TRUE)

  ids_all <- seq_len(N_FORMAL)
  batches <- split(ids_all, ceiling(ids_all / BATCH_SIZE))

  cat("\n------------------------------------------------------------\n")
  cat("FORMAL K =", K_now, "| baseline =", sprintf("%.6f", 1 / K_now), "\n")
  cat("------------------------------------------------------------\n")

  cl <- NULL
  if (N_WORKERS > 1L) {
    cl <- parallel::makeCluster(min(N_WORKERS, N_FORMAL))
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
    batch_file <- file.path(k_dir, sprintf("batch_%04d.rds", b))

    if (file.exists(batch_file)) {
      old <- tryCatch(readRDS(batch_file), error = function(e) NULL)
      if (!is.null(old) && nrow(old) == length(ids) && all(as.integer(old$rep_id) == ids)) {
        already_completed <- already_completed + length(ids)
        cat(sprintf("K=%d batch %04d/%04d already exists; skipped.\n", K_now, b, length(batches)))
        next
      } else {
        stop("Incomplete/corrupt formal batch: ", batch_file)
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
    eta <- if (newly_completed > 0L) elapsed / newly_completed * (N_FORMAL - done) else NA_real_

    cat(sprintf(
      "K=%3d batch %04d/%04d | reps %d-%d | %s | %.1f%% | ETA %s\n",
      K_now, b, length(batches), min(ids), max(ids), format_seconds(bsec),
      100 * done / N_FORMAL, format_seconds(eta)
    ))
  }

  if (!is.null(cl)) parallel::stopCluster(cl)

  expected <- file.path(k_dir, sprintf("batch_%04d.rds", seq_along(batches)))
  if (any(!file.exists(expected))) stop("Formal K=", K_now, " incomplete.")

  dfs <- lapply(expected, readRDS)
  results <- do.call(rbind, dfs)
  results <- results[order(as.integer(results$rep_id)), , drop = FALSE]

  if (nrow(results) != N_FORMAL || !all(as.integer(results$rep_id) == seq_len(N_FORMAL))) {
    stop("Formal repetition-ID audit failed for K=", K_now)
  }

  out_csv <- file.path(OUTPUT_DIR, sprintf("FORMAL_K%03d_ALL_RESULTS.csv", K_now))
  out_rds <- file.path(OUTPUT_DIR, sprintf("FORMAL_K%03d_ALL_RESULTS.rds", K_now))
  utils::write.csv(results, out_csv, row.names = FALSE)
  saveRDS(results, out_rds)

  k_sec <- as.numeric(difftime(Sys.time(), k_start, units = "secs"))
  n_valid <- sum(results$status == "ok", na.rm = TRUE)
  n_fail <- N_FORMAL - n_valid

  all_timing[[length(all_timing) + 1L]] <- data.frame(
    K = K_now,
    N = N_FORMAL,
    runtime_seconds = k_sec,
    valid = n_valid,
    failed = n_fail,
    failure_rate = n_fail / N_FORMAL,
    stringsAsFactors = FALSE
  )

  cat("K=", K_now, " complete | runtime ", format_seconds(k_sec),
      " | valid ", n_valid, " | failed ", n_fail, "\n", sep = "")
}

timing <- do.call(rbind, all_timing)
utils::write.csv(timing, file.path(OUTPUT_DIR, "FORMAL_RUNTIME_BY_K.csv"), row.names = FALSE)

cat("\n============================================================\n")
cat("FORMAL K-GRID COMPLETE\n")
cat("Output:", normalizePath(OUTPUT_DIR), "\n")
cat("Next: run 06_SUMMARIZE_TASKA_N10000.R.\n")
cat("============================================================\n")
