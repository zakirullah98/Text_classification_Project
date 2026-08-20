# ============================================================
# FORMAL RUN: Frozen Protocol v1.1
# N = 10,000 PRE-NUMBERED independent Monte Carlo draws.
# Checkpointed and resumable by batch.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

source("01_FedFDR_SourceAttack_functions.R")
check_required_packages()

# ============================================================
# USER-ADJUSTABLE COMPUTING SETTINGS ONLY
# These settings change runtime, NOT the statistical protocol.
# ============================================================
N_WORKERS <- 8L
BATCH_SIZE <- 100L
OUTPUT_DIR <- "results_formal_R10000"
# ============================================================

cfg <- make_config(n_reps = 10000L)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

saveRDS(cfg, file.path(OUTPUT_DIR, "config.rds"))

cat("\n============================================================\n")
cat("FORMAL Fed-FDR Source-Site Privacy Experiment\n")
cat("Protocol:", cfg$protocol_version, "\n")
cat("K=2 | n1=n2=50 | d=50 | N=10,000\n")
cat("workers =", N_WORKERS, "| batch size =", BATCH_SIZE, "\n")
cat("Output =", normalizePath(OUTPUT_DIR), "\n")
cat("============================================================\n\n")

common_file <- normalizePath("01_FedFDR_SourceAttack_functions.R")
all_ids <- seq_len(cfg$n_reps)
batch_ids <- split(all_ids, ceiling(all_ids / BATCH_SIZE))

cl <- NULL
if (N_WORKERS > 1L) {
  cl <- parallel::makeCluster(N_WORKERS)
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

start <- Sys.time()
newly_completed <- 0L
already_completed <- 0L

for (b in seq_along(batch_ids)) {
  ids <- batch_ids[[b]]
  batch_file <- file.path(OUTPUT_DIR, sprintf("batch_%03d.rds", b))

  # Resume support: do not overwrite an already completed batch.
  if (file.exists(batch_file)) {
    old <- tryCatch(readRDS(batch_file), error = function(e) NULL)
    if (!is.null(old) && nrow(old) == length(ids) && all(as.integer(old$rep_id) == ids)) {
      already_completed <- already_completed + length(ids)
      cat(sprintf("Batch %03d/%03d already exists; skipped.\n", b, length(batch_ids)))
      next
    } else {
      stop("Existing batch file is incomplete/corrupt: ", batch_file,
           ". Move or delete it after inspection before resuming.")
    }
  }

  b0 <- Sys.time()
  if (is.null(cl)) {
    rows <- lapply(ids, function(i) run_one_rep(i, cfg, return_debug = FALSE))
  } else {
    rows <- parallel::parLapplyLB(
      cl,
      ids,
      function(i, cfg_local) run_one_rep(i, cfg_local, return_debug = FALSE),
      cfg_local = cfg
    )
  }

  batch_df <- lists_to_data_frame(rows)
  saveRDS(batch_df, batch_file)
  newly_completed <- newly_completed + length(ids)

  bsec <- as.numeric(difftime(Sys.time(), b0, units = "secs"))
  elapsed_new <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  effective_done <- already_completed + newly_completed

  if (newly_completed > 0L) {
    sec_per_new_rep <- elapsed_new / newly_completed
    remaining <- cfg$n_reps - effective_done
    eta <- sec_per_new_rep * remaining
  } else {
    eta <- NA_real_
  }

  valid_b <- sum(batch_df$status == "ok", na.rm = TRUE)
  fail_b <- sum(batch_df$status == "failed", na.rm = TRUE)

  cat(sprintf(
    paste0("Batch %03d/%03d | reps %d-%d | %.1f sec | valid %d fail %d | ",
           "overall %.1f%% | ETA %s\n"),
    b, length(batch_ids), min(ids), max(ids), bsec, valid_b, fail_b,
    100 * effective_done / cfg$n_reps, format_seconds(eta)
  ))
}

if (!is.null(cl)) {
  parallel::stopCluster(cl)
  cl <- NULL
}

# Combine all checkpoint files only after every batch exists.
expected_files <- file.path(OUTPUT_DIR, sprintf("batch_%03d.rds", seq_along(batch_ids)))
missing_files <- expected_files[!file.exists(expected_files)]
if (length(missing_files) > 0L) {
  stop("Formal run is incomplete. Missing batch files: ", paste(basename(missing_files), collapse = ", "))
}

all_batches <- lapply(expected_files, readRDS)
results <- do.call(rbind, all_batches)
results <- results[order(as.integer(results$rep_id)), , drop = FALSE]

if (nrow(results) != cfg$n_reps || !all(as.integer(results$rep_id) == seq_len(cfg$n_reps))) {
  stop("Combined result failed the 1:10,000 repetition-ID audit.")
}

saveRDS(results, file.path(OUTPUT_DIR, "FORMAL_R10000_ALL_RESULTS.rds"))
utils::write.csv(results, file.path(OUTPUT_DIR, "FORMAL_R10000_ALL_RESULTS.csv"), row.names = FALSE)

valid_n <- sum(results$status == "ok", na.rm = TRUE)
fail_n <- sum(results$status == "failed", na.rm = TRUE)
manifest <- list(
  protocol_version = cfg$protocol_version,
  completed_at = as.character(Sys.time()),
  attempted = cfg$n_reps,
  valid = valid_n,
  failed = fail_n,
  failure_rate = fail_n / cfg$n_reps,
  n_workers = N_WORKERS,
  batch_size = BATCH_SIZE
)
saveRDS(manifest, file.path(OUTPUT_DIR, "run_manifest.rds"))

cat("\n============================================================\n")
cat("FORMAL RUN COMPLETE\n")
cat("Attempted:", cfg$n_reps, "\n")
cat("Valid:    ", valid_n, "\n")
cat("Failed:   ", fail_n, "\n")
cat("Failure rate:", round(fail_n / cfg$n_reps, 6), "\n")
cat("============================================================\n")

if (fail_n > 0L) {
  cat("\nFailure stages:\n")
  print(sort(table(results$failure_stage[results$status == "failed"]), decreasing = TRUE))
}

cat("\nNext: run 05_ANALYZE_RESULTS.R\n")
