# ============================================================
# PILOT: Frozen Protocol v1.1
# N = 100. Estimates failure rate and runtime before formal run.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

source("01_FedFDR_SourceAttack_functions.R")
check_required_packages()

# ---- User-adjustable computing settings ----
N_WORKERS <- 4L
BATCH_SIZE <- 20L
# -------------------------------------------

cfg <- make_config(n_reps = 100L)
out_dir <- file.path("results_pilot_R100")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n========================================\n")
cat("Fed-FDR Source-Site Attack PILOT\n")
cat("Protocol:", cfg$protocol_version, "\n")
cat("N =", cfg$n_reps, "| workers =", N_WORKERS, "| batch size =", BATCH_SIZE, "\n")
cat("========================================\n\n")

common_file <- normalizePath("01_FedFDR_SourceAttack_functions.R")

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

all_rows <- list()
start <- Sys.time()
completed <- 0L
batch_ids <- split(seq_len(cfg$n_reps), ceiling(seq_len(cfg$n_reps) / BATCH_SIZE))

for (b in seq_along(batch_ids)) {
  ids <- batch_ids[[b]]
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
  all_rows[[b]] <- batch_df
  saveRDS(batch_df, file.path(out_dir, sprintf("batch_%03d.rds", b)))

  completed <- completed + length(ids)
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  rate <- elapsed / completed
  eta <- rate * (cfg$n_reps - completed)
  bsec <- as.numeric(difftime(Sys.time(), b0, units = "secs"))

  cat(sprintf(
    "Batch %d/%d | reps %d-%d | %.1f sec | progress %.1f%% | ETA %s\n",
    b, length(batch_ids), min(ids), max(ids), bsec,
    100 * completed / cfg$n_reps, format_seconds(eta)
  ))
}

if (!is.null(cl)) {
  parallel::stopCluster(cl)
  cl <- NULL
}

results <- do.call(rbind, all_rows)
saveRDS(results, file.path(out_dir, "pilot_R100_results.rds"))
utils::write.csv(results, file.path(out_dir, "pilot_R100_results.csv"), row.names = FALSE)
saveRDS(cfg, file.path(out_dir, "config.rds"))

elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
valid_n <- sum(results$status == "ok", na.rm = TRUE)
fail_n <- sum(results$status == "failed", na.rm = TRUE)

cat("\nFinished in", format_seconds(elapsed), "\n")
cat("Valid:", valid_n, "| Failed:", fail_n,
    "| Failure rate:", round(fail_n / cfg$n_reps, 4), "\n")

if (valid_n > 0L) {
  sec_per_rep <- elapsed / cfg$n_reps
  estimated_formal <- sec_per_rep * 10000L
  cat("Crude N=10,000 runtime projection at this worker setting:",
      format_seconds(estimated_formal), "\n")
}

if (fail_n > 0L) {
  cat("\nFailure stages:\n")
  print(sort(table(results$failure_stage[results$status == "failed"]), decreasing = TRUE))
  cat("\nFailure reasons (top):\n")
  print(head(sort(table(results$failure_reason[results$status == "failed"]), decreasing = TRUE), 10))
}

cat("\nSaved to:", normalizePath(out_dir), "\n")
cat("Inspect failure rate before formal N=10,000.\n")
