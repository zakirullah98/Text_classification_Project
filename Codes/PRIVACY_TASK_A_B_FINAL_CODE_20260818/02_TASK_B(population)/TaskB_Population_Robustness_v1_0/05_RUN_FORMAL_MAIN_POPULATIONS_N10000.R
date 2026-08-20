# ============================================================
# Task B FORMAL: N=10,000 attempted repetitions PER MAIN population
# Checkpointed and resumable by population and batch.
# IMPORTANT: oracle diagnostic sample is OFF in formal runs; it is a Pilot
# diagnostic only and is not needed for the attack estimands.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)
source("01_TaskB_PopulationAttack_functions.R")
check_required_packages()

N_FORMAL <- 10000L
N_WORKERS <- 8L
BATCH_SIZE <- 100L
BASE_SEED <- 20260815L
OUTPUT_DIR <- "results_taskB_formal_N10000"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
common_file <- normalizePath("01_TaskB_PopulationAttack_functions.R")

cat("\n============================================================\n")
cat("TASK B FORMAL -- FIVE MAIN POPULATIONS\n")
cat("K=2 | n1=n2=50 | p=50 | N=10,000 per population\n")
cat("Total attempted repetitions = 50,000\n")
cat("workers=", N_WORKERS, " | batch size=", BATCH_SIZE, "\n", sep = "")
cat("ATTACK FIXED; ONLY X POPULATION CHANGES\n")
cat("============================================================\n")

cl <- NULL
if (N_WORKERS > 1L) {
  cl <- parallel::makeCluster(N_WORKERS)
  parallel::clusterEvalQ(cl, { library(glmnet); library(mvtnorm); NULL })
  parallel::clusterCall(cl, function(path) { source(path, local = .GlobalEnv); NULL }, common_file)
}

formal_manifest <- list()

for (pop in TASKB_MAIN_POPULATIONS) {
  cfg <- make_taskb_config(pop, n_reps = N_FORMAL, base_seed = BASE_SEED,
                           diagnostic_n = 0L)
  pdir <- file.path(OUTPUT_DIR, pop)
  dir.create(pdir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(cfg, file.path(pdir, "config.rds"))

  cat("\n============================================================\n")
  cat("FORMAL:", cfg$population_label, "\n")
  cat("============================================================\n")

  all_ids <- seq_len(N_FORMAL)
  batch_ids <- split(all_ids, ceiling(all_ids / BATCH_SIZE))
  start <- Sys.time()
  newly_completed <- 0L
  already_completed <- 0L

  for (b in seq_along(batch_ids)) {
    ids <- batch_ids[[b]]
    batch_file <- file.path(pdir, sprintf("batch_%03d.rds", b))

    if (file.exists(batch_file)) {
      old <- tryCatch(readRDS(batch_file), error = function(e) NULL)
      if (!is.null(old) && nrow(old) == length(ids) &&
          all(as.integer(old$rep_id) == ids)) {
        already_completed <- already_completed + length(ids)
        cat(sprintf("[%s] batch %03d/%03d exists; skipped.\n", pop, b, length(batch_ids)))
        next
      } else {
        stop("Incomplete/corrupt checkpoint: ", batch_file)
      }
    }

    b0 <- Sys.time()
    if (is.null(cl)) {
      rows <- lapply(ids, function(i) run_one_rep_taskb(i, cfg, FALSE))
    } else {
      rows <- parallel::parLapplyLB(
        cl, ids,
        function(i, cfg_local) run_one_rep_taskb(i, cfg_local, FALSE),
        cfg_local = cfg
      )
    }
    batch_df <- lists_to_data_frame(rows)
    saveRDS(batch_df, batch_file)
    newly_completed <- newly_completed + length(ids)

    bsec <- as.numeric(difftime(Sys.time(), b0, units = "secs"))
    elapsed_new <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    effective_done <- already_completed + newly_completed
    eta <- if (newly_completed > 0L) {
      (elapsed_new / newly_completed) * (N_FORMAL - effective_done)
    } else NA_real_

    cat(sprintf(
      "[%s] batch %03d/%03d | reps %d-%d | %.1fs | valid %d fail %d | %.1f%% | ETA %s\n",
      pop, b, length(batch_ids), min(ids), max(ids), bsec,
      sum(batch_df$status == "ok", na.rm = TRUE),
      sum(batch_df$status == "failed", na.rm = TRUE),
      100 * effective_done / N_FORMAL, format_seconds(eta)
    ))
  }

  expected_files <- file.path(pdir, sprintf("batch_%03d.rds", seq_along(batch_ids)))
  if (any(!file.exists(expected_files))) {
    stop("Population run incomplete: ", pop)
  }

  results <- do.call(rbind, lapply(expected_files, readRDS))
  results <- results[order(as.integer(results$rep_id)), , drop = FALSE]
  if (nrow(results) != N_FORMAL || !all(as.integer(results$rep_id) == seq_len(N_FORMAL))) {
    stop("1:10,000 repetition-ID audit failed for ", pop)
  }

  saveRDS(results, file.path(pdir, "FORMAL_N10000_ALL_RESULTS.rds"))
  utils::write.csv(results, file.path(pdir, "FORMAL_N10000_ALL_RESULTS.csv"), row.names = FALSE)

  valid_n <- sum(results$status == "ok", na.rm = TRUE)
  fail_n <- sum(results$status == "failed", na.rm = TRUE)
  man <- list(
    population = pop,
    population_label = cfg$population_label,
    attempted = N_FORMAL,
    valid = valid_n,
    failed = fail_n,
    failure_rate = fail_n / N_FORMAL,
    completed_at = as.character(Sys.time()),
    n_workers = N_WORKERS,
    batch_size = BATCH_SIZE
  )
  saveRDS(man, file.path(pdir, "run_manifest.rds"))
  formal_manifest[[pop]] <- man

  cat("Completed:", pop, "| valid", valid_n, "| failed", fail_n, "\n")
}

if (!is.null(cl)) parallel::stopCluster(cl)

saveRDS(formal_manifest, file.path(OUTPUT_DIR, "TASKB_FORMAL_MANIFEST.rds"))
cfg0 <- make_taskb_config(TASKB_MAIN_POPULATIONS[1], n_reps = N_FORMAL)
write_taskb_p2_permutation(cfg0, file.path(OUTPUT_DIR, "P2_FIXED_PERMUTATION.csv"))

cat("\n============================================================\n")
cat("TASK B FORMAL COMPLETE\n")
cat("Next: run 06_SUMMARIZE_TASKB_FORMAL.R\n")
cat("============================================================\n")
