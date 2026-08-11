# ============================================================
# FORMAL AUDIT: paper-faithful Fed-FDR central step, R = 10,000
# Checkpointed and resumable.
# ============================================================
rm(list = ls())
options(stringsAsFactors = FALSE)

source("07_PAPER_FAITHFUL_CENTRAL_AUDIT_FUNCTIONS.R")
check_required_packages()

# USER-ADJUSTABLE COMPUTING SETTINGS ONLY
N_WORKERS <- 8L
BATCH_SIZE <- 100L
OUTPUT_DIR <- "results_central_audit_R10000"

cfg <- make_config(n_reps = 10000L)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
saveRDS(cfg, file.path(OUTPUT_DIR, "config.rds"))

cat("\n============================================================\n")
cat("FORMAL Paper-Faithful Fed-FDR Central-Step Audit\n")
cat("Same seeds/upstream fitting as Frozen Protocol v1.1\n")
cat("A = CODE_THETA     : theta_mirror + public-code central logic\n")
cat("B = CODE_RAWBETA   : raw refined beta + public-code central logic\n")
cat("C = PAPER_RAWBETA  : raw refined beta + paper Algorithm 1\n")
cat("K=2 | n1=n2=50 | d=50 | N=10,000 | alpha=0.1\n")
cat("workers =", N_WORKERS, "| batch size =", BATCH_SIZE, "\n")
cat("============================================================\n\n")

all_ids <- seq_len(cfg$n_reps)
batch_ids <- split(all_ids, ceiling(all_ids / BATCH_SIZE))
common_files <- c(
  normalizePath("01_FedFDR_SourceAttack_functions.R"),
  normalizePath("07_PAPER_FAITHFUL_CENTRAL_AUDIT_FUNCTIONS.R")
)

cl <- NULL
if (N_WORKERS > 1L) {
  cl <- parallel::makeCluster(N_WORKERS)
  parallel::clusterEvalQ(cl, { library(glmnet); library(mvtnorm); NULL })
  parallel::clusterCall(cl, function(files) {
    source(files[1L], local = .GlobalEnv)
    source(files[2L], local = .GlobalEnv)
    NULL
  }, common_files)
}

start <- Sys.time()
new_done <- 0L
old_done <- 0L

for (b in seq_along(batch_ids)) {
  ids <- batch_ids[[b]]
  bf <- file.path(OUTPUT_DIR, sprintf("batch_%03d.rds", b))

  if (file.exists(bf)) {
    old <- tryCatch(readRDS(bf), error = function(e) NULL)
    if (!is.null(old) && nrow(old) == length(ids) && all(as.integer(old$rep_id) == ids)) {
      old_done <- old_done + length(ids)
      cat(sprintf("Batch %03d/%03d exists; skipped.\n", b, length(batch_ids)))
      next
    }
    stop("Incomplete/corrupt checkpoint: ", bf)
  }

  b0 <- Sys.time()
  if (is.null(cl)) {
    rows <- lapply(ids, function(i) run_central_audit_rep(i, cfg))
  } else {
    rows <- parallel::parLapplyLB(
      cl, ids,
      function(i, cfg0) run_central_audit_rep(i, cfg0),
      cfg0 = cfg
    )
  }
  df <- audit_rows_to_df(rows)
  saveRDS(df, bf)
  new_done <- new_done + length(ids)

  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  total_done <- old_done + new_done
  eta <- if (new_done > 0L) elapsed / new_done * (cfg$n_reps - total_done) else NA_real_
  cat(sprintf(
    "Batch %03d/%03d | reps %d-%d | valid %d fail %d | %.1f%% | ETA %s\n",
    b, length(batch_ids), min(ids), max(ids),
    sum(df$status == "ok"), sum(df$status == "failed"),
    100 * total_done / cfg$n_reps, fmt_sec(eta)
  ))
}

if (!is.null(cl)) parallel::stopCluster(cl)

expected <- file.path(OUTPUT_DIR, sprintf("batch_%03d.rds", seq_along(batch_ids)))
if (any(!file.exists(expected))) stop("Audit incomplete: missing checkpoint batches.")

res <- do.call(rbind, lapply(expected, readRDS))
res <- res[order(as.integer(res$rep_id)), , drop = FALSE]
if (nrow(res) != 10000L || !all(as.integer(res$rep_id) == seq_len(10000L))) {
  stop("1:10,000 repetition-ID audit failed.")
}

saveRDS(res, file.path(OUTPUT_DIR, "AUDIT_R10000_ALL_RESULTS.rds"))
utils::write.csv(res, file.path(OUTPUT_DIR, "AUDIT_R10000_ALL_RESULTS.csv"), row.names = FALSE)

manifest <- list(
  completed_at = as.character(Sys.time()),
  attempted = nrow(res),
  valid = sum(res$status == "ok"),
  failed = sum(res$status == "failed"),
  n_workers = N_WORKERS,
  batch_size = BATCH_SIZE,
  variants = c("A_CODE_THETA", "B_CODE_RAWBETA", "C_PAPER_RAWBETA")
)
saveRDS(manifest, file.path(OUTPUT_DIR, "run_manifest.rds"))

cat("\n============================================================\n")
cat("AUDIT COMPLETE\n")
cat("Attempted:", nrow(res), "\n")
cat("Valid:    ", sum(res$status == "ok"), "\n")
cat("Failed:   ", sum(res$status == "failed"), "\n")
cat("Next: run 10_ANALYZE_CENTRAL_AUDIT.R\n")
cat("============================================================\n")
