# ============================================================
# QUICK AUDIT: paper-faithful Fed-FDR central step, R = 20
# ============================================================
rm(list = ls())
options(stringsAsFactors = FALSE)

source("07_PAPER_FAITHFUL_CENTRAL_AUDIT_FUNCTIONS.R")
check_required_packages()

N_WORKERS <- 4L
OUTPUT_DIR <- "results_central_audit_R20"

cfg <- make_config(n_reps = 20L)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
saveRDS(cfg, file.path(OUTPUT_DIR, "config.rds"))

ids <- seq_len(cfg$n_reps)
common_files <- c(
  normalizePath("01_FedFDR_SourceAttack_functions.R"),
  normalizePath("07_PAPER_FAITHFUL_CENTRAL_AUDIT_FUNCTIONS.R")
)

cat("\nPaper-faithful central-step QUICK audit\n")
cat("R=20 | K=2 | n=50/site | d=50 | alpha=0.1\n")
cat("Comparing A=CODE_THETA, B=CODE_RAWBETA, C=PAPER_RAWBETA\n\n")

t0 <- Sys.time()
if (N_WORKERS > 1L) {

  # Keep cluster creation/cleanup inside a function so on.exit() is
  # scoped correctly. This prevents the cluster from being stopped twice.
  run_parallel_quick <- function(ids0, cfg0, files0, n_workers0) {
    cl <- parallel::makeCluster(n_workers0)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

    parallel::clusterEvalQ(cl, {
      library(glmnet)
      library(mvtnorm)
      NULL
    })

    parallel::clusterCall(cl, function(files) {
      source(files[1L], local = .GlobalEnv)
      source(files[2L], local = .GlobalEnv)
      NULL
    }, files0)

    parallel::parLapplyLB(
      cl, ids0,
      function(i, cfg1) run_central_audit_rep(i, cfg1),
      cfg1 = cfg0
    )
  }

  rows <- run_parallel_quick(ids, cfg, common_files, N_WORKERS)

} else {
  rows <- lapply(ids, function(i) run_central_audit_rep(i, cfg))
}

res <- audit_rows_to_df(rows)
saveRDS(res, file.path(OUTPUT_DIR, "AUDIT_R20_RESULTS.rds"))
utils::write.csv(res, file.path(OUTPUT_DIR, "AUDIT_R20_RESULTS.csv"), row.names = FALSE)

cat("Completed in", fmt_sec(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "\n")
cat("Valid:", sum(res$status == "ok"), " Failed:", sum(res$status == "failed"), "\n")
cat("Output:", normalizePath(OUTPUT_DIR), "\n")
cat("Next: inspect R20; if clean, run 09_RUN_AUDIT_FORMAL_R10000.R\n")
