# ============================================================
# OPTIONAL Task-B sensitivity template
# Do NOT run automatically. Use only if P1 or P4 main result needs explanation.
# Scenarios:
#   P1_IID_UNIFORM
#   P4_POISSON_RAW
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)
source("01_TaskB_PopulationAttack_functions.R")
check_required_packages()

# USER CHOICE ------------------------------------------------
STAGE <- "QUICK"  # "QUICK", "PILOT", or "FORMAL"
POPULATIONS <- TASKB_SENSITIVITY_POPULATIONS
N_WORKERS <- 8L
BASE_SEED <- 20260815L
# ------------------------------------------------------------

if (STAGE == "QUICK") {
  N_REPS <- 5L; DIAGNOSTIC_N <- 1000L; BATCH_SIZE <- 5L
} else if (STAGE == "PILOT") {
  N_REPS <- 100L; DIAGNOSTIC_N <- 2000L; BATCH_SIZE <- 50L
} else if (STAGE == "FORMAL") {
  N_REPS <- 10000L; DIAGNOSTIC_N <- 0L; BATCH_SIZE <- 100L
} else stop("Unknown STAGE")

OUTPUT_DIR <- paste0("results_taskB_sensitivity_", tolower(STAGE))
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
common_file <- normalizePath("01_TaskB_PopulationAttack_functions.R")

cl <- NULL
if (N_WORKERS > 1L && STAGE != "QUICK") {
  cl <- parallel::makeCluster(N_WORKERS)
  parallel::clusterEvalQ(cl, { library(glmnet); library(mvtnorm); NULL })
  parallel::clusterCall(cl, function(path) { source(path, local = .GlobalEnv); NULL }, common_file)
}

for (pop in POPULATIONS) {
  cfg <- make_taskb_config(pop, n_reps = N_REPS, base_seed = BASE_SEED,
                           diagnostic_n = DIAGNOSTIC_N)
  pdir <- file.path(OUTPUT_DIR, pop)
  dir.create(pdir, showWarnings = FALSE, recursive = TRUE)

  ids <- seq_len(N_REPS)
  rows <- if (is.null(cl)) {
    lapply(ids, function(i) run_one_rep_taskb(i, cfg, FALSE))
  } else {
    parallel::parLapplyLB(cl, ids,
      function(i, cfg_local) run_one_rep_taskb(i, cfg_local, FALSE),
      cfg_local = cfg)
  }
  df <- lists_to_data_frame(rows)
  saveRDS(df, file.path(pdir, paste0(STAGE, "_results.rds")))
  utils::write.csv(df, file.path(pdir, paste0(STAGE, "_results.csv")), row.names = FALSE)
  utils::write.csv(summarize_taskb_attacks(df), file.path(pdir, "attack_summary.csv"), row.names = FALSE)
  utils::write.csv(summarize_taskb_dgp(df), file.path(pdir, "dgp_summary.csv"), row.names = FALSE)
  cat(pop, ": valid", sum(df$status == "ok"), "failed", sum(df$status == "failed"), "\n")
}

if (!is.null(cl)) parallel::stopCluster(cl)
cat("Sensitivity run complete. These scenarios are secondary, not part of the frozen five-population main table.\n")
