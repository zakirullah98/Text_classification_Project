# ============================================================
# QUICK TEST: Frozen Protocol v1.1
# N = 20, single worker by default for transparent debugging.
# Run this FIRST.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

source("01_FedFDR_SourceAttack_functions.R")
check_required_packages()

cfg <- make_config(n_reps = 20L)

cat("\n========================================\n")
cat("Fed-FDR Source-Site Attack QUICK TEST\n")
cat("Protocol:", cfg$protocol_version, "\n")
cat("K=", cfg$K, ", n/site=", cfg$n_per_site,
    ", d=", cfg$d, ", s=", cfg$s,
    ", N=", cfg$n_reps, "\n", sep = "")
cat("========================================\n\n")

start <- Sys.time()
rows <- vector("list", cfg$n_reps)

for (r in seq_len(cfg$n_reps)) {
  t0 <- Sys.time()
  rows[[r]] <- run_one_rep(r, cfg, return_debug = FALSE)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  status <- rows[[r]]$status %||% "unknown"
  cat(sprintf("Rep %3d/%d | %-6s | %.2f sec\n", r, cfg$n_reps, status, dt))
}

results <- lists_to_data_frame(rows)

out_dir <- file.path("results_quick_R20")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(results, file.path(out_dir, "quick_R20_results.rds"))
utils::write.csv(results, file.path(out_dir, "quick_R20_results.csv"), row.names = FALSE)
saveRDS(cfg, file.path(out_dir, "config.rds"))

elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
cat("\nFinished in", format_seconds(elapsed), "\n")
cat("Successful:", sum(results$status == "ok", na.rm = TRUE), "\n")
cat("Failed:    ", sum(results$status == "failed", na.rm = TRUE), "\n")

if (any(results$status == "failed", na.rm = TRUE)) {
  cat("\nFailure table:\n")
  print(sort(table(results$failure_stage[results$status == "failed"]), decreasing = TRUE))
}

ok <- results[results$status == "ok", , drop = FALSE]
if (nrow(ok) > 0L) {
  cat("\nQuick member accuracies (NOT formal results):\n")
  conds <- c("R0_X", "R0_XY", "R1_X", "R1_XY", "R2_X", "R2_XY", "R3_X", "R3_XY")
  quick <- do.call(rbind, lapply(conds, function(cc) summarize_condition(ok, cc)))
  print(quick[, c("condition", "member_accuracy", "fresh_accuracy", "delta_train")], row.names = FALSE)

  cat("\nProtocol checks:\n")
  cat("Mean |S1hat| =", mean(ok$support1_size), "\n")
  cat("Mean |S2hat| =", mean(ok$support2_size), "\n")
  cat("Mean Fed-FDR FDP =", mean(ok$fed_fdp), "\n")
  cat("Mean Fed-FDR Power =", mean(ok$fed_power), "\n")
  cat("Hessian regularization rate site1 =", mean(ok$hessian_regularized_site1), "\n")
  cat("Hessian regularization rate site2 =", mean(ok$hessian_regularized_site2), "\n")
}

cat("\nSaved to:", normalizePath(out_dir), "\n")
cat("DO NOT launch N=10,000 until this quick test is inspected.\n")
