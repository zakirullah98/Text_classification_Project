# ============================================================
# Analysis for Frozen Protocol v1.1 formal results
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

source("01_FedFDR_SourceAttack_functions.R")

INPUT_DIR <- "results_formal_R10000"
input_file <- file.path(INPUT_DIR, "FORMAL_R10000_ALL_RESULTS.rds")
if (!file.exists(input_file)) {
  stop("Cannot find ", input_file, ". Finish 04_RUN_FORMAL_R10000.R first.")
}

res <- readRDS(input_file)
valid <- res[res$status == "ok", , drop = FALSE]
failed <- res[res$status == "failed", , drop = FALSE]

analysis_dir <- file.path(INPUT_DIR, "analysis")
dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n========================================\n")
cat("ANALYSIS -- Frozen Protocol v1.1\n")
cat("Attempted:", nrow(res), "\n")
cat("Valid:    ", nrow(valid), "\n")
cat("Failed:   ", nrow(failed), "\n")
cat("Failure rate:", round(nrow(failed) / nrow(res), 6), "\n")
cat("========================================\n\n")

# ------------------------------------------------------------
# Primary table
# ------------------------------------------------------------
conditions <- c("R0_X", "R0_XY", "R1_X", "R1_XY", "R2_X", "R2_XY", "R3_X", "R3_XY")
primary <- do.call(rbind, lapply(conditions, function(cc) summarize_condition(valid, cc)))

# Add B0 exact no-information baseline.
B0 <- data.frame(
  condition = "B0_ZERO_INFORMATION",
  member_accuracy = 0.5,
  attack_advantage = 0,
  member_ci_lower = 0.5,
  member_ci_upper = 0.5,
  one_sided_p_vs_0_5 = 1,
  fresh_accuracy = 0.5,
  fresh_ci_lower = 0.5,
  fresh_ci_upper = 0.5,
  delta_train = 0,
  delta_train_ci_lower = 0,
  delta_train_ci_upper = 0,
  member_auc = 0.5,
  fresh_auc = 0.5,
  n_valid = nrow(valid),
  stringsAsFactors = FALSE
)
primary <- rbind(B0, primary)

utils::write.csv(primary, file.path(analysis_dir, "01_primary_attack_results.csv"), row.names = FALSE)
saveRDS(primary, file.path(analysis_dir, "01_primary_attack_results.rds"))

cat("PRIMARY RESULTS:\n")
print(primary, row.names = FALSE, digits = 4)

# ------------------------------------------------------------
# Full-regime component decomposition
# ------------------------------------------------------------
component_summary_one <- function(label, member_score, member_credit, fresh_score, fresh_credit) {
  mci <- mean_ci(member_credit)
  fci <- mean_ci(fresh_credit)
  dci <- paired_delta_ci(member_credit, fresh_credit)
  data.frame(
    component = label,
    member_accuracy = mci["mean"],
    fresh_accuracy = fci["mean"],
    delta_train = dci["delta"],
    member_auc = auc_rank(member_score, valid$source_site),
    fresh_auc = auc_rank(fresh_score, valid$source_site),
    stringsAsFactors = FALSE
  )
}

full_components <- rbind(
  component_summary_one(
    "X: support dS",
    valid$member_dS,
    valid$member_credit_R1_X,
    valid$fresh_dS,
    valid$fresh_credit_R1_X
  ),
  component_summary_one(
    "X: coefficient dBeta",
    valid$member_dBeta_X,
    valid$member_credit_R2_X,
    valid$fresh_dBeta_X,
    valid$fresh_credit_R2_X
  ),
  component_summary_one(
    "X: equal-weight dS+dBeta",
    valid$member_score_R3_X,
    valid$member_credit_R3_X,
    valid$fresh_score_R3_X,
    valid$fresh_credit_R3_X
  ),
  component_summary_one(
    "XY: support dS",
    valid$member_dS,
    valid$member_credit_R1_XY,
    valid$fresh_dS,
    valid$fresh_credit_R1_XY
  ),
  component_summary_one(
    "XY: coefficient dBeta",
    valid$member_dBeta_XY,
    valid$member_credit_R2_XY,
    valid$fresh_dBeta_XY,
    valid$fresh_credit_R2_XY
  ),
  component_summary_one(
    "XY: equal-weight dS+dBeta",
    valid$member_score_R3_XY,
    valid$member_credit_R3_XY,
    valid$fresh_score_R3_XY,
    valid$fresh_credit_R3_XY
  )
)
utils::write.csv(full_components, file.path(analysis_dir, "02_full_component_results.csv"), row.names = FALSE)

# Agreement / conflict / component ties
relation_member_x <- summarize_relation(valid$member_relation_X)
relation_member_xy <- summarize_relation(valid$member_relation_XY)
relation_fresh_x <- summarize_relation(valid$fresh_relation_X)
relation_fresh_xy <- summarize_relation(valid$fresh_relation_XY)

relations <- data.frame(
  sample_type = c("member", "member", "fresh", "fresh"),
  candidate_info = c("X", "XY", "X", "XY"),
  agreement = c(relation_member_x["agreement"], relation_member_xy["agreement"],
                relation_fresh_x["agreement"], relation_fresh_xy["agreement"]),
  conflict = c(relation_member_x["conflict"], relation_member_xy["conflict"],
               relation_fresh_x["conflict"], relation_fresh_xy["conflict"]),
  component_tie = c(relation_member_x["component_tie"], relation_member_xy["component_tie"],
                    relation_fresh_x["component_tie"], relation_fresh_xy["component_tie"]),
  n = c(relation_member_x["n"], relation_member_xy["n"], relation_fresh_x["n"], relation_fresh_xy["n"]),
  stringsAsFactors = FALSE
)
utils::write.csv(relations, file.path(analysis_dir, "03_component_agreement_conflict.csv"), row.names = FALSE)

# ------------------------------------------------------------
# Support-only mechanism diagnostics
# ------------------------------------------------------------
support_diag <- data.frame(
  metric = c(
    "mean_D1_size", "mean_D2_size",
    "mean_D1_signal", "mean_D2_signal",
    "mean_D1_null", "mean_D2_null",
    "prop_D1_signal_among_D1", "prop_D2_signal_among_D2"
  ),
  value = c(
    mean(valid$D1_size), mean(valid$D2_size),
    mean(valid$D1_signal), mean(valid$D2_signal),
    mean(valid$D1_null), mean(valid$D2_null),
    sum(valid$D1_signal) / max(sum(valid$D1_size), 1),
    sum(valid$D2_signal) / max(sum(valid$D2_size), 1)
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(support_diag, file.path(analysis_dir, "04_support_mechanism_diagnostics.csv"), row.names = FALSE)

# ------------------------------------------------------------
# Fed-FDR sanity checks
# ------------------------------------------------------------
fed_diag <- data.frame(
  metric = c(
    "mean_FDP", "mean_Power", "mean_final_selected_size",
    "mean_support1_size", "mean_support2_size",
    "mean_external1_size", "mean_external2_size",
    "hessian_regularization_rate_site1", "hessian_regularization_rate_site2",
    "mean_y_site1", "mean_y_site2"
  ),
  value = c(
    mean(valid$fed_fdp), mean(valid$fed_power), mean(valid$final_selected_size),
    mean(valid$support1_size), mean(valid$support2_size),
    mean(valid$external1_size), mean(valid$external2_size),
    mean(valid$hessian_regularized_site1), mean(valid$hessian_regularized_site2),
    mean(valid$y_mean_site1), mean(valid$y_mean_site2)
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(fed_diag, file.path(analysis_dir, "05_fedfdr_sanity_checks.csv"), row.names = FALSE)

# ------------------------------------------------------------
# Failure audit
# ------------------------------------------------------------
if (nrow(failed) > 0L) {
  fail_stage <- as.data.frame(sort(table(failed$failure_stage), decreasing = TRUE), stringsAsFactors = FALSE)
  names(fail_stage) <- c("failure_stage", "count")
  utils::write.csv(fail_stage, file.path(analysis_dir, "06_failure_stages.csv"), row.names = FALSE)

  fail_reason <- as.data.frame(sort(table(failed$failure_reason), decreasing = TRUE), stringsAsFactors = FALSE)
  names(fail_reason) <- c("failure_reason", "count")
  utils::write.csv(fail_reason, file.path(analysis_dir, "07_failure_reasons.csv"), row.names = FALSE)
}

# ------------------------------------------------------------
# Simple plots (base R only)
# ------------------------------------------------------------
pdf(file.path(analysis_dir, "attack_accuracy_plot.pdf"), width = 10, height = 6)
plot_df <- primary[primary$condition != "B0_ZERO_INFORMATION", , drop = FALSE]
ord <- seq_len(nrow(plot_df))
ylim <- range(c(0.45, plot_df$member_ci_lower, plot_df$member_ci_upper), na.rm = TRUE)
plot(ord, plot_df$member_accuracy,
     ylim = ylim, xaxt = "n", xlab = "Attack condition", ylab = "Source-site accuracy",
     pch = 19, main = "Fed-FDR Source-Site Attack Accuracy")
abline(h = 0.5, lty = 2)
arrows(ord, plot_df$member_ci_lower, ord, plot_df$member_ci_upper,
       angle = 90, code = 3, length = 0.04)
axis(1, at = ord, labels = plot_df$condition, las = 2)
dev.off()

pdf(file.path(analysis_dir, "member_vs_fresh_accuracy.pdf"), width = 10, height = 6)
mat <- rbind(plot_df$member_accuracy, plot_df$fresh_accuracy)
barplot(mat, beside = TRUE, names.arg = plot_df$condition, las = 2,
        ylab = "Accuracy", main = "Member vs Fresh Source-Site Accuracy",
        legend.text = c("Member", "Fresh"), args.legend = list(x = "topright"))
abline(h = 0.5, lty = 2)
dev.off()

cat("\nFull component decomposition:\n")
print(full_components, row.names = FALSE, digits = 4)
cat("\nAgreement/conflict diagnostics:\n")
print(relations, row.names = FALSE, digits = 4)
cat("\nSupport mechanism diagnostics:\n")
print(support_diag, row.names = FALSE, digits = 4)
cat("\nFed-FDR sanity checks:\n")
print(fed_diag, row.names = FALSE, digits = 4)

cat("\nAnalysis files saved to:\n", normalizePath(analysis_dir), "\n")
