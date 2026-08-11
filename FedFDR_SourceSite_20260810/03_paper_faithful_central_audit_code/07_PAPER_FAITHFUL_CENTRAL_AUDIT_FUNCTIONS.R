# ============================================================
# Paper-faithful Fed-FDR central-step audit
# ============================================================
# Purpose
# -------
# Keep ALL upstream steps identical to Frozen Protocol v1.1:
#   * same parameter draws and deterministic seed streams
#   * same K=2, n1=n2=50, d=50 data
#   * same first-stage ordinary GLM Lasso (5-fold CV)
#   * same between-site support swap
#   * same second-stage lower-dimensional GLM Lasso + refined debiasing
#
# Change ONLY the central Fed-FDR selection step and compare three variants:
#
# A. CODE_THETA
#    Current formal experiment / public-repository-compatible implementation:
#       mirror input = theta_mirror = beta_refined / sqrt(SE)
#       central threshold + inclusion aggregation = existing mfdr_public_logic()
#
# B. CODE_RAWBETA
#    Input-only correction:
#       mirror input = raw zero-padded refined beta
#       central threshold + inclusion aggregation = existing mfdr_public_logic()
#    This isolates the effect of replacing theta_mirror by the paper's beta-hat.
#
# C. PAPER_RAWBETA
#    Paper-text central step:
#       mirror input = raw zero-padded refined beta
#       pairwise threshold uses the strict inequalities printed in Algorithm 1
#       final support uses {j: I_j > I_(m)} exactly as printed in Algorithm 1
#
# IMPORTANT K=2 NOTE
# ------------------
# With K=2 there is only one site pair. Hence all features selected by that pair
# have the same inclusion rate 1/|S^(12)|. The final Algorithm-1 inclusion step
# can therefore be degenerate because of ties. We record BOTH the pairwise
# selection and the final inclusion-aggregated support.
# ============================================================

source("01_FedFDR_SourceAttack_functions.R")

# ------------------------------------------------------------
# Paper Algorithm 1: pairwise mirror threshold
# ------------------------------------------------------------

fdr_control_pair_paper <- function(beta1, beta2, q) {
  M <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))

  # The paper writes
  #   tau = min{tau > 0 : #{M < -tau}/(#{M > tau} v 1) <= q}
  # and S = {j: M_j > tau}.
  # The selection only changes when tau crosses an observed |M_j|, so the
  # positive observed magnitudes are used as the finite candidate grid.
  Tpos <- sort(unique(abs(M[M != 0])))

  if (length(Tpos) == 0L) {
    return(list(
      selection = integer(length(M)),
      mirror = M,
      tau_hat = Inf,
      fdp_hat_at_tau = 0,
      tau_convention = "all_mirror_zero"
    ))
  }

  # Finite-sample boundary convention for the paper's displayed rule.
  # The paper writes a minimum over tau > 0 with strict inequalities. If the
  # estimated FDP is already <= q for arbitrarily small positive tau, the set
  # has infimum 0 but may have no strictly positive minimum. This can occur not
  # only when there are no negative mirror statistics, but more generally when
  #   #{M < 0} / (#{M > 0} v 1) <= q.
  # We therefore use tau=0 as the explicitly recorded boundary-limit convention
  # whenever the right-limit at 0 is already feasible. Otherwise, we search the
  # positive observed |M_j| values using the paper's strict inequalities.
  fdp_zero_plus <- sum(M < 0) / max(sum(M > 0), 1L)

  if (fdp_zero_plus <= q) {
    tau_hat <- 0
    fdp_hat <- fdp_zero_plus
    convention <- "tau0_boundary_infimum_feasible"
  } else {
    tau_hat <- Inf
    fdp_hat <- NA_real_
    convention <- "paper_positive_observed_threshold_grid"

    for (tau in Tpos) {
      cur <- sum(M < (-tau)) / max(sum(M > tau), 1L)
      if (cur <= q) {
        tau_hat <- tau
        fdp_hat <- cur
        break
      }
    }
  }

  if (is.infinite(tau_hat)) {
    S <- integer(length(M))
  } else {
    S <- as.integer(M > tau_hat)
  }

  list(
    selection = S,
    mirror = M,
    tau_hat = tau_hat,
    fdp_hat_at_tau = fdp_hat,
    tau_convention = convention
  )
}

# ------------------------------------------------------------
# Paper Algorithm 1: exact inclusion-rate aggregation
# ------------------------------------------------------------

mfdr_paper_algorithm1 <- function(Beta, q) {
  K <- nrow(Beta)
  p <- ncol(Beta)
  pair_selections <- NULL
  pair_taus <- numeric(0)
  pair_fdp_hats <- numeric(0)
  pair_conventions <- character(0)

  if (K < 2L) stop("Algorithm 1 requires K >= 2.")

  for (s in seq_len(K - 1L)) {
    for (t in (s + 1L):K) {
      fit <- fdr_control_pair_paper(Beta[s, ], Beta[t, ], q)
      pair_selections <- rbind(pair_selections, fit$selection)
      pair_taus <- c(pair_taus, fit$tau_hat)
      pair_fdp_hats <- c(pair_fdp_hats, fit$fdp_hat_at_tau)
      pair_conventions <- c(pair_conventions, fit$tau_convention)
    }
  }

  n_pairs <- K * (K - 1) / 2
  pair_sizes <- rowSums(pair_selections)
  denom <- pmax(pair_sizes, 1L)

  # Algorithm 1, line 7:
  # I_j = 2/[K(K-1)] sum_{s<t} I(j in S^(st))/(|S^(st)| v 1)
  inclusion <- colSums(pair_selections / denom) / n_pairs

  # Algorithm 1, lines 8-10.
  ord <- order(inclusion)
  I_sorted <- inclusion[ord]
  cum_I <- cumsum(I_sorted)
  valid_m <- which(cum_I <= q)

  if (length(valid_m) == 0L) {
    # The printed algorithm does not specify this finite-sample edge case.
    # With d=50 and q=.1 this is not expected in our K=2 experiment because
    # zero-inclusion coordinates or sufficiently small 1/p values normally
    # make at least one m feasible. We nevertheless record it transparently.
    m <- 0L
    threshold_I <- -Inf
    final <- rep(1L, p)
    m_convention <- "no_m_exists_select_all_by_threshold_minus_inf"
  } else {
    m <- max(valid_m)
    threshold_I <- I_sorted[m]
    # EXACT paper line 10: S = {j: I_j > I_(m)}
    final <- as.integer(inclusion > threshold_I)
    m_convention <- "paper_exact_strict_greater_than_I_m"
  }

  list(
    final_selection = final,
    inclusion = inclusion,
    pair_selection = pair_selections,
    pair_tau = pair_taus,
    pair_fdp_hat = pair_fdp_hats,
    pair_tau_convention = pair_conventions,
    inclusion_order = ord,
    inclusion_sorted = I_sorted,
    inclusion_cumsum = cum_I,
    m = m,
    inclusion_threshold = threshold_I,
    m_convention = m_convention
  )
}

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

selection_metrics <- function(sel01, true_support) {
  selected <- which(sel01 != 0)
  tp <- sum(selected %in% true_support)
  fp <- sum(!(selected %in% true_support))
  list(
    size = length(selected),
    fdp = fp / max(length(selected), 1L),
    power = tp / length(true_support),
    support = selected
  )
}

jaccard_index <- function(a, b) {
  a <- as.integer(a)
  b <- as.integer(b)
  U <- union(a, b)
  if (length(U) == 0L) return(1)
  length(intersect(a, b)) / length(U)
}

as_support_string <- function(sel01) {
  idx <- which(sel01 != 0)
  if (length(idx) == 0L) return("")
  paste(idx, collapse = ",")
}

# ------------------------------------------------------------
# One audit repetition
# ------------------------------------------------------------

run_central_audit_rep <- function(rep_id, cfg) {
  stage <- "start"

  tryCatch({
    stage <- "parameter_generation"
    params <- generate_parameters(rep_id, cfg)

    stage <- "site_data_generation"
    site_data <- lapply(seq_len(cfg$K), function(k) {
      generate_site_data(rep_id, k, params, cfg)
    })

    stage <- "upstream_fedfdr_fit"
    # This reproduces the exact same upstream fit used in the formal privacy
    # experiment and also gives us the original CODE_THETA central result.
    fedfit <- fit_fedfdr_one_rep(site_data, rep_id, cfg)

    stage <- "central_variants"
    # A: original formal/public-code-compatible central result
    central_code_theta <- fedfit$central

    # B: isolate mirror-input difference only
    central_code_rawbeta <- mfdr_public_logic(fedfit$beta_refined, cfg$alpha)

    # C: paper-text Algorithm 1 with raw refined beta
    central_paper_rawbeta <- mfdr_paper_algorithm1(fedfit$beta_refined, cfg$alpha)

    stage <- "diagnostics"
    A_final <- selection_metrics(central_code_theta$final_selection, cfg$true_support)
    B_final <- selection_metrics(central_code_rawbeta$final_selection, cfg$true_support)
    C_final <- selection_metrics(central_paper_rawbeta$final_selection, cfg$true_support)

    # With K=2 there is exactly one pair, so pairwise diagnostics are direct.
    A_pair_sel <- as.integer(central_code_theta$pair_selection[1L, ])
    B_pair_sel <- as.integer(central_code_rawbeta$pair_selection[1L, ])
    C_pair_sel <- as.integer(central_paper_rawbeta$pair_selection[1L, ])

    A_pair <- selection_metrics(A_pair_sel, cfg$true_support)
    B_pair <- selection_metrics(B_pair_sel, cfg$true_support)
    C_pair <- selection_metrics(C_pair_sel, cfg$true_support)

    data.frame(
      rep_id = as.integer(rep_id),
      status = "ok",
      failure_stage = "",
      failure_reason = "",

      support1_size = length(fedfit$support_local[[1L]]),
      support2_size = length(fedfit$support_local[[2L]]),
      external1_size = length(fedfit$second[[1L]]$external_support),
      external2_size = length(fedfit$second[[2L]]$external_support),

      # A. Original code-compatible Theta path
      A_code_theta_pair_tau = central_code_theta$pair_tau[1L],
      A_code_theta_pair_size = A_pair$size,
      A_code_theta_pair_fdp = A_pair$fdp,
      A_code_theta_pair_power = A_pair$power,
      A_code_theta_final_size = A_final$size,
      A_code_theta_final_fdp = A_final$fdp,
      A_code_theta_final_power = A_final$power,
      A_code_theta_final_support = as_support_string(central_code_theta$final_selection),

      # B. Same public-code central logic, raw beta input
      B_code_rawbeta_pair_tau = central_code_rawbeta$pair_tau[1L],
      B_code_rawbeta_pair_size = B_pair$size,
      B_code_rawbeta_pair_fdp = B_pair$fdp,
      B_code_rawbeta_pair_power = B_pair$power,
      B_code_rawbeta_final_size = B_final$size,
      B_code_rawbeta_final_fdp = B_final$fdp,
      B_code_rawbeta_final_power = B_final$power,
      B_code_rawbeta_final_support = as_support_string(central_code_rawbeta$final_selection),

      # C. Paper-text Algorithm 1, raw beta input
      C_paper_rawbeta_pair_tau = central_paper_rawbeta$pair_tau[1L],
      C_paper_rawbeta_pair_fdp_hat = central_paper_rawbeta$pair_fdp_hat[1L],
      C_paper_rawbeta_tau_convention = central_paper_rawbeta$pair_tau_convention[1L],
      C_paper_rawbeta_pair_size = C_pair$size,
      C_paper_rawbeta_pair_fdp = C_pair$fdp,
      C_paper_rawbeta_pair_power = C_pair$power,
      C_paper_rawbeta_inclusion_m = central_paper_rawbeta$m,
      C_paper_rawbeta_inclusion_threshold = central_paper_rawbeta$inclusion_threshold,
      C_paper_rawbeta_m_convention = central_paper_rawbeta$m_convention,
      C_paper_rawbeta_final_size = C_final$size,
      C_paper_rawbeta_final_fdp = C_final$fdp,
      C_paper_rawbeta_final_power = C_final$power,
      C_paper_rawbeta_final_support = as_support_string(central_paper_rawbeta$final_selection),

      # Decomposition comparisons
      jaccard_A_vs_B_pair = jaccard_index(A_pair$support, B_pair$support),
      jaccard_B_vs_C_pair = jaccard_index(B_pair$support, C_pair$support),
      jaccard_A_vs_C_pair = jaccard_index(A_pair$support, C_pair$support),
      jaccard_A_vs_B_final = jaccard_index(A_final$support, B_final$support),
      jaccard_B_vs_C_final = jaccard_index(B_final$support, C_final$support),
      jaccard_A_vs_C_final = jaccard_index(A_final$support, C_final$support),
      same_A_B_pair = as.integer(identical(A_pair_sel, B_pair_sel)),
      same_B_C_pair = as.integer(identical(B_pair_sel, C_pair_sel)),
      same_A_C_pair = as.integer(identical(A_pair_sel, C_pair_sel)),
      same_A_B_final = as.integer(identical(as.integer(central_code_theta$final_selection), as.integer(central_code_rawbeta$final_selection))),
      same_B_C_final = as.integer(identical(as.integer(central_code_rawbeta$final_selection), as.integer(central_paper_rawbeta$final_selection))),
      same_A_C_final = as.integer(identical(as.integer(central_code_theta$final_selection), as.integer(central_paper_rawbeta$final_selection))),

      hessian_regularized_site1 = as.integer(fedfit$second[[1L]]$hessian_regularized),
      hessian_regularized_site2 = as.integer(fedfit$second[[2L]]$hessian_regularized),
      singleton_patch_site1 = as.integer(fedfit$second[[1L]]$singleton_glmnet_patch),
      singleton_patch_site2 = as.integer(fedfit$second[[2L]]$singleton_glmnet_patch),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      rep_id = as.integer(rep_id),
      status = "failed",
      failure_stage = stage,
      failure_reason = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
}

# Bind rows with different failure/success columns safely.
audit_rows_to_df <- function(rows) {
  nms <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows2 <- lapply(rows, function(x) {
    miss <- setdiff(nms, names(x))
    for (nm in miss) x[[nm]] <- NA
    x[nms]
  })
  out <- do.call(rbind, rows2)
  rownames(out) <- NULL
  out
}

fmt_sec <- function(x) {
  if (!is.finite(x)) return("NA")
  x <- max(0, x)
  h <- floor(x / 3600)
  m <- floor((x %% 3600) / 60)
  s <- round(x %% 60)
  sprintf("%02dh:%02dm:%02ds", h, m, s)
}
