# ============================================================
# Task B: population robustness of the frozen K=2 source-site attack
# Frozen Protocol v1.0 -- common functions
# ============================================================
# DESIGN PRINCIPLE:
#   ATTACK FIXED, POPULATION CHANGED.
#
# This file deliberately reuses the exact Fed-FDR fitting and R0--R3 attack
# implementation from Frozen source-site Protocol v1.1.1.  Only the X
# population generator is modularized.
#
# Main populations:
#   P0_TOEPLITZ_GAUSSIAN
#   P1_COPULA_UNIFORM
#   P2_BLOCK_SPARSE_GAUSSIAN
#   P3_TRIDIAG_NEARSING_GAUSSIAN
#   P4_POISSON_STANDARDIZED
#
# Optional sensitivity populations:
#   P1_IID_UNIFORM
#   P4_POISSON_RAW
#
# K=2, n1=n2=50, p=50, true support={1,...,5}; response model,
# Fed-FDR, candidate sampling, matched-fresh control, and attack scores are
# unchanged from Frozen v1.1.1.
# ============================================================

source(file.path("reference_old", "00_REFERENCE_FROZEN_K2_functions.R"))

TASKB_MAIN_POPULATIONS <- c(
  "P0_TOEPLITZ_GAUSSIAN",
  "P1_COPULA_UNIFORM",
  "P2_BLOCK_SPARSE_GAUSSIAN",
  "P3_TRIDIAG_NEARSING_GAUSSIAN",
  "P4_POISSON_STANDARDIZED"
)

TASKB_SENSITIVITY_POPULATIONS <- c(
  "P1_IID_UNIFORM",
  "P4_POISSON_RAW"
)

population_label <- function(population) {
  switch(
    population,
    P0_TOEPLITZ_GAUSSIAN = "P0 Toeplitz Gaussian",
    P1_COPULA_UNIFORM = "P1 Gaussian-copula correlated Uniform",
    P2_BLOCK_SPARSE_GAUSSIAN = "P2 fixed-permuted block-sparse Gaussian",
    P3_TRIDIAG_NEARSING_GAUSSIAN = "P3 tridiagonal near-singular Gaussian",
    P4_POISSON_STANDARDIZED = "P4 theoretically standardized spatial-Poisson grid counts",
    P1_IID_UNIFORM = "P1 sensitivity: iid Uniform",
    P4_POISSON_RAW = "P4 sensitivity: raw spatial-Poisson grid counts",
    population
  )
}

make_taskb_config <- function(
    population,
    n_reps = 5L,
    base_seed = 20260815L,
    diagnostic_n = 0L,
    p2_rho_range = c(0.75, 0.90),
    p3_rho_range = c(0.45, 0.49),
    p4_lambda_range = c(1, 4)
) {
  allowed <- c(TASKB_MAIN_POPULATIONS, TASKB_SENSITIVITY_POPULATIONS)
  if (!(population %in% allowed)) {
    stop("Unknown Task-B population: ", population)
  }

  cfg <- make_config(
    n_reps = as.integer(n_reps),
    base_seed = as.integer(base_seed),
    K = 2L,
    n_per_site = 50L,
    d = 50L,
    s = 5L,
    true_support = 1:5,
    alpha = 0.10,
    gamma_range = c(-0.5, 0.5),
    rho_range = c(0.3, 0.5),
    sigma_range = c(7, 11),
    nfolds = 5L,
    hessian_delta = 1e-5,
    score_eps = 1e-12,
    tie_tol = 1e-12
  )

  cfg$taskb_protocol_version <- "Task B Frozen Population Protocol v1.0"
  cfg$population <- population
  cfg$population_label <- population_label(population)
  cfg$diagnostic_n <- as.integer(diagnostic_n)
  cfg$p2_rho_range <- p2_rho_range
  cfg$p3_rho_range <- p3_rho_range
  cfg$p4_lambda_range <- p4_lambda_range

  # Fixed covariance permutation for P2.
  # Original variables 1:5 (the active support) are mapped to base block
  # positions 1,6,11,16,21, hence they lie in five different 5-variable blocks.
  # The remaining positions are a once-frozen permutation; the SAME placement
  # is used for both sites and every repetition.
  cfg$p2_base_position_for_variable <- c(
     1,  6, 11, 16, 21, 29, 30, 25, 28, 15,
    40, 20, 10,  9, 38, 48, 33, 27, 19, 35,
    12, 13, 45, 46, 14, 49,  2, 18, 36, 26,
    31, 17, 43, 23, 42, 44,  5, 32, 50,  4,
    47, 22, 37,  7, 24, 34, 41, 39,  3,  8
  )

  stopifnot(length(cfg$p2_base_position_for_variable) == cfg$d)
  stopifnot(length(unique(cfg$p2_base_position_for_variable)) == cfg$d)
  stopifnot(all(sort(cfg$p2_base_position_for_variable) == seq_len(cfg$d)))

  active_base_positions <- cfg$p2_base_position_for_variable[cfg$true_support]
  active_blocks <- ceiling(active_base_positions / 5)
  if (length(unique(active_blocks)) != cfg$s) {
    stop("P2 frozen permutation failed: active variables are not in distinct blocks.")
  }

  cfg
}

# ------------------------------------------------------------
# Population-specific parameter generation
# ------------------------------------------------------------

make_block_toeplitz_cov <- function(d, block_size, rho) {
  if (d %% block_size != 0L) stop("d must be divisible by block_size")
  B <- toeplitz(rho^(0:(block_size - 1L)))
  nb <- d / block_size
  out <- matrix(0, nrow = d, ncol = d)
  for (b in seq_len(nb)) {
    idx <- ((b - 1L) * block_size + 1L):(b * block_size)
    out[idx, idx] <- B
  }
  out
}

make_tridiag_cov <- function(d, rho) {
  S <- diag(1, d)
  if (d >= 2L) {
    idx <- seq_len(d - 1L)
    S[cbind(idx, idx + 1L)] <- rho
    S[cbind(idx + 1L, idx)] <- rho
  }
  S
}

covariance_metrics <- function(Sigma) {
  ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  lambda_min <- min(ev)
  lambda_max <- max(ev)
  kappa <- if (lambda_min > 0) lambda_max / lambda_min else Inf
  det_root <- if (all(ev > 0)) exp(mean(log(ev))) else NA_real_
  list(
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    kappa = kappa,
    det_root = det_root
  )
}

generate_parameters_taskb <- function(rep_id, cfg) {
  # P0 uses the exact old generator to make the baseline compatibility test
  # genuinely exact, not merely approximate.
  if (identical(cfg$population, "P0_TOEPLITZ_GAUSSIAN")) {
    old <- generate_parameters(rep_id, cfg)
    old$population <- cfg$population
    old$population_parameter <- old$rho
    old$population_parameter_name <- "rho"
    old$cov_metrics <- lapply(old$Sigma, covariance_metrics)
    return(old)
  }

  set.seed(seed_stream(cfg, rep_id, 1L))

  gamma <- stats::runif(cfg$K, cfg$gamma_range[1], cfg$gamma_range[2])

  if (cfg$population %in% c("P1_COPULA_UNIFORM", "P1_IID_UNIFORM")) {
    poppar <- stats::runif(cfg$K, cfg$rho_range[1], cfg$rho_range[2])
    poppar_name <- "rho"
  } else if (cfg$population == "P2_BLOCK_SPARSE_GAUSSIAN") {
    poppar <- stats::runif(cfg$K, cfg$p2_rho_range[1], cfg$p2_rho_range[2])
    poppar_name <- "rho"
  } else if (cfg$population == "P3_TRIDIAG_NEARSING_GAUSSIAN") {
    poppar <- stats::runif(cfg$K, cfg$p3_rho_range[1], cfg$p3_rho_range[2])
    poppar_name <- "rho"
  } else if (cfg$population %in% c("P4_POISSON_STANDARDIZED", "P4_POISSON_RAW")) {
    poppar <- stats::runif(cfg$K, cfg$p4_lambda_range[1], cfg$p4_lambda_range[2])
    poppar_name <- "lambda_cell"
  } else {
    stop("Unhandled population in parameter generation.")
  }

  signal_sign <- sample(c(-1, 1), cfg$s, replace = TRUE)
  sigma_active <- matrix(
    stats::runif(cfg$K * cfg$s, cfg$sigma_range[1], cfg$sigma_range[2]),
    nrow = cfg$K,
    ncol = cfg$s,
    byrow = TRUE
  )

  beta <- matrix(0, nrow = cfg$K, ncol = cfg$d)
  scale_factor <- sqrt(log(cfg$d) / cfg$n_per_site)
  for (k in seq_len(cfg$K)) {
    beta[k, cfg$true_support] <- signal_sign * sigma_active[k, ] * scale_factor
  }

  Sigma <- vector("list", cfg$K)
  cov_metrics <- vector("list", cfg$K)

  if (cfg$population == "P1_COPULA_UNIFORM") {
    for (k in seq_len(cfg$K)) {
      Sigma[[k]] <- toeplitz(poppar[k]^(0:(cfg$d - 1L)))
      cov_metrics[[k]] <- covariance_metrics(Sigma[[k]])
    }
  } else if (cfg$population == "P1_IID_UNIFORM") {
    for (k in seq_len(cfg$K)) {
      Sigma[[k]] <- diag(1, cfg$d)
      cov_metrics[[k]] <- covariance_metrics(Sigma[[k]])
    }
  } else if (cfg$population == "P2_BLOCK_SPARSE_GAUSSIAN") {
    pos <- cfg$p2_base_position_for_variable
    for (k in seq_len(cfg$K)) {
      Sbase <- make_block_toeplitz_cov(cfg$d, block_size = 5L, rho = poppar[k])
      # Sigma for the ORIGINAL variable labels 1,...,50.
      Sigma[[k]] <- Sbase[pos, pos, drop = FALSE]
      cov_metrics[[k]] <- covariance_metrics(Sigma[[k]])
    }
  } else if (cfg$population == "P3_TRIDIAG_NEARSING_GAUSSIAN") {
    for (k in seq_len(cfg$K)) {
      Sigma[[k]] <- make_tridiag_cov(cfg$d, poppar[k])
      cov_metrics[[k]] <- covariance_metrics(Sigma[[k]])
      if (!is.finite(cov_metrics[[k]]$lambda_min) || cov_metrics[[k]]$lambda_min <= 0) {
        stop("P3 covariance is not positive definite; lambda_min <= 0.")
      }
    }
  } else {
    # P4 has no Gaussian covariance matrix; disjoint equal-area PPP cell counts
    # are independent Poisson counts. lambda_cell is the mean count per cell.
    # On the unit square with 50 equal cells, the corresponding homogeneous
    # point-process intensity per unit area is 50 * lambda_cell.
    Sigma <- replicate(cfg$K, NULL, simplify = FALSE)
    cov_metrics <- replicate(cfg$K, list(
      lambda_min = NA_real_, lambda_max = NA_real_, kappa = NA_real_, det_root = NA_real_
    ), simplify = FALSE)
  }

  list(
    gamma = gamma,
    rho = if (poppar_name == "rho") poppar else rep(NA_real_, cfg$K),
    lambda_cell = if (poppar_name == "lambda_cell") poppar else rep(NA_real_, cfg$K),
    population_parameter = poppar,
    population_parameter_name = poppar_name,
    sigma_active = sigma_active,
    signal_sign = signal_sign,
    beta = beta,
    Sigma = Sigma,
    cov_metrics = cov_metrics,
    scale_factor = scale_factor,
    population = cfg$population
  )
}

# ------------------------------------------------------------
# X generation under each population
# ------------------------------------------------------------

generate_X_taskb <- function(n, site_k, params, cfg, seed_value) {
  set.seed(seed_value)

  pop <- cfg$population

  if (pop == "P0_TOEPLITZ_GAUSSIAN") {
    X <- mvtnorm::rmvnorm(n = n, mean = rep(0, cfg$d), sigma = params$Sigma[[site_k]])

  } else if (pop == "P1_COPULA_UNIFORM") {
    Z <- mvtnorm::rmvnorm(n = n, mean = rep(0, cfg$d), sigma = params$Sigma[[site_k]])
    U <- stats::pnorm(Z)
    X <- 2 * sqrt(3) * (U - 0.5)

  } else if (pop == "P1_IID_UNIFORM") {
    X <- matrix(
      stats::runif(n * cfg$d, min = -sqrt(3), max = sqrt(3)),
      nrow = n,
      ncol = cfg$d,
      byrow = TRUE
    )

  } else if (pop %in% c("P2_BLOCK_SPARSE_GAUSSIAN", "P3_TRIDIAG_NEARSING_GAUSSIAN")) {
    X <- mvtnorm::rmvnorm(n = n, mean = rep(0, cfg$d), sigma = params$Sigma[[site_k]])

  } else if (pop %in% c("P4_POISSON_STANDARDIZED", "P4_POISSON_RAW")) {
    lam <- params$lambda_cell[site_k]
    counts <- matrix(
      stats::rpois(n * cfg$d, lambda = lam),
      nrow = n,
      ncol = cfg$d,
      byrow = TRUE
    )
    if (pop == "P4_POISSON_STANDARDIZED") {
      # Theoretical standardization using the TRUE realized lambda_k.
      # This is performed by the data generator and does NOT release lambda_k
      # to the attacker.
      X <- (counts - lam) / sqrt(lam)
    } else {
      X <- counts
    }

  } else {
    stop("Unhandled Task-B population in X generation: ", pop)
  }

  X <- as.matrix(X)
  storage.mode(X) <- "double"
  X
}

generate_site_data_taskb <- function(rep_id, site_k, params, cfg) {
  # Exact P0 branch preserves the old baseline generator byte-for-byte in logic.
  if (identical(cfg$population, "P0_TOEPLITZ_GAUSSIAN")) {
    return(generate_site_data(rep_id, site_k, params, cfg))
  }

  X <- generate_X_taskb(
    n = cfg$n_per_site,
    site_k = site_k,
    params = params,
    cfg = cfg,
    seed_value = seed_stream(cfg, rep_id, 10L + as.integer(site_k))
  )

  eta <- as.vector(params$gamma[site_k] + X %*% params$beta[site_k, ])
  prob <- safe_plogis(eta)
  y <- stats::rbinom(cfg$n_per_site, size = 1L, prob = prob)

  list(X = X, y = as.integer(y), eta = eta, prob = prob)
}

generate_fresh_observation_taskb <- function(rep_id, source_site, params, cfg) {
  if (identical(cfg$population, "P0_TOEPLITZ_GAUSSIAN")) {
    return(generate_fresh_observation(rep_id, source_site, params, cfg))
  }

  X <- generate_X_taskb(
    n = 1L,
    site_k = source_site,
    params = params,
    cfg = cfg,
    seed_value = seed_stream(cfg, rep_id, 51L)
  )
  x <- as.numeric(X[1L, ])
  eta <- params$gamma[source_site] + sum(x * params$beta[source_site, ])
  p <- safe_plogis(eta)
  y <- stats::rbinom(1L, size = 1L, prob = p)
  list(x = x, y = as.integer(y), eta = eta, prob = p)
}

# ------------------------------------------------------------
# Oracle DGP diagnostics (Quick/Pilot only)
# ------------------------------------------------------------

auc_binary_rank_taskb <- function(score, y) {
  ok <- is.finite(score) & !is.na(y)
  s <- score[ok]
  yy <- as.integer(y[ok])
  n1 <- sum(yy == 1L)
  n0 <- sum(yy == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(s, ties.method = "average")
  (sum(r[yy == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

oracle_diagnostics_one_site <- function(rep_id, site_k, params, cfg) {
  M <- cfg$diagnostic_n
  if (M <= 0L) {
    return(list(
      diag_sd_xbeta = NA_real_,
      diag_sd_eta = NA_real_,
      diag_extreme_pi_prop = NA_real_,
      diag_mean_pi = NA_real_,
      diag_ybar = NA_real_,
      diag_oracle_auc = NA_real_
    ))
  }

  X <- generate_X_taskb(
    n = M,
    site_k = site_k,
    params = params,
    cfg = cfg,
    seed_value = seed_stream(cfg, rep_id, 70L + as.integer(site_k))
  )

  xbeta <- as.vector(X %*% params$beta[site_k, ])
  eta <- params$gamma[site_k] + xbeta
  pi <- safe_plogis(eta)

  # A separate deterministic Bernoulli stream avoids reusing training/fresh RNG.
  set.seed(seed_stream(cfg, rep_id, 80L + as.integer(site_k)))
  y <- stats::rbinom(M, size = 1L, prob = pi)

  list(
    diag_sd_xbeta = stats::sd(xbeta),
    diag_sd_eta = stats::sd(eta),
    diag_extreme_pi_prop = mean(pi < 0.05 | pi > 0.95),
    diag_mean_pi = mean(pi),
    diag_ybar = mean(y),
    diag_oracle_auc = auc_binary_rank_taskb(eta, y)
  )
}

# ------------------------------------------------------------
# One Task-B repetition
# ------------------------------------------------------------

run_one_rep_taskb <- function(rep_id, cfg, return_debug = FALSE) {
  stage <- "start"

  ans <- tryCatch({
    stage <- "parameter_generation"
    params <- generate_parameters_taskb(rep_id, cfg)

    stage <- "site_data_generation"
    site_data <- lapply(seq_len(cfg$K), function(k) {
      generate_site_data_taskb(rep_id, k, params, cfg)
    })

    stage <- "fedfdr_fit"
    # EXACT frozen fitting function.
    fedfit <- fit_fedfdr_one_rep(site_data, rep_id, cfg)

    stage <- "member_candidate"
    set.seed(seed_stream(cfg, rep_id, 41L))
    source_site <- sample(seq_len(cfg$K), size = 1L)
    member_index <- sample(seq_len(cfg$n_per_site), size = 1L)
    member <- list(
      x = as.numeric(site_data[[source_site]]$X[member_index, ]),
      y = as.integer(site_data[[source_site]]$y[member_index])
    )

    stage <- "fresh_candidate"
    fresh <- generate_fresh_observation_taskb(rep_id, source_site, params, cfg)

    stage <- "member_attack"
    member_attack <- attack_one_candidate(
      x = member$x, y = member$y, source_site = source_site,
      fedfit = fedfit, cfg = cfg
    )

    stage <- "fresh_attack"
    fresh_attack <- attack_one_candidate(
      x = fresh$x, y = fresh$y, source_site = source_site,
      fedfit = fedfit, cfg = cfg
    )

    stage <- "oracle_diagnostics"
    diag1 <- oracle_diagnostics_one_site(rep_id, 1L, params, cfg)
    diag2 <- oracle_diagnostics_one_site(rep_id, 2L, params, cfg)

    stage <- "diagnostics"
    fed_diag <- calc_fdp_power(
      selection01 = fedfit$central$final_selection,
      true_support = cfg$true_support,
      d = cfg$d
    )

    support1 <- fedfit$support_local[[1L]]
    support2 <- fedfit$support_local[[2L]]
    ext1 <- fedfit$second[[1L]]$external_support
    ext2 <- fedfit$second[[2L]]$external_support

    cm1 <- params$cov_metrics[[1L]]
    cm2 <- params$cov_metrics[[2L]]

    base <- list(
      population = cfg$population,
      population_label = cfg$population_label,
      rep_id = as.integer(rep_id),
      seed = seed_stream(cfg, rep_id, 0L),
      status = "ok",
      failure_stage = "",
      failure_reason = "",
      source_site = as.integer(source_site),
      member_index = as.integer(member_index),
      member_y = as.integer(member$y),
      fresh_y = as.integer(fresh$y),
      y_mean_site1 = mean(site_data[[1L]]$y),
      y_mean_site2 = mean(site_data[[2L]]$y),
      gamma1 = params$gamma[1L],
      gamma2 = params$gamma[2L],
      population_parameter1 = params$population_parameter[1L],
      population_parameter2 = params$population_parameter[2L],
      population_parameter_name = params$population_parameter_name,
      rho1 = params$rho[1L] %||% NA_real_,
      rho2 = params$rho[2L] %||% NA_real_,
      lambda_cell1 = params$lambda_cell[1L] %||% NA_real_,
      lambda_cell2 = params$lambda_cell[2L] %||% NA_real_,
      p4_process_intensity1 = if (is.finite(params$lambda_cell[1L] %||% NA_real_)) 50 * params$lambda_cell[1L] else NA_real_,
      p4_process_intensity2 = if (is.finite(params$lambda_cell[2L] %||% NA_real_)) 50 * params$lambda_cell[2L] else NA_real_,
      cov_lambda_min_site1 = cm1$lambda_min,
      cov_lambda_min_site2 = cm2$lambda_min,
      cov_lambda_max_site1 = cm1$lambda_max,
      cov_lambda_max_site2 = cm2$lambda_max,
      cov_kappa_site1 = cm1$kappa,
      cov_kappa_site2 = cm2$kappa,
      cov_det_root_site1 = cm1$det_root,
      cov_det_root_site2 = cm2$det_root,
      diag_sd_xbeta_site1 = diag1$diag_sd_xbeta,
      diag_sd_xbeta_site2 = diag2$diag_sd_xbeta,
      diag_sd_eta_site1 = diag1$diag_sd_eta,
      diag_sd_eta_site2 = diag2$diag_sd_eta,
      diag_extreme_pi_prop_site1 = diag1$diag_extreme_pi_prop,
      diag_extreme_pi_prop_site2 = diag2$diag_extreme_pi_prop,
      diag_mean_pi_site1 = diag1$diag_mean_pi,
      diag_mean_pi_site2 = diag2$diag_mean_pi,
      diag_ybar_site1 = diag1$diag_ybar,
      diag_ybar_site2 = diag2$diag_ybar,
      diag_oracle_auc_site1 = diag1$diag_oracle_auc,
      diag_oracle_auc_site2 = diag2$diag_oracle_auc,
      signal_scale = params$scale_factor,
      sign_1 = params$signal_sign[1L],
      sign_2 = params$signal_sign[2L],
      sign_3 = params$signal_sign[3L],
      sign_4 = params$signal_sign[4L],
      sign_5 = params$signal_sign[5L],
      sigma11 = params$sigma_active[1L, 1L],
      sigma12 = params$sigma_active[1L, 2L],
      sigma13 = params$sigma_active[1L, 3L],
      sigma14 = params$sigma_active[1L, 4L],
      sigma15 = params$sigma_active[1L, 5L],
      sigma21 = params$sigma_active[2L, 1L],
      sigma22 = params$sigma_active[2L, 2L],
      sigma23 = params$sigma_active[2L, 3L],
      sigma24 = params$sigma_active[2L, 4L],
      sigma25 = params$sigma_active[2L, 5L],
      support1_size = length(support1),
      support2_size = length(support2),
      external1_size = length(ext1),
      external2_size = length(ext2),
      support1 = indices_to_string(support1),
      support2 = indices_to_string(support2),
      external1 = indices_to_string(ext1),
      external2 = indices_to_string(ext2),
      final_support = indices_to_string(which(fedfit$central$final_selection != 0)),
      final_selected_size = fed_diag$selected_size,
      fed_fdp = fed_diag$fdp,
      fed_power = fed_diag$power,
      pair_tau = fedfit$central$pair_tau[1L] %||% NA_real_,
      hessian_regularized_site1 = as.integer(fedfit$second[[1L]]$hessian_regularized),
      hessian_regularized_site2 = as.integer(fedfit$second[[2L]]$hessian_regularized),
      singleton_glmnet_patch_site1 = as.integer(fedfit$second[[1L]]$singleton_glmnet_patch),
      singleton_glmnet_patch_site2 = as.integer(fedfit$second[[2L]]$singleton_glmnet_patch),
      D1_signal = member_attack$support$D1_signal,
      D1_null = member_attack$support$D1_null,
      D2_signal = member_attack$support$D2_signal,
      D2_null = member_attack$support$D2_null,
      D1_size = length(member_attack$support$D1),
      D2_size = length(member_attack$support$D2),
      B0_member_credit = 0.5,
      B0_fresh_credit = 0.5
    )

    row <- c(
      base,
      flatten_attack("member", member_attack),
      flatten_attack("fresh", fresh_attack)
    )

    if (return_debug) {
      list(summary = row, debug = list(
        cfg = cfg, params = params, site_data = site_data, fedfit = fedfit,
        member = member, fresh = fresh,
        member_attack = member_attack, fresh_attack = fresh_attack,
        oracle_diag_site1 = diag1, oracle_diag_site2 = diag2
      ))
    } else {
      row
    }
  }, error = function(e) {
    fail <- list(
      population = cfg$population,
      population_label = cfg$population_label,
      rep_id = as.integer(rep_id),
      seed = seed_stream(cfg, rep_id, 0L),
      status = "failed",
      failure_stage = stage,
      failure_reason = conditionMessage(e)
    )
    if (return_debug) list(summary = fail, debug = NULL) else fail
  })

  ans
}

# ------------------------------------------------------------
# Task-B summary helpers
# ------------------------------------------------------------

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}

summarize_taskb_dgp <- function(df) {
  ok <- df[df$status == "ok", , drop = FALSE]
  if (nrow(ok) == 0L) return(data.frame())

  data.frame(
    population = ok$population[1L],
    population_label = ok$population_label[1L],
    n_valid = nrow(ok),
    poppar1_mean = safe_mean(ok$population_parameter1),
    poppar2_mean = safe_mean(ok$population_parameter2),
    train_ybar_site1 = safe_mean(ok$y_mean_site1),
    train_ybar_site2 = safe_mean(ok$y_mean_site2),
    oracle_sd_xbeta_site1 = safe_mean(ok$diag_sd_xbeta_site1),
    oracle_sd_xbeta_site2 = safe_mean(ok$diag_sd_xbeta_site2),
    oracle_extreme_pi_site1 = safe_mean(ok$diag_extreme_pi_prop_site1),
    oracle_extreme_pi_site2 = safe_mean(ok$diag_extreme_pi_prop_site2),
    oracle_mean_pi_site1 = safe_mean(ok$diag_mean_pi_site1),
    oracle_mean_pi_site2 = safe_mean(ok$diag_mean_pi_site2),
    oracle_auc_site1 = safe_mean(ok$diag_oracle_auc_site1),
    oracle_auc_site2 = safe_mean(ok$diag_oracle_auc_site2),
    cov_lambda_min_site1 = safe_mean(ok$cov_lambda_min_site1),
    cov_lambda_min_site2 = safe_mean(ok$cov_lambda_min_site2),
    cov_kappa_site1 = safe_mean(ok$cov_kappa_site1),
    cov_kappa_site2 = safe_mean(ok$cov_kappa_site2),
    cov_det_root_site1 = safe_mean(ok$cov_det_root_site1),
    cov_det_root_site2 = safe_mean(ok$cov_det_root_site2),
    stringsAsFactors = FALSE
  )
}

summarize_taskb_attacks <- function(df) {
  ok <- df[df$status == "ok", , drop = FALSE]
  if (nrow(ok) == 0L) return(data.frame())
  conds <- c("R0_X", "R0_XY", "R1_X", "R1_XY", "R2_X", "R2_XY", "R3_X", "R3_XY")
  out <- do.call(rbind, lapply(conds, function(cc) summarize_condition(ok, cc)))
  out$population <- ok$population[1L]
  out$population_label <- ok$population_label[1L]
  out <- out[, c("population", "population_label", setdiff(names(out), c("population", "population_label")))]
  out
}

write_taskb_p2_permutation <- function(cfg, file) {
  pos <- cfg$p2_base_position_for_variable
  tab <- data.frame(
    variable = seq_len(cfg$d),
    is_active = seq_len(cfg$d) %in% cfg$true_support,
    base_block_position = pos,
    block_id = ceiling(pos / 5),
    within_block_position = ((pos - 1L) %% 5L) + 1L
  )
  utils::write.csv(tab, file, row.names = FALSE)
  invisible(tab)
}
