# ============================================================
# Task A: Source-Site Attack Sensitivity to the Number of Sites
# K-General Source Attack v1.0 -- common functions
# ============================================================
# Purpose:
#   Extend the frozen K=2 source-site privacy attacks to K >= 2 while
#   preserving the original K=2 decision rules exactly.
#
# Frozen statistical choices:
#   - Known training member; infer source site J in {1,...,K}.
#   - J is sampled uniformly, so the random baseline is 1/K.
#   - Candidate information: X-only or (X,Y).
#   - Information regimes:
#       R0 = final site-agnostic output baseline;
#       R1 = site-specific first-stage supports;
#       R2 = site-specific refined coefficient objects;
#       R3 = supports + refined coefficients.
#   - R3 uses fixed equal weight; no formal-data tuning.
#   - Exact score ties receive expected uniform tie credit.
#   - Intercepts are NOT released to the attacker.
#   - Privacy attacks use RAW refined slope coefficients, not Theta.
#
# Task-A default experiment:
#   K in {2,3,5,10,20,50,100}; n_k=50; p=50; p1=5; alpha=0.1.
#   Only K changes in Task A. The baseline Gaussian-Toeplitz DGP is retained.
#
# Backward compatibility:
#   For K=2, RNG streams preserve Frozen Protocol v1.1.1 and the new
#   multi-site scores algebraically reduce to the old D-scores.
# ============================================================

required_pkgs <- c("glmnet", "mvtnorm")

check_required_packages <- function() {
  missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Missing R package(s): ", paste(missing, collapse = ", "),
      ". Install them before running the experiment."
    )
  }
  invisible(TRUE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

make_config <- function(
    n_reps = 5L,
    base_seed = 20260807L,
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
    tie_tol = 1e-12,
    compute_central = TRUE
) {
  stopifnot(K >= 2L)
  stopifnot(n_per_site >= 2L)
  stopifnot(d >= 2L)
  stopifnot(length(true_support) == s)
  stopifnot(all(true_support >= 1L & true_support <= d))

  list(
    protocol_version = "Task A Frozen Attack Definition v1.0",
    implementation_version = "K-General Source Attack v1.0",
    n_reps = as.integer(n_reps),
    base_seed = as.integer(base_seed),
    K = as.integer(K),
    n_per_site = as.integer(n_per_site),
    d = as.integer(d),
    s = as.integer(s),
    true_support = as.integer(true_support),
    alpha = alpha,
    gamma_range = gamma_range,
    rho_range = rho_range,
    sigma_range = sigma_range,
    nfolds = as.integer(nfolds),
    hessian_delta = hessian_delta,
    score_eps = score_eps,
    tie_tol = tie_tol,
    compute_central = isTRUE(compute_central)
  )
}

# ------------------------------------------------------------
# Deterministic RNG streams
# ------------------------------------------------------------

seed_stream_v11 <- function(cfg, rep_id, stream_id) {
  # Exact legacy formula from Frozen Protocol v1.1.1.
  z <- as.double(cfg$base_seed) + as.double(rep_id) * 1000 + as.double(stream_id) * 37
  z <- z %% (.Machine$integer.max - 1)
  as.integer(z + 1)
}

seed_for <- function(cfg, rep_id, phase, site_k = 0L) {
  # For K=2, preserve the exact old streams so a same-repetition regression
  # test can reproduce the Frozen v1.1.1 run.
  if (cfg$K == 2L) {
    stream_id <- switch(
      phase,
      rep = 0L,
      parameters = 1L,
      site_data = 10L + as.integer(site_k),
      first_stage = 20L + as.integer(site_k),
      second_stage = 30L + as.integer(site_k),
      member = 41L,
      fresh = 51L,
      stop("Unknown RNG phase: ", phase)
    )
    return(seed_stream_v11(cfg, rep_id, stream_id))
  }

  # For K>2, separate phase and site identifiers to avoid stream collisions.
  phase_code <- switch(
    phase,
    rep = 0L,
    parameters = 1L,
    site_data = 2L,
    first_stage = 3L,
    second_stage = 4L,
    member = 5L,
    fresh = 6L,
    stop("Unknown RNG phase: ", phase)
  )

  z <- as.double(cfg$base_seed) +
    as.double(rep_id) * 100000 +
    as.double(phase_code) * 1000 +
    as.double(site_k)
  z <- z %% (.Machine$integer.max - 1)
  as.integer(z + 1)
}

safe_plogis <- function(x) stats::plogis(x)

# ------------------------------------------------------------
# Data generation -- unchanged baseline DGP except K may vary
# ------------------------------------------------------------

generate_parameters <- function(rep_id, cfg) {
  set.seed(seed_for(cfg, rep_id, "parameters"))

  gamma <- stats::runif(cfg$K, cfg$gamma_range[1], cfg$gamma_range[2])
  rho <- stats::runif(cfg$K, cfg$rho_range[1], cfg$rho_range[2])

  # Common active-feature direction across all sites; site-specific magnitudes.
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

  Sigma <- lapply(seq_len(cfg$K), function(k) {
    toeplitz(rho[k]^(0:(cfg$d - 1L)))
  })

  list(
    gamma = gamma,
    rho = rho,
    sigma_active = sigma_active,
    signal_sign = signal_sign,
    beta = beta,
    Sigma = Sigma,
    scale_factor = scale_factor
  )
}

generate_site_data <- function(rep_id, site_k, params, cfg) {
  set.seed(seed_for(cfg, rep_id, "site_data", site_k))

  X <- mvtnorm::rmvnorm(
    n = cfg$n_per_site,
    mean = rep(0, cfg$d),
    sigma = params$Sigma[[site_k]]
  )
  X <- as.matrix(X)
  storage.mode(X) <- "double"

  eta <- as.vector(params$gamma[site_k] + X %*% params$beta[site_k, ])
  prob <- safe_plogis(eta)
  y <- stats::rbinom(cfg$n_per_site, size = 1L, prob = prob)

  list(X = X, y = as.integer(y), eta = eta, prob = prob)
}

generate_fresh_observation <- function(rep_id, source_site, params, cfg) {
  # Independent draw from the SAME realized source-site DGP as the member.
  set.seed(seed_for(cfg, rep_id, "fresh"))

  x <- as.numeric(mvtnorm::rmvnorm(
    n = 1L,
    mean = rep(0, cfg$d),
    sigma = params$Sigma[[source_site]]
  ))
  eta <- params$gamma[source_site] + sum(x * params$beta[source_site, ])
  p <- safe_plogis(eta)
  y <- stats::rbinom(1L, size = 1L, prob = p)

  list(x = x, y = as.integer(y), eta = eta, prob = p)
}

# ------------------------------------------------------------
# Fed-FDR fitting helpers -- preserved from Frozen v1.1.1
# ------------------------------------------------------------

validate_binary_outcome <- function(y, label = "outcome") {
  if (length(unique(y)) < 2L) {
    stop(label, " contains only one class.")
  }
  invisible(TRUE)
}

fit_first_stage_lasso <- function(X, y, rep_id, site_k, cfg) {
  validate_binary_outcome(y, paste0("site ", site_k, " first-stage outcome"))
  set.seed(seed_for(cfg, rep_id, "first_stage", site_k))

  cvfit <- glmnet::cv.glmnet(
    x = X,
    y = y,
    family = "binomial",
    nfolds = cfg$nfolds,
    alpha = 1
  )

  lambda <- cvfit$lambda.min
  fit <- glmnet::glmnet(
    x = X,
    y = y,
    family = "binomial",
    lambda = lambda,
    alpha = 1
  )

  coef_all <- as.vector(stats::coef(fit))
  slope <- coef_all[-1L]
  support <- which(slope != 0)

  list(
    support = as.integer(support),
    lambda = as.numeric(lambda),
    lasso_intercept = coef_all[1L],
    lasso_slope = slope
  )
}

refined_debiased_lasso <- function(X, y, lasso_est, delta = 1e-5) {
  nn <- length(y)
  X_aug <- cbind(1, X)
  if (ncol(X_aug) != length(lasso_est)) {
    stop("refined_debiased_lasso: lasso_est length does not match design matrix.")
  }

  eta <- as.vector(X_aug %*% lasso_est)
  mu <- safe_plogis(eta)

  grad <- -as.vector(crossprod(X_aug, y - mu)) / nn
  W <- mu * (1 - mu)
  hess <- crossprod(X_aug, W * X_aug) / nn

  hessian_regularized <- FALSE
  theta_inv <- tryCatch(
    solve(hess),
    error = function(e) NULL
  )

  if (is.null(theta_inv) || any(!is.finite(theta_inv))) {
    hessian_regularized <- TRUE
    theta_inv <- tryCatch(
      solve(hess + delta * diag(ncol(X_aug))),
      error = function(e) NULL
    )
  }

  if (is.null(theta_inv) || any(!is.finite(theta_inv))) {
    stop("Refined debiasing Hessian could not be inverted, even after delta regularization.")
  }

  diag_theta <- diag(theta_inv)
  if (any(!is.finite(diag_theta)) || any(diag_theta <= 0)) {
    stop("Refined debiasing produced non-positive/invalid inverse-Hessian diagonal.")
  }

  est <- as.vector(lasso_est - theta_inv %*% grad)
  se <- sqrt(diag_theta) / sqrt(nn)

  if (any(!is.finite(est)) || any(!is.finite(se)) || any(se <= 0)) {
    stop("Refined debiasing produced invalid estimates or standard errors.")
  }

  list(
    est = est,
    se = se,
    theta_inv = theta_inv,
    hessian_regularized = hessian_regularized
  )
}

fit_second_stage_refined <- function(X, y, external_support, rep_id, site_k, cfg) {
  if (length(external_support) == 0L) {
    stop("External support is empty; second-stage Fed-FDR fit is undefined.")
  }

  X_sub <- X[, external_support, drop = FALSE]
  validate_binary_outcome(y, paste0("site ", site_k, " second-stage outcome"))

  # Implementation-only patch retained from v1.1.1:
  # glmnet rejects a one-column input matrix, so append an excluded all-zero
  # dummy column, then remove it before refined debiasing.
  singleton_glmnet_patch <- (ncol(X_sub) == 1L)

  if (singleton_glmnet_patch) {
    X_glmnet <- cbind(X_sub, .glmnet_dummy = 0)
    glmnet_exclude <- 2L
  } else {
    X_glmnet <- X_sub
    glmnet_exclude <- NULL
  }

  set.seed(seed_for(cfg, rep_id, "second_stage", site_k))

  cv_args <- list(
    x = X_glmnet,
    y = y,
    family = "binomial",
    nfolds = cfg$nfolds,
    alpha = 1
  )
  if (!is.null(glmnet_exclude)) cv_args$exclude <- glmnet_exclude
  cvfit <- do.call(glmnet::cv.glmnet, cv_args)
  lambda <- cvfit$lambda.min

  fit_args <- list(
    x = X_glmnet,
    y = y,
    family = "binomial",
    lambda = lambda,
    alpha = 1
  )
  if (!is.null(glmnet_exclude)) fit_args$exclude <- glmnet_exclude
  fit <- do.call(glmnet::glmnet, fit_args)

  coef_all <- as.vector(stats::coef(fit))

  if (singleton_glmnet_patch) {
    lasso_est <- coef_all[c(1L, 2L)]
  } else {
    lasso_est <- coef_all
  }

  refined <- refined_debiased_lasso(
    X = X_sub,
    y = y,
    lasso_est = lasso_est,
    delta = cfg$hessian_delta
  )

  # RAW refined slopes for privacy attacks.
  beta_refined_full <- numeric(cfg$d)
  beta_refined_full[external_support] <- refined$est[-1L]

  se_refined_full <- rep(NA_real_, cfg$d)
  se_refined_full[external_support] <- refined$se[-1L]

  # Public-repository-compatible mirror input retained only for Fed-FDR
  # central sanity diagnostics; privacy attacks do NOT use this quantity.
  theta_mirror_full <- numeric(cfg$d)
  theta_active <- refined$est[-1L] / sqrt(refined$se[-1L])
  if (any(!is.finite(theta_active))) {
    stop("Official-style Theta mirror transformation produced non-finite values.")
  }
  theta_mirror_full[external_support] <- theta_active

  list(
    external_support = as.integer(external_support),
    lambda = as.numeric(lambda),
    refined_intercept = refined$est[1L],
    beta_refined_full = beta_refined_full,
    se_refined_full = se_refined_full,
    theta_mirror_full = theta_mirror_full,
    hessian_regularized = refined$hessian_regularized,
    singleton_glmnet_patch = singleton_glmnet_patch
  )
}

# ------------------------------------------------------------
# Fed-FDR central selection -- repository-compatible diagnostic
# ------------------------------------------------------------

fdr_control_pair <- function(beta1, beta2, q) {
  M <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
  tau_seq <- abs(M[M != 0])

  if (!any(M < 0)) {
    tau_hat <- 0
  } else {
    tau_seq <- sort(tau_seq)
    tau_hat <- 0
    for (tau in c(0, tau_seq)) {
      fdp_hat <- sum(M <= (-tau)) / max(sum(M >= tau), 1)
      tau_hat <- tau
      if (fdp_hat <= q) break
    }
  }

  S <- integer(length(M))
  S[which(M > tau_hat)] <- 1L
  list(selection = S, mirror = M, tau_hat = tau_hat)
}

mfdr_public_logic <- function(Theta, q) {
  K <- nrow(Theta)
  p <- ncol(Theta)
  pair_selections <- NULL
  pair_taus <- numeric(0)

  for (i in seq_len(K - 1L)) {
    for (j in (i + 1L):K) {
      fit <- fdr_control_pair(Theta[i, ], Theta[j, ], q)
      pair_selections <- rbind(pair_selections, fit$selection)
      pair_taus <- c(pair_taus, fit$tau_hat)
    }
  }

  a <- rowSums(pair_selections)
  denom <- a * as.integer(a >= 1) + as.integer(a < 1)
  inclusion <- colSums(pair_selections / denom) / (K * (K - 1) / 2)

  index <- order(inclusion)
  I_seq <- inclusion[index]
  index1 <- which(cumsum(I_seq) > q)
  supp <- index[index1]

  final <- integer(p)
  if (length(supp) > 0L) final[supp] <- 1L

  list(
    final_selection = final,
    inclusion = inclusion,
    pair_selection = pair_selections,
    pair_tau = pair_taus
  )
}

fit_fedfdr_one_rep <- function(site_data, rep_id, cfg) {
  K <- cfg$K

  first <- vector("list", K)
  for (k in seq_len(K)) {
    first[[k]] <- fit_first_stage_lasso(
      X = site_data[[k]]$X,
      y = site_data[[k]]$y,
      rep_id = rep_id,
      site_k = k,
      cfg = cfg
    )
  }

  support_local <- lapply(first, `[[`, "support")

  second <- vector("list", K)
  for (k in seq_len(K)) {
    others <- setdiff(seq_len(K), k)
    external_support <- sort(unique(unlist(support_local[others], use.names = FALSE)))
    second[[k]] <- fit_second_stage_refined(
      X = site_data[[k]]$X,
      y = site_data[[k]]$y,
      external_support = external_support,
      rep_id = rep_id,
      site_k = k,
      cfg = cfg
    )
  }

  beta_refined <- do.call(rbind, lapply(second, `[[`, "beta_refined_full"))
  theta_mirror <- do.call(rbind, lapply(second, `[[`, "theta_mirror_full"))

  central <- NULL
  if (isTRUE(cfg$compute_central)) {
    central <- mfdr_public_logic(theta_mirror, cfg$alpha)
  }

  list(
    first = first,
    second = second,
    support_local = support_local,
    beta_refined = beta_refined,
    theta_mirror = theta_mirror,
    central = central
  )
}

# ------------------------------------------------------------
# K-general source-site attack scores
# ------------------------------------------------------------

support_scores_K <- function(x, support_local, d, eps = 1e-12) {
  K <- length(support_local)
  stopifnot(K >= 2L)
  stopifnot(length(x) == d)

  S <- matrix(0, nrow = K, ncol = d)
  for (k in seq_len(K)) {
    if (length(support_local[[k]]) > 0L) {
      S[k, support_local[[k]]] <- 1
    }
  }

  total_selected_by_sites <- colSums(S)
  W <- matrix(0, nrow = K, ncol = d)

  for (k in seq_len(K)) {
    other_rate <- (total_selected_by_sites - S[k, ]) / (K - 1)
    W[k, ] <- S[k, ] * (1 - other_rate)
  }

  weight_sum <- rowSums(W)
  x2 <- x^2
  uS <- numeric(K)
  for (k in seq_len(K)) {
    if (weight_sum[k] > 0) {
      uS[k] <- sum(W[k, ] * x2) / weight_sum[k]
    } else {
      uS[k] <- 0
    }
  }

  QS <- uS / (sum(uS) + eps)

  list(
    Q = QS,
    u = uS,
    weights = W,
    weight_sum = weight_sum,
    support_matrix = S,
    zero_weight_sites = sum(weight_sum <= eps)
  )
}

coefficient_scores_K <- function(x, y, beta_refined, eps = 1e-12) {
  beta_refined <- as.matrix(beta_refined)
  K <- nrow(beta_refined)
  d <- ncol(beta_refined)
  stopifnot(length(x) == d)

  r <- as.vector(beta_refined %*% x)

  uX <- abs(r)
  QX <- uX / (sum(uX) + eps)

  ystar <- 2 * as.numeric(y) - 1
  margin <- ystar * r
  uXY <- safe_plogis(margin)
  QXY <- uXY / (sum(uXY) + eps)

  list(
    r = r,
    margin = margin,
    uX = uX,
    uXY = uXY,
    QX = QX,
    QXY = QXY
  )
}

expected_tie_decision <- function(scores, source_site, tie_tol = 1e-12) {
  scores <- as.numeric(scores)
  if (length(scores) < 2L || any(!is.finite(scores))) {
    return(list(
      pred = NA_integer_,
      credit = NA_real_,
      tie_size = NA_integer_,
      true_score = NA_real_,
      max_score = NA_real_,
      score_gap = NA_real_,
      tied_sites = integer(0)
    ))
  }

  max_score <- max(scores)
  tied <- which(abs(scores - max_score) <= tie_tol)
  tie_size <- length(tied)
  pred <- if (tie_size == 1L) as.integer(tied) else 0L
  credit <- if (source_site %in% tied) 1 / tie_size else 0

  sorted_scores <- sort(scores, decreasing = TRUE)
  score_gap <- if (length(sorted_scores) >= 2L) sorted_scores[1L] - sorted_scores[2L] else NA_real_

  list(
    pred = pred,
    credit = as.numeric(credit),
    tie_size = as.integer(tie_size),
    true_score = scores[source_site],
    max_score = max_score,
    score_gap = score_gap,
    tied_sites = as.integer(tied)
  )
}

attack_one_candidate_K <- function(x, y, source_site, fedfit, cfg) {
  K <- cfg$K

  scomp <- support_scores_K(
    x = x,
    support_local = fedfit$support_local,
    d = cfg$d,
    eps = cfg$score_eps
  )

  bcomp <- coefficient_scores_K(
    x = x,
    y = y,
    beta_refined = fedfit$beta_refined,
    eps = cfg$score_eps
  )

  Q0 <- rep(1 / K, K)
  Q_R1 <- scomp$Q
  Q_R2_X <- bcomp$QX
  Q_R2_XY <- bcomp$QXY
  Q_R3_X <- Q_R1 + Q_R2_X
  Q_R3_XY <- Q_R1 + Q_R2_XY

  score_list <- list(
    R0_X = Q0,
    R0_XY = Q0,
    R1_X = Q_R1,
    R1_XY = Q_R1,
    R2_X = Q_R2_X,
    R2_XY = Q_R2_XY,
    R3_X = Q_R3_X,
    R3_XY = Q_R3_XY
  )

  decisions <- lapply(score_list, function(z) {
    expected_tie_decision(z, source_site = source_site, tie_tol = cfg$tie_tol)
  })

  list(
    scores = score_list,
    decisions = decisions,
    support = scomp,
    coefficient = bcomp
  )
}

# ------------------------------------------------------------
# One repetition
# ------------------------------------------------------------

indices_to_string <- function(x) {
  if (length(x) == 0L) return("")
  paste(as.integer(x), collapse = ",")
}

calc_fdp_power <- function(selection01, true_support, d) {
  if (is.null(selection01)) {
    return(list(fdp = NA_real_, power = NA_real_, selected_size = NA_integer_))
  }
  selected <- which(selection01 != 0)
  true <- as.integer(true_support)
  false_count <- sum(!(selected %in% true))
  true_count <- sum(selected %in% true)
  list(
    fdp = false_count / max(length(selected), 1L),
    power = true_count / length(true),
    selected_size = length(selected)
  )
}

flatten_attack_K <- function(prefix, attack) {
  out <- list()
  for (nm in names(attack$decisions)) {
    dd <- attack$decisions[[nm]]
    out[[paste0(prefix, "_credit_", nm)]] <- dd$credit
    out[[paste0(prefix, "_pred_", nm)]] <- dd$pred
    out[[paste0(prefix, "_tie_size_", nm)]] <- dd$tie_size
    out[[paste0(prefix, "_true_score_", nm)]] <- dd$true_score
    out[[paste0(prefix, "_max_score_", nm)]] <- dd$max_score
    out[[paste0(prefix, "_score_gap_", nm)]] <- dd$score_gap
  }

  out
}

classify_failure_reason <- function(msg) {
  if (grepl("External support is empty", msg, fixed = TRUE)) return("empty_external_support")
  if (grepl("only one class", msg, fixed = TRUE)) return("single_class_outcome")
  if (grepl("Hessian", msg, ignore.case = TRUE)) return("refined_hessian_failure")
  if (grepl("non-finite", msg, ignore.case = TRUE) || grepl("invalid estimates", msg, ignore.case = TRUE)) {
    return("nonfinite_refined_coefficient")
  }
  if (grepl("glmnet", msg, ignore.case = TRUE)) return("glmnet_failure")
  "other"
}

run_one_rep_K <- function(rep_id, cfg, return_debug = FALSE) {
  stage <- "start"

  ans <- tryCatch({
    stage <- "parameter_generation"
    params <- generate_parameters(rep_id, cfg)

    stage <- "site_data_generation"
    site_data <- lapply(seq_len(cfg$K), function(k) {
      generate_site_data(rep_id, k, params, cfg)
    })

    stage <- "fedfdr_fit"
    fedfit <- fit_fedfdr_one_rep(site_data, rep_id, cfg)

    stage <- "member_candidate"
    set.seed(seed_for(cfg, rep_id, "member"))
    source_site <- sample(seq_len(cfg$K), size = 1L)
    member_index <- sample(seq_len(cfg$n_per_site), size = 1L)
    member <- list(
      x = as.numeric(site_data[[source_site]]$X[member_index, ]),
      y = as.integer(site_data[[source_site]]$y[member_index])
    )

    stage <- "fresh_candidate"
    fresh <- generate_fresh_observation(rep_id, source_site, params, cfg)

    stage <- "member_attack"
    member_attack <- attack_one_candidate_K(
      x = member$x,
      y = member$y,
      source_site = source_site,
      fedfit = fedfit,
      cfg = cfg
    )

    stage <- "fresh_attack"
    fresh_attack <- attack_one_candidate_K(
      x = fresh$x,
      y = fresh$y,
      source_site = source_site,
      fedfit = fedfit,
      cfg = cfg
    )

    stage <- "diagnostics"
    fed_diag <- calc_fdp_power(
      selection01 = if (!is.null(fedfit$central)) fedfit$central$final_selection else NULL,
      true_support = cfg$true_support,
      d = cfg$d
    )

    local_sizes <- vapply(fedfit$support_local, length, integer(1))
    external_sizes <- vapply(fedfit$second, function(z) length(z$external_support), integer(1))
    hess_reg <- vapply(fedfit$second, function(z) as.integer(z$hessian_regularized), integer(1))
    singleton_patch <- vapply(fedfit$second, function(z) as.integer(z$singleton_glmnet_patch), integer(1))

    pair_tau <- if (!is.null(fedfit$central)) fedfit$central$pair_tau else numeric(0)

    base <- list(
      rep_id = as.integer(rep_id),
      K = as.integer(cfg$K),
      seed = seed_for(cfg, rep_id, "rep"),
      status = "ok",
      failure_stage = "",
      failure_class = "",
      failure_reason = "",
      source_site = as.integer(source_site),
      member_index = as.integer(member_index),
      member_y = as.integer(member$y),
      fresh_y = as.integer(fresh$y),
      random_baseline = 1 / cfg$K,
      source_gamma = params$gamma[source_site],
      source_rho = params$rho[source_site],
      mean_gamma = mean(params$gamma),
      mean_rho = mean(params$rho),
      signal_scale = params$scale_factor,
      signal_signs = paste(params$signal_sign, collapse = ","),
      mean_local_support_size = mean(local_sizes),
      min_local_support_size = min(local_sizes),
      max_local_support_size = max(local_sizes),
      mean_external_support_size = mean(external_sizes),
      min_external_support_size = min(external_sizes),
      max_external_support_size = max(external_sizes),
      prop_external_full_p = mean(external_sizes >= cfg$d),
      n_external_singleton = sum(external_sizes == 1L),
      n_hessian_regularized = sum(hess_reg),
      prop_hessian_regularized = mean(hess_reg),
      n_singleton_glmnet_patch = sum(singleton_patch),
      final_selected_size = fed_diag$selected_size,
      fed_fdp = fed_diag$fdp,
      fed_power = fed_diag$power,
      mean_pair_tau = if (length(pair_tau) > 0L) mean(pair_tau) else NA_real_,
      min_pair_tau = if (length(pair_tau) > 0L) min(pair_tau) else NA_real_,
      max_pair_tau = if (length(pair_tau) > 0L) max(pair_tau) else NA_real_,
      member_R1_true_weight_sum = member_attack$support$weight_sum[source_site],
      member_R1_zero_weight_sites = member_attack$support$zero_weight_sites,
      fresh_R1_true_weight_sum = fresh_attack$support$weight_sum[source_site],
      fresh_R1_zero_weight_sites = fresh_attack$support$zero_weight_sites,
      member_R2_XY_true_margin = member_attack$coefficient$margin[source_site],
      member_R2_XY_max_margin = max(member_attack$coefficient$margin),
      fresh_R2_XY_true_margin = fresh_attack$coefficient$margin[source_site],
      fresh_R2_XY_max_margin = max(fresh_attack$coefficient$margin)
    )

    row <- c(
      base,
      flatten_attack_K("member", member_attack),
      flatten_attack_K("fresh", fresh_attack)
    )

    if (return_debug) {
      list(summary = row, debug = list(
        cfg = cfg,
        params = params,
        site_data = site_data,
        fedfit = fedfit,
        member = member,
        fresh = fresh,
        member_attack = member_attack,
        fresh_attack = fresh_attack
      ))
    } else {
      row
    }
  }, error = function(e) {
    msg <- conditionMessage(e)
    fail <- list(
      rep_id = as.integer(rep_id),
      K = as.integer(cfg$K),
      seed = seed_for(cfg, rep_id, "rep"),
      status = "failed",
      failure_stage = stage,
      failure_class = classify_failure_reason(msg),
      failure_reason = msg,
      random_baseline = 1 / cfg$K
    )
    if (return_debug) list(summary = fail, debug = NULL) else fail
  })

  ans
}

# ------------------------------------------------------------
# Data-frame / runtime helpers
# ------------------------------------------------------------

lists_to_data_frame <- function(rows) {
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  out <- lapply(rows, function(x) {
    miss <- setdiff(all_names, names(x))
    if (length(miss) > 0L) x[miss] <- NA
    x[all_names]
  })

  df <- as.data.frame(
    do.call(rbind, lapply(out, function(x) as.data.frame(x, stringsAsFactors = FALSE))),
    stringsAsFactors = FALSE
  )

  char_keep <- c(
    "status", "failure_stage", "failure_class", "failure_reason", "signal_signs"
  )

  for (nm in setdiff(names(df), char_keep)) {
    suppressWarnings(z <- as.numeric(df[[nm]]))
    original_nonmissing <- !is.na(df[[nm]]) & df[[nm]] != ""
    converted_ok <- !is.na(z) | !original_nonmissing
    if (all(converted_ok)) df[[nm]] <- z
  }

  df
}

format_seconds <- function(x) {
  if (!is.finite(x)) return("NA")
  x <- max(0, x)
  h <- floor(x / 3600)
  m <- floor((x %% 3600) / 60)
  s <- round(x %% 60)
  sprintf("%02dh:%02dm:%02ds", h, m, s)
}

# ------------------------------------------------------------
# Analysis helpers generalized to random baseline 1/K
# ------------------------------------------------------------

mean_ci_unbounded <- function(x, level = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0L) return(c(mean = NA, se = NA, lower = NA, upper = NA, n = 0))
  m <- mean(x)
  se <- if (n > 1L) stats::sd(x) / sqrt(n) else NA_real_
  z <- stats::qnorm(1 - (1 - level) / 2)
  c(
    mean = m,
    se = se,
    lower = if (is.finite(se)) m - z * se else NA_real_,
    upper = if (is.finite(se)) m + z * se else NA_real_,
    n = n
  )
}

mean_ci_probability <- function(x, level = 0.95) {
  z <- mean_ci_unbounded(x, level)
  z["lower"] <- max(0, z["lower"])
  z["upper"] <- min(1, z["upper"])
  z
}

paired_delta_ci <- function(x, y, level = 0.95) {
  ok <- is.finite(x) & is.finite(y)
  mean_ci_unbounded(x[ok] - y[ok], level = level)
}

one_sided_p_vs_baseline <- function(credit, baseline) {
  x <- credit[is.finite(credit)]
  if (length(x) < 2L) return(NA_real_)
  se <- stats::sd(x) / sqrt(length(x))
  m <- mean(x)
  if (!is.finite(se) || se == 0) {
    return(if (m > baseline) 0 else 1)
  }
  z <- (m - baseline) / se
  stats::pnorm(z, lower.tail = FALSE)
}

summarize_condition_K <- function(df, condition, K_value) {
  member_col <- paste0("member_credit_", condition)
  fresh_col <- paste0("fresh_credit_", condition)

  mc <- df[[member_col]]
  fc <- df[[fresh_col]]
  mci <- mean_ci_probability(mc)
  fci <- mean_ci_probability(fc)
  dci <- paired_delta_ci(mc, fc)
  baseline <- 1 / K_value

  data.frame(
    K = K_value,
    condition = condition,
    member_accuracy = unname(mci["mean"]),
    random_baseline = baseline,
    member_minus_random = unname(mci["mean"]) - baseline,
    member_mcse = unname(mci["se"]),
    member_ci_lower = unname(mci["lower"]),
    member_ci_upper = unname(mci["upper"]),
    one_sided_p_vs_random = one_sided_p_vs_baseline(mc, baseline),
    fresh_accuracy = unname(fci["mean"]),
    fresh_mcse = unname(fci["se"]),
    fresh_ci_lower = unname(fci["lower"]),
    fresh_ci_upper = unname(fci["upper"]),
    member_minus_fresh = unname(dci["mean"]),
    member_minus_fresh_mcse = unname(dci["se"]),
    member_minus_fresh_ci_lower = unname(dci["lower"]),
    member_minus_fresh_ci_upper = unname(dci["upper"]),
    n_valid = unname(mci["n"]),
    stringsAsFactors = FALSE
  )
}

summarize_K_results <- function(df, K_value) {
  ok <- df$status == "ok"
  valid <- df[ok, , drop = FALSE]
  conditions <- c("R0_X", "R0_XY", "R1_X", "R1_XY", "R2_X", "R2_XY", "R3_X", "R3_XY")

  attack <- do.call(rbind, lapply(conditions, function(cc) {
    summarize_condition_K(valid, cc, K_value)
  }))

  diagnostics <- data.frame(
    K = K_value,
    attempted = nrow(df),
    valid = sum(ok, na.rm = TRUE),
    failed = sum(!ok, na.rm = TRUE),
    failure_rate = mean(!ok, na.rm = TRUE),
    mean_local_support_size = mean(valid$mean_local_support_size, na.rm = TRUE),
    mean_external_support_size = mean(valid$mean_external_support_size, na.rm = TRUE),
    mean_max_external_support_size = mean(valid$max_external_support_size, na.rm = TRUE),
    prop_external_full_p = mean(valid$prop_external_full_p, na.rm = TRUE),
    mean_hessian_regularization_rate = mean(valid$prop_hessian_regularized, na.rm = TRUE),
    mean_singleton_patch_count = mean(valid$n_singleton_glmnet_patch, na.rm = TRUE),
    mean_fed_fdp = mean(valid$fed_fdp, na.rm = TRUE),
    mean_fed_power = mean(valid$fed_power, na.rm = TRUE),
    mean_final_selected_size = mean(valid$final_selected_size, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  list(attack = attack, diagnostics = diagnostics)
}
