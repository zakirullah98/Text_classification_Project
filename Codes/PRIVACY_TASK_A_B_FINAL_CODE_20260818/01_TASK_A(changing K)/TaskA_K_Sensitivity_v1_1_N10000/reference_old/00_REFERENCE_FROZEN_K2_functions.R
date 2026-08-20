# ============================================================
# Fed-FDR source-site privacy attack experiment
# Frozen Protocol v1.1 -- common functions
# ============================================================
# This file implements the protocol frozen in the project discussion:
#   K = 2, n1 = n2 = 50, d = 50, N = 10,000 formal repetitions.
#   Target: source-site inference for a known training observation.
#   Candidate information: X-only or (X,Y).
#   Information regimes: R0 final aggregated, R1 support-only,
#                        R2 site-specific coefficients, R3 full 2K objects.
#   Intercepts are NOT released to the attacker.
#   R3 combines support and coefficient scores with fixed 1:1 weight.
#   Fresh controls are independently generated from the same source-site DGP.
#
# IMPORTANT:
# - The Fed-FDR selection path preserves the public repository's mirror-input
#   transformation Theta = refined_est / sqrt(refined_se) on active coordinates.
# - Privacy attacks use the RAW refined slope coefficients, not Theta.
# - The first-stage screening uses ordinary GLM-Lasso (weight = FALSE),
#   matching Algorithm 2 in the paper rather than the repository's weighted
#   default option.
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
    n_reps = 20L,
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
    tie_tol = 1e-12
) {
  stopifnot(K == 2L)
  stopifnot(length(true_support) == s)
  stopifnot(all(true_support >= 1L & true_support <= d))

  list(
    protocol_version = "Frozen Protocol v1.1",
    implementation_version = "v1.1.1-singleton-glmnet-patch",
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
    tie_tol = tie_tol
  )
}

seed_stream <- function(cfg, rep_id, stream_id) {
  # Deterministic sub-streams make the simulation reproducible even if the
  # parallel execution order changes.
  z <- as.double(cfg$base_seed) + as.double(rep_id) * 1000 + as.double(stream_id) * 37
  z <- z %% (.Machine$integer.max - 1)
  as.integer(z + 1)
}

safe_plogis <- function(x) stats::plogis(x)

# ------------------------------------------------------------
# Data generation
# ------------------------------------------------------------

generate_parameters <- function(rep_id, cfg) {
  set.seed(seed_stream(cfg, rep_id, 1L))

  gamma <- stats::runif(cfg$K, cfg$gamma_range[1], cfg$gamma_range[2])
  rho <- stats::runif(cfg$K, cfg$rho_range[1], cfg$rho_range[2])

  # Common sign across sites; site-specific signal magnitudes.
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
  set.seed(seed_stream(cfg, rep_id, 10L + as.integer(site_k)))

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
  # Independent draw from exactly the SAME source-site DGP used in the current
  # repetition. It is never included in training, CV, screening, Hessian, etc.
  set.seed(seed_stream(cfg, rep_id, 51L))
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
# Fed-FDR fitting helpers
# ------------------------------------------------------------

validate_binary_outcome <- function(y, label = "outcome") {
  if (length(unique(y)) < 2L) {
    stop(label, " contains only one class.")
  }
  invisible(TRUE)
}

fit_first_stage_lasso <- function(X, y, rep_id, site_k, cfg) {
  validate_binary_outcome(y, paste0("site ", site_k, " first-stage outcome"))
  set.seed(seed_stream(cfg, rep_id, 20L + as.integer(site_k)))

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
    # This remains a genuine finite-sample Fed-FDR boundary event.
    # Do NOT repair it by force-selecting a variable.
    stop("External support is empty; second-stage Fed-FDR fit is undefined.")
  }

  X_sub <- X[, external_support, drop = FALSE]
  validate_binary_outcome(y, paste0("site ", site_k, " second-stage outcome"))

  # ----------------------------------------------------------
  # IMPLEMENTATION-ONLY PATCH FOR |S^(-k)| = 1
  # glmnet requires an input matrix with at least two columns even though
  # the one-predictor penalized logistic problem is mathematically valid.
  # We therefore append an all-zero dummy column and explicitly exclude it
  # from penalization/fitting. The refined debiasing step is still performed
  # on the original one-column X_sub, so the statistical target is unchanged.
  # ----------------------------------------------------------
  singleton_glmnet_patch <- (ncol(X_sub) == 1L)

  if (singleton_glmnet_patch) {
    X_glmnet <- cbind(X_sub, .glmnet_dummy = 0)
    glmnet_exclude <- 2L
  } else {
    X_glmnet <- X_sub
    glmnet_exclude <- NULL
  }

  set.seed(seed_stream(cfg, rep_id, 30L + as.integer(site_k)))

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
    # coef_all = (intercept, real predictor, excluded dummy predictor).
    # refined_debiased_lasso() must receive only intercept + REAL predictor.
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

  # Privacy attack object: RAW refined slope coefficients, padded to d.
  beta_refined_full <- numeric(cfg$d)
  beta_refined_full[external_support] <- refined$est[-1L]

  se_refined_full <- rep(NA_real_, cfg$d)
  se_refined_full[external_support] <- refined$se[-1L]

  # Fed-FDR mirror input: intentionally reproduce the public repository's
  # expression fit0$est / sqrt(fit0$se), excluding the intercept.
  theta_mirror_full <- numeric(cfg$d)
  theta_active <- refined$est[-1L] / sqrt(refined$se[-1L])
  if (any(!is.finite(theta_active))) {
    stop("Official-style Theta mirror transformation produced non-finite values.")
  }
  theta_mirror_full[external_support] <- theta_active

  list(
    external_support = as.integer(external_support),
    lambda = as.numeric(lambda),
    refined_intercept = refined$est[1L],      # kept internally; NOT released to attack functions
    beta_refined_full = beta_refined_full,
    se_refined_full = se_refined_full,
    theta_mirror_full = theta_mirror_full,
    hessian_regularized = refined$hessian_regularized,
    singleton_glmnet_patch = singleton_glmnet_patch
  )
}

# ------------------------------------------------------------
# Fed-FDR central selection -- public-repository-compatible logic
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

  central <- mfdr_public_logic(theta_mirror, cfg$alpha)

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
# Attack scores
# ------------------------------------------------------------

support_component_details <- function(x, support1, support2, true_support, eps) {
  D1 <- setdiff(support1, support2)
  D2 <- setdiff(support2, support1)

  u1 <- if (length(D1) > 0L) mean(x[D1]^2) else 0
  u2 <- if (length(D2) > 0L) mean(x[D2]^2) else 0
  dS <- (u1 - u2) / (u1 + u2 + eps)

  list(
    dS = dS,
    uS1 = u1,
    uS2 = u2,
    D1 = as.integer(D1),
    D2 = as.integer(D2),
    D1_signal = sum(D1 %in% true_support),
    D1_null = sum(!(D1 %in% true_support)),
    D2_signal = sum(D2 %in% true_support),
    D2_null = sum(!(D2 %in% true_support))
  )
}

coefficient_scores <- function(x, y, beta1, beta2, eps) {
  r1 <- sum(x * beta1)
  r2 <- sum(x * beta2)

  uX1 <- abs(r1)
  uX2 <- abs(r2)
  d_beta_x <- (uX1 - uX2) / (uX1 + uX2 + eps)

  ystar <- 2 * as.numeric(y) - 1
  m1 <- ystar * r1
  m2 <- ystar * r2
  uXY1 <- safe_plogis(m1)
  uXY2 <- safe_plogis(m2)
  d_beta_xy <- (uXY1 - uXY2) / (uXY1 + uXY2 + eps)

  list(
    r1 = r1,
    r2 = r2,
    m1 = m1,
    m2 = m2,
    uBetaX1 = uX1,
    uBetaX2 = uX2,
    uBetaXY1 = uXY1,
    uBetaXY2 = uXY2,
    d_beta_x = d_beta_x,
    d_beta_xy = d_beta_xy
  )
}

score_to_credit <- function(score, source_site, tie_tol = 1e-12) {
  if (is.na(score) || !is.finite(score)) return(NA_real_)
  if (abs(score) <= tie_tol) return(0.5)
  pred <- if (score > 0) 1L else 2L
  as.numeric(pred == source_site)
}

score_to_pred <- function(score, tie_tol = 1e-12) {
  if (is.na(score) || !is.finite(score)) return(NA_integer_)
  if (abs(score) <= tie_tol) return(0L) # 0 denotes tie
  if (score > 0) 1L else 2L
}

component_relation <- function(dS, dB, tie_tol = 1e-12) {
  if (!is.finite(dS) || !is.finite(dB)) return(NA_character_)
  if (abs(dS) <= tie_tol || abs(dB) <= tie_tol) return("component_tie")
  if (sign(dS) == sign(dB)) "agreement" else "conflict"
}

attack_one_candidate <- function(x, y, source_site, fedfit, cfg) {
  support1 <- fedfit$support_local[[1L]]
  support2 <- fedfit$support_local[[2L]]
  beta1 <- fedfit$beta_refined[1L, ]
  beta2 <- fedfit$beta_refined[2L, ]

  scomp <- support_component_details(
    x = x,
    support1 = support1,
    support2 = support2,
    true_support = cfg$true_support,
    eps = cfg$score_eps
  )
  bcomp <- coefficient_scores(
    x = x,
    y = y,
    beta1 = beta1,
    beta2 = beta2,
    eps = cfg$score_eps
  )

  d_full_x <- scomp$dS + bcomp$d_beta_x
  d_full_xy <- scomp$dS + bcomp$d_beta_xy

  # R0: structural no-site-specific score. Credit is fixed at 0.5.
  # This is distinct conceptually from B0, but has the same numerical credit
  # under the frozen attack rule.
  scores <- list(
    R0_X = NA_real_,
    R0_XY = NA_real_,
    R1_X = scomp$dS,
    R1_XY = scomp$dS,
    R2_X = bcomp$d_beta_x,
    R2_XY = bcomp$d_beta_xy,
    R3_X = d_full_x,
    R3_XY = d_full_xy
  )

  credits <- list(
    R0_X = 0.5,
    R0_XY = 0.5,
    R1_X = score_to_credit(scores$R1_X, source_site, cfg$tie_tol),
    R1_XY = score_to_credit(scores$R1_XY, source_site, cfg$tie_tol),
    R2_X = score_to_credit(scores$R2_X, source_site, cfg$tie_tol),
    R2_XY = score_to_credit(scores$R2_XY, source_site, cfg$tie_tol),
    R3_X = score_to_credit(scores$R3_X, source_site, cfg$tie_tol),
    R3_XY = score_to_credit(scores$R3_XY, source_site, cfg$tie_tol)
  )

  preds <- list(
    R0_X = 0L,
    R0_XY = 0L,
    R1_X = score_to_pred(scores$R1_X, cfg$tie_tol),
    R1_XY = score_to_pred(scores$R1_XY, cfg$tie_tol),
    R2_X = score_to_pred(scores$R2_X, cfg$tie_tol),
    R2_XY = score_to_pred(scores$R2_XY, cfg$tie_tol),
    R3_X = score_to_pred(scores$R3_X, cfg$tie_tol),
    R3_XY = score_to_pred(scores$R3_XY, cfg$tie_tol)
  )

  list(
    scores = scores,
    credits = credits,
    preds = preds,
    support = scomp,
    coefficient = bcomp,
    relation_x = component_relation(scomp$dS, bcomp$d_beta_x, cfg$tie_tol),
    relation_xy = component_relation(scomp$dS, bcomp$d_beta_xy, cfg$tie_tol)
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

flatten_attack <- function(prefix, attack) {
  out <- list()
  for (nm in names(attack$scores)) {
    out[[paste0(prefix, "_score_", nm)]] <- attack$scores[[nm]]
    out[[paste0(prefix, "_credit_", nm)]] <- attack$credits[[nm]]
    out[[paste0(prefix, "_pred_", nm)]] <- attack$preds[[nm]]
  }
  out[[paste0(prefix, "_dS")]] <- attack$support$dS
  out[[paste0(prefix, "_dBeta_X")]] <- attack$coefficient$d_beta_x
  out[[paste0(prefix, "_dBeta_XY")]] <- attack$coefficient$d_beta_xy
  out[[paste0(prefix, "_relation_X")]] <- attack$relation_x
  out[[paste0(prefix, "_relation_XY")]] <- attack$relation_xy
  out
}

run_one_rep <- function(rep_id, cfg, return_debug = FALSE) {
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
    set.seed(seed_stream(cfg, rep_id, 41L))
    source_site <- sample(seq_len(cfg$K), size = 1L)
    member_index <- sample(seq_len(cfg$n_per_site), size = 1L)
    member <- list(
      x = as.numeric(site_data[[source_site]]$X[member_index, ]),
      y = as.integer(site_data[[source_site]]$y[member_index])
    )

    stage <- "fresh_candidate"
    fresh <- generate_fresh_observation(rep_id, source_site, params, cfg)

    stage <- "member_attack"
    member_attack <- attack_one_candidate(
      x = member$x,
      y = member$y,
      source_site = source_site,
      fedfit = fedfit,
      cfg = cfg
    )

    stage <- "fresh_attack"
    fresh_attack <- attack_one_candidate(
      x = fresh$x,
      y = fresh$y,
      source_site = source_site,
      fedfit = fedfit,
      cfg = cfg
    )

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

    base <- list(
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
      rho1 = params$rho[1L],
      rho2 = params$rho[2L],
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
    fail <- list(
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
# Data-frame / batch helpers
# ------------------------------------------------------------

lists_to_data_frame <- function(rows) {
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  out <- lapply(rows, function(x) {
    miss <- setdiff(all_names, names(x))
    if (length(miss) > 0L) x[miss] <- NA
    x[all_names]
  })
  df <- as.data.frame(do.call(rbind, lapply(out, function(x) as.data.frame(x, stringsAsFactors = FALSE))),
                      stringsAsFactors = FALSE)

  # Convert obvious numeric/logical fields back from character if rbind coerced them.
  char_keep <- c(
    "status", "failure_stage", "failure_reason",
    "support1", "support2", "external1", "external2", "final_support",
    "member_relation_X", "member_relation_XY", "fresh_relation_X", "fresh_relation_XY"
  )
  for (nm in setdiff(names(df), char_keep)) {
    suppressWarnings({
      z <- as.numeric(df[[nm]])
    })
    # Only replace if every non-missing original value converted successfully.
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
# Analysis helpers
# ------------------------------------------------------------

mean_ci <- function(x, level = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0L) return(c(mean = NA, se = NA, lower = NA, upper = NA, n = 0))
  m <- mean(x)
  se <- if (n > 1L) stats::sd(x) / sqrt(n) else NA_real_
  z <- stats::qnorm(1 - (1 - level) / 2)
  lo <- if (is.finite(se)) m - z * se else NA_real_
  hi <- if (is.finite(se)) m + z * se else NA_real_
  c(mean = m, se = se, lower = max(0, lo), upper = min(1, hi), n = n)
}

paired_delta_ci <- function(x, y, level = 0.95) {
  ok <- is.finite(x) & is.finite(y)
  d <- x[ok] - y[ok]
  n <- length(d)
  if (n == 0L) return(c(delta = NA, se = NA, lower = NA, upper = NA, n = 0))
  m <- mean(d)
  se <- if (n > 1L) stats::sd(d) / sqrt(n) else NA_real_
  z <- stats::qnorm(1 - (1 - level) / 2)
  c(
    delta = m,
    se = se,
    lower = if (is.finite(se)) m - z * se else NA_real_,
    upper = if (is.finite(se)) m + z * se else NA_real_,
    n = n
  )
}

one_sided_p_vs_half <- function(credit) {
  x <- credit[is.finite(credit)]
  if (length(x) < 2L) return(NA_real_)
  se <- stats::sd(x) / sqrt(length(x))
  m <- mean(x)
  if (!is.finite(se) || se == 0) {
    return(if (m > 0.5) 0 else 1)
  }
  z <- (m - 0.5) / se
  stats::pnorm(z, lower.tail = FALSE)
}

auc_rank <- function(score, source_site) {
  ok <- is.finite(score) & source_site %in% c(1, 2)
  s <- score[ok]
  lab <- source_site[ok]
  n1 <- sum(lab == 1)
  n0 <- sum(lab == 2)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(s, ties.method = "average")
  (sum(r[lab == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

summarize_condition <- function(df, condition) {
  member_credit_col <- paste0("member_credit_", condition)
  fresh_credit_col <- paste0("fresh_credit_", condition)
  member_score_col <- paste0("member_score_", condition)
  fresh_score_col <- paste0("fresh_score_", condition)

  mc <- df[[member_credit_col]]
  fc <- df[[fresh_credit_col]]
  mci <- mean_ci(mc)
  fci <- mean_ci(fc)
  dci <- paired_delta_ci(mc, fc)

  member_auc <- if (member_score_col %in% names(df)) auc_rank(df[[member_score_col]], df$source_site) else NA_real_
  fresh_auc <- if (fresh_score_col %in% names(df)) auc_rank(df[[fresh_score_col]], df$source_site) else NA_real_

  data.frame(
    condition = condition,
    member_accuracy = unname(mci["mean"]),
    attack_advantage = unname(mci["mean"]) - 0.5,
    member_ci_lower = unname(mci["lower"]),
    member_ci_upper = unname(mci["upper"]),
    one_sided_p_vs_0_5 = one_sided_p_vs_half(mc),
    fresh_accuracy = unname(fci["mean"]),
    fresh_ci_lower = unname(fci["lower"]),
    fresh_ci_upper = unname(fci["upper"]),
    delta_train = unname(dci["delta"]),
    delta_train_ci_lower = unname(dci["lower"]),
    delta_train_ci_upper = unname(dci["upper"]),
    member_auc = member_auc,
    fresh_auc = fresh_auc,
    n_valid = unname(mci["n"]),
    stringsAsFactors = FALSE
  )
}

summarize_relation <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n == 0L) return(c(agreement = NA, conflict = NA, component_tie = NA, n = 0))
  c(
    agreement = mean(x == "agreement"),
    conflict = mean(x == "conflict"),
    component_tie = mean(x == "component_tie"),
    n = n
  )
}
