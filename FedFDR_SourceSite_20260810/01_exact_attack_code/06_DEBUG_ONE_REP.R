# ============================================================
# Deep debug of ONE repetition.
# Use this if QUICK R20 reports a failure or a suspicious result.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

source("01_FedFDR_SourceAttack_functions.R")
check_required_packages()

REP_ID <- 1L
cfg <- make_config(n_reps = 1L)

obj <- run_one_rep(REP_ID, cfg, return_debug = TRUE)
print(obj$summary)

if (!is.null(obj$debug)) {
  cat("\n--- TRUE PARAMETERS ---\n")
  print(obj$debug$params[c("gamma", "rho", "signal_sign", "sigma_active")])

  cat("\n--- LOCAL SUPPORTS ---\n")
  print(obj$debug$fedfit$support_local)

  cat("\n--- EXTERNAL SUPPORTS ---\n")
  print(lapply(obj$debug$fedfit$second, `[[`, "external_support"))

  cat("\n--- RAW REFINED BETA, nonzero coordinates ---\n")
  for (k in 1:2) {
    b <- obj$debug$fedfit$beta_refined[k, ]
    print(data.frame(site = k, j = which(b != 0), beta = b[b != 0]))
  }

  cat("\n--- THETA MIRROR, nonzero coordinates ---\n")
  for (k in 1:2) {
    th <- obj$debug$fedfit$theta_mirror[k, ]
    print(data.frame(site = k, j = which(th != 0), theta = th[th != 0]))
  }

  cat("\n--- MEMBER ATTACK ---\n")
  print(obj$debug$member_attack)

  cat("\n--- FRESH ATTACK ---\n")
  print(obj$debug$fresh_attack)
}
