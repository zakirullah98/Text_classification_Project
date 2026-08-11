Fed-FDR Source-Site Privacy Attack -- Frozen Protocol v1.1
==========================================================

WHAT TO RUN FIRST
-----------------
1) Put all files in one folder.
2) Open RStudio and set the working directory to this folder.
3) Install packages if needed:
     install.packages(c("glmnet", "mvtnorm"))
4) Run:
     source("02_RUN_QUICK_R20.R")
5) Send the resulting folder "results_quick_R20" back for inspection.

DO NOT start R=10,000 before the quick test has been checked.

FILES
-----
01_FedFDR_SourceAttack_functions.R
  All simulation, Fed-FDR, attack, and analysis helper functions.

02_RUN_QUICK_R20.R
  20-repetition debugging run. Run this first.

03_RUN_PILOT_R100.R
  100-repetition pilot with configurable parallel workers and ETA.

04_RUN_FORMAL_R10000.R
  Formal N=10,000 run. Checkpointed, resumable, deterministic repetition IDs.

05_ANALYZE_RESULTS.R
  Produces primary attack table, Full-regime decomposition, support mechanism
  diagnostics, fresh-control comparison, Fed-FDR sanity checks, and plots.

06_DEBUG_ONE_REP.R
  Prints one repetition in detail for debugging.

FROZEN STATISTICAL SETTINGS
---------------------------
K = 2
n1 = n2 = 50
d = 50
true support = {1,2,3,4,5}
N formal = 10,000
alpha = 0.10
gamma_k ~ U(-0.5, 0.5)
rho_k   ~ U(0.3, 0.5)
sigma_kj ~ U(7, 11)
common random sign per active feature across sites
beta_kj = sign_j * sigma_kj * sqrt(log(d)/n)
ordinary logistic Lasso with 5-fold CV in stage 1
refined debiased Lasso in stage 2
intercept is NOT released to the attacker

ATTACK CONDITIONS
-----------------
B0: exact no-information baseline, accuracy = 0.500
R0-X / R0-XY: final aggregated output; structural no-site-specific score
R1-X / R1-XY: support-only score dS
R2-X: coefficient score dBeta_X
R2-XY: coefficient score dBeta_XY
R3-X: dS + dBeta_X, fixed 1:1 weight
R3-XY: dS + dBeta_XY, fixed 1:1 weight

FULL REGIME REPORTING
---------------------
The analysis always reports:
- support component accuracy
- coefficient component accuracy
- 1:1 combined accuracy
- agreement / conflict / component-tie rates

FRESH CONTROL
-------------
For every repetition, one fresh observation is independently generated from
exactly the same true source-site DGP (same gamma, rho, Sigma, beta), and is
never used in training, CV, screening, Hessian estimation, or model fitting.

IMPORTANT IMPLEMENTATION NOTES
------------------------------
- Privacy attacks use RAW refined slope coefficients.
- Fed-FDR feature selection uses the repository-compatible Theta mirror input.
- The intercept is retained internally only for data/model construction and is
  never passed to attack scoring functions.
- In K=2, coefficient coordinates can reveal the other site's local support;
  therefore R2 vs R3 is mainly a comparison of information-use rules, not a
  clean comparison of additional information content.
- Formal repetitions are pre-numbered. Failed draws are logged; they are not
  silently replaced by easier draws.
