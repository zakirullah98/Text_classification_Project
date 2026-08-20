============================================================
Task A: Source-Site Attack Sensitivity to the Number of Sites
K-General Source Attack v1.0
============================================================

THIS PACKAGE IS THE NEW-WEEK EXTENSION OF THE FROZEN K=2 SOURCE-SITE
ATTACK EXPERIMENT. IT DOES NOT OVERWRITE THE OLD N=10,000 CODE.

Main question
-------------
For K in {3,5,10,20,50,100}, can a known training member's source site
still be inferred better than the 1/K random-guessing baseline?

Fixed settings for Task A
-------------------------
- n_k = 50 for every site
- p = 50
- p1 = 5, true support = {1,...,5}
- alpha = 0.1
- Same Gaussian-Toeplitz DGP as last week's K=2 experiment
- Same first-stage and second-stage logistic Lasso fitting
- Same refined debiasing implementation
- Same singleton glmnet patch
- Same matched-fresh control idea
- Only K and the source-site decision rule are generalized

Frozen K-general attack rule
----------------------------
For each candidate, every site k gets a score Q_k.
Predict the source by argmax_k Q_k.
If multiple sites tie for the maximum, use EXPECTED uniform tie credit:
  credit = 1 / number_of_tied_sites if the true site is tied, else 0.

R0:
  Q_0k = 1/K for every site (structural baseline).

R1 support score:
  s_kj = 1{j in S_hat^(k)}
  w_kj = s_kj * [1 - (1/(K-1))*sum_{l != k} s_lj]
  u_Sk = weighted mean of x_j^2 using w_kj
  Q_Sk = u_Sk / (sum_l u_Sl + eps)

R2 X-only:
  r_k = x' beta_hat^(k)
  u_k = |r_k|
  Q_betaX,k = u_k / (sum_l u_l + eps)

R2 (X,Y):
  m_k = (2y-1) x' beta_hat^(k)
  u_k = expit(m_k)
  Q_betaXY,k = u_k / (sum_l u_l + eps)
  This has the same argmax as choosing the largest signed margin /
  smallest slope-only logistic loss.

R3:
  Q_Full,k = Q_Sk + Q_betak
  Fixed equal weight; no tuning.

IMPORTANT K=2 property
----------------------
At K=2, the new definitions algebraically reduce to the old D-scores:
  Q_1 - Q_2 = D.
Therefore the new implementation MUST reproduce the old K=2 decisions.

RUN ORDER
---------
1) 02_TEST_K2_BACKWARD_COMPATIBILITY.R
   - compares the new code directly with the exact old Frozen v1.1.1 code
   - MUST print PASS
   - if it fails, STOP and do not run K>2

2) 03_RUN_QUICK_KGRID.R
   - tiny N only
   - checks K = 2,3,5,10,20,50,100 can run
   - inspect runtime, external-support growth, Hessian regularization,
     and failure classes
   - DO NOT interpret quick accuracies scientifically

3) 04_RUN_PILOT_KGRID.R
   - default N=20 per K
   - designed for runtime / numerical-stability / failure diagnostics
   - resumable by K and batch

4) 06_SUMMARIZE_KGRID.R
   - default INPUT_STAGE="PILOT"
   - produces attack summary, diagnostics, and figures

5) 05_RUN_FORMAL_KGRID_TEMPLATE.R
   - intentionally has N_FORMAL = NA
   - DO NOT run until the pilot tells us what runtime is feasible and
     we choose a justified formal Monte Carlo size

Output metrics
--------------
For every K and every regime/candidate condition, report:
- member source-site accuracy
- random baseline = 1/K
- member accuracy - 1/K
- matched-fresh accuracy
- member - fresh difference
- MCSE / confidence intervals

Diagnostics by K
----------------
- failure rate and failure class
- mean/min/max local support size
- mean/min/max external support size
- proportion external support reaches all p variables
- Hessian regularization rate
- singleton glmnet patch count
- Fed-FDR sanity metrics (repository-compatible central step; diagnostic only)

The two main Task-A questions are:
1) Does Acc_member(K) - 1/K remain positive as K increases?
2) Does Acc_member(K) - Acc_fresh(K) remain positive as K increases?

============================================================
