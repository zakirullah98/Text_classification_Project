# Task B Frozen Population Protocol v1.0

## Question
With K=2, does the previously observed source-site attack remain detectable when the covariate population changes?

## Invariant components
- K=2, n1=n2=50, p=50.
- Formal N=10,000 attempted repetitions **per main population**.
- True support S1={1,...,5}.
- Common active signs across sites; site-specific active magnitudes as in Frozen source-site Protocol v1.1.1.
- Logistic response model, Fed-FDR fitting, R0--R3 attacks, candidate sampling, and matched-fresh control are unchanged.
- Intercepts and population parameters are not released to the attacker.
- Matched fresh is generated from the same realized source-site population parameters as the member, but is not used in training.

## Main populations

### P0: Toeplitz Gaussian baseline
X_k ~ N(0,Sigma_k), (Sigma_k)_{rs}=rho_k^{|r-s|}, rho_k ~ U(0.3,0.5), independently by site.

### P1: Gaussian-copula correlated Uniform
Generate latent Z_k ~ N(0,Sigma_k) with the same Toeplitz latent structure and rho_k ~ U(0.3,0.5). Set
U_j=Phi(Z_j), X_j=2*sqrt(3)*(U_j-1/2).
Thus each marginal is U(-sqrt(3),sqrt(3)) with mean 0 and variance 1. The realized Pearson correlation of X is not asserted to equal the latent Gaussian correlation.

### P2: Fixed-permuted block-sparse Gaussian
Start from ten 5x5 Toeplitz blocks with within-block parameter rho_k ~ U(0.75,0.90), zero covariance between blocks, and use one fixed variable permutation for both sites and all repetitions. The five active variables 1,...,5 are deliberately placed in five different covariance blocks, preventing artificial alignment of the true support with one strong-correlation block.

### P3: Tridiagonal near-singular Gaussian
Sigma has unit diagonal, first off-diagonal rho_k, and all entries with |r-s|>1 equal zero. rho_k ~ U(0.45,0.49). Quick/Pilot explicitly records lambda_min(Sigma), condition number, and det(Sigma)^(1/p). Any non-positive minimum eigenvalue is treated as invalid rather than silently repaired.

### P4: Theoretically standardized spatial-Poisson grid counts
Operationalization: each observation corresponds to an independent homogeneous spatial Poisson point pattern on [0,1]^2, discretized into 5x10=50 equal cells. Let lambda_k ~ U(1,4) denote the expected count per cell, corresponding to spatial intensity 50*lambda_k per unit area. For cell count C_j ~ Poisson(lambda_k), define
X_j=(C_j-lambda_k)/sqrt(lambda_k).
This is theoretical standardization using the true realized lambda_k in the data generator; lambda_k is not released to the attacker. Standardization leaves site-dependent higher-moment/shape differences when lambda_1 != lambda_2.

## Optional sensitivities
- P1-IID: iid U(-sqrt(3),sqrt(3)); removes dependence structure.
- P4-RAW: raw Poisson cell counts C_j; preserves mean/intensity scale differences.
These are secondary and should be run only if needed to interpret P1/P4 main results.

## DGP sanity diagnostics (Quick/Pilot only)
Use an independent oracle diagnostic sample that never enters fitting or attacks. Record by site:
- SD(X^T beta*),
- SD(eta),
- proportion with pi<0.05 or pi>0.95,
- mean pi,
- diagnostic y prevalence,
- oracle AUC using the true eta as score.
For covariance scenarios also record relevant eigenvalue/conditioning metrics, especially P3.

These diagnostics are used to detect genuine DGP degeneration, not to post-hoc tune parameters until populations have identical difficulty. If Quick/Pilot is numerically and statistically sane, parameters are frozen for Formal N=10,000.

## Primary estimands
For R2-(X,Y) and R3-(X,Y), report:
- member source-site accuracy,
- member accuracy - 0.5,
- matched-fresh accuracy,
- member accuracy - matched-fresh accuracy,
- Monte Carlo uncertainty / 95% intervals.
All R0--R3 conditions remain available in the full table.

## Interpretation rule
- Member > 0.5 and Member > Fresh: attack signal persists and is stronger for actual training members under this population.
- Member > 0.5 but Member approximately Fresh: source identification may be driven substantially by population-level site differences.
- Member approximately 0.5 and Fresh approximately 0.5: the pre-specified attack does not detect clear source-site signal in this population; this is not a proof of no leakage.
- Do not post-hoc reverse below-chance attack directions.
