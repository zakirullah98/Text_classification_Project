import numpy as np
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')
# Disable numpy scientific notation, display up to 4 decimal places
np.set_printoptions(
    suppress=True,  # Core parameter: disable scientific notation
    precision=4,    # Keep 4 decimal places
    floatmode='fixed'
)
# ==========================================
# 0. Dataset Generation (Section 5 Simulation)
# ==========================================
np.random.seed(42)

K = 3          # Number of sites
p = 5          # Feature dimension
n = 200        # Sample size per site
alpha = 0.2    # Target FDR control level

# True coefficients: features 0 and 1 are signals, 2-4 are noise
true_beta = np.array([2.0, -1.5, 0.0, 0.0, 0.0])
true_support = {0, 1}

# Generate distributed data
data_sites = []
for k in range(K):
    X_k = np.random.randn(n, p)
    logits = X_k @ true_beta
    probs = 1 / (1 + np.exp(-logits))
    y_k = np.random.binomial(1, probs)
    data_sites.append((X_k, y_k))

# Inspect data
print("\n--- Inspecting Site 1 Data ---")
X_0, y_0 = data_sites[0]
print(f"Site 1: {X_0.shape[0]} samples, {X_0.shape[1]} features")
print("\nFirst 5 rows of X:")
print(np.round(X_0[:5], 2))
print("\nFirst 5 labels y:")
print(y_0[:5])
print("-" * 50, "\n")

# ==========================================
# 1. Algorithm 2: Local Site Task
# ==========================================
print("--- Algorithm 2: Local Support Extraction and Swap ---")
local_supports = []

# Step 1: Local Lasso screening
for k in range(K):
    X_k, y_k = data_sites[k]
    lasso = LogisticRegressionCV(
        cv=5, penalty='l1', solver='liblinear', max_iter=1000, random_state=42
    )
    lasso.fit(X_k, y_k)
    support_k = set(np.where(np.abs(lasso.coef_[0]) > 1e-5)[0])
    local_supports.append(support_k)
    print(f"Site {k+1} Lasso support: {support_k}")

# ==============================================
# 教学演示开关：注释掉下面这行使用真实Lasso结果
# ==============================================
print("\n⚠️  [Teaching Demo] Manually setting divergent supports to demonstrate noise conflict")
local_supports = [{1, 2}, {0, 3}, {0, 1, 4}]
print(f"\nLocal supports: {local_supports}")
# ==============================================

# Step 2: Between-site swap and refined desparsified Lasso
beta_hat_lower = np.zeros((K, p))
print("\n--- Step 2: Support Swap and Debiased Estimation ---")

for k in range(K):
    X_k, y_k = data_sites[k]
    
    # Construct S^(-k): union of supports from all other sites
    S_minus_k = set()
    for j in range(K):
        if j != k:
            S_minus_k.update(local_supports[j])
    S_minus_k = sorted(list(S_minus_k))
    print(f"\nSite {k+1} S^(-k) (borrowed from others): {S_minus_k}")
    
    if len(S_minus_k) > 0:
        # Section 3.2 of the paper: Refined desparsified Lasso
        lasso_sub = LogisticRegressionCV(
            cv=5, penalty='l1', solver='liblinear', max_iter=1000, random_state=42
        )
        lasso_sub.fit(X_k[:, S_minus_k], y_k)
        beta_lasso_sub = lasso_sub.coef_[0]
        
        # Compute gradient and Hessian
        z = X_k[:, S_minus_k] @ beta_lasso_sub
        prob = 1 / (1 + np.exp(-z))
        grad = (1/n) * X_k[:, S_minus_k].T @ (prob - y_k)
        W = prob * (1 - prob)
        Hessian = (1/n) * X_k[:, S_minus_k].T @ (W[:, np.newaxis] * X_k[:, S_minus_k])
        
        # Invert Hessian (low-dimensional guarantee)
        try:
            Theta = np.linalg.inv(Hessian)
        except np.linalg.LinAlgError:
            Theta = np.linalg.inv(Hessian + 1e-4 * np.eye(len(S_minus_k)))
        
        # Debias
        beta_desparsified_sub = beta_lasso_sub - Theta @ grad
        beta_hat_lower[k, S_minus_k] = beta_desparsified_sub
        
        print(f"    ↳ Debiased coefficients: {np.round(beta_desparsified_sub, 3)}")
        print(f"    ↳ Sparse payload: indices {S_minus_k}, values {np.round(beta_desparsified_sub, 3)}")

print("\n" + "-" * 50)
print("🚀 All sites sent sparse payloads to central server")
print("📥 Central server assembled coefficient matrix (zero-padded):")
print(np.round(beta_hat_lower, 3))
print("-" * 50, "\n")

# ==========================================
# 2. Algorithm 1: Central Site FDR Control
# ==========================================
print("--- Algorithm 1: Central Site Generalized Mirror Statistics ---")
I_hat = np.zeros(p)
pair_count = 0
all_pairs = list(combinations(range(K), 2))

for s, t in all_pairs:
    pair_count += 1
    print(f"\n=== Processing pair (Site {s+1}, Site {t+1}) ===")
    
    # Paper formula (2): Generalized mirror statistic
    M_st = np.sign(beta_hat_lower[s] * beta_hat_lower[t]) * (np.abs(beta_hat_lower[s]) + np.abs(beta_hat_lower[t]))
    print(f"Mirror statistics M: {np.round(M_st, 3)}")
    
    positive_M = M_st[M_st > 0]
    if len(positive_M) == 0:
        print(f"    ↳ No positive mirror statistics, selected set is empty")
        continue
    
    candidate_taus = np.sort(positive_M) # Sort ascending to find the smallest valid τ
    tau_alpha = np.inf
    selected_S_st = np.array([])
    
    # Paper formula (3): Find the smallest τ such that estimated FDP ≤ α
    for tau in candidate_taus:
        num_neg = np.sum(M_st < -tau)
        num_pos = max(np.sum(M_st > tau), 1)
        fdp_est = num_neg / num_pos
        
        if fdp_est <= alpha:
            tau_alpha = tau
            selected_S_st = np.where(M_st > tau)[0]
            break
    
    print(f"    ↳ Selected threshold τ = {tau_alpha:.3f}")
    print(f"    ↳ FDP estimate: {np.sum(M_st < -tau_alpha)} / {max(np.sum(M_st > tau_alpha), 1)} = {fdp_est:.3f} ≤ {alpha}")
    print(f"    ↳ Selected features this round: {selected_S_st}")
    
    # Calculate inclusion rate scores
    num_selected = len(selected_S_st)
    if num_selected > 0:
        score = 1.0 / num_selected
        print(f"    ↳ Each selected feature gets {score:.3f} points")
        for j in selected_S_st:
            I_hat[j] += score
    else:
        print(f"    ↳ No features selected, no points allocated")

# Paper formula: Normalize inclusion rates (2/(K(K-1)) = 1/pair_count)
I_hat = I_hat / pair_count
print("\n" + "-" * 50)
print("--- Algorithm 1 Final Step: Inclusion Rate Truncation ---")
print(f"Raw inclusion rates: {np.round(I_hat, 3)}")

# Step 8: Sort inclusion rates in ascending order
sorted_indices = np.argsort(I_hat)
sorted_I = I_hat[sorted_indices]
print(f"\n1. Ascending sorted inclusion rates: {np.round(sorted_I, 3)}")
print(f"   Corresponding feature indices: {sorted_indices}")

# Step 9: Find the largest m such that the sum of the first m rates ≤ α
m = 0
current_sum = 0.0
print(f"\n2. Finding maximum m with sum of first m rates ≤ {alpha}:")
for i, val in enumerate(sorted_I):
    if current_sum + val <= alpha:
        current_sum += val
        m += 1
        print(f"   - Add {val:.3f}, sum = {current_sum:.3f} ≤ {alpha}, m = {m}")
    else:
        print(f"   - Adding {val:.3f} would make sum = {current_sum + val:.3f} > {alpha}, stop")
        break

# Step 10: Final feature selection
print("\n3. Final selection:")
if m == 0:
    final_S = np.array([])
    print(f"   No features meet FDR control requirement")
else:
    threshold_I = sorted_I[m-1]
    final_S = np.where(I_hat > threshold_I)[0]
    print(f"   Threshold = {threshold_I:.3f} (the {m}-th smallest rate)")
    print(f"   Selected features (I_hat > {threshold_I:.3f}): {final_S}")

# Result evaluation
print("\n" + "=" * 50)
print("=== Final Result Evaluation ===")
print(f"True support: {sorted(list(true_support))}")
print(f"Fed-FDR selected: {sorted(list(final_S))}")

tp = len(set(final_S) & true_support)
fp = len(set(final_S) - true_support)
fdp = fp / max(len(final_S), 1)
power = tp / len(true_support)

print(f"\nTrue Positives (TP): {tp}")
print(f"False Positives (FP): {fp}")
print(f"Empirical FDP: {fdp:.3f} (target ≤ {alpha})")
print(f"Empirical Power: {power:.3f}")
print("=" * 50)