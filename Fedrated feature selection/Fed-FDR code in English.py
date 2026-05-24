import numpy as np
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')
# ==========================================
# 0. Dataset Generation (Corresponds to Section 5 Simulation in the paper)
# ==========================================
np.random.seed(42)

K = 3          # Number of sites
p = 5          # Feature dimension
n = 200        # Sample size per site
alpha = 0.2    # Target FDR control level (relaxed in this toy model for demonstration)

# Set true coefficients: first two features are true signals, last three are pure noise (Null features)
true_beta = np.array([2.0, -1.5, 0.0, 0.0, 0.0])

# Generate local data distributed across 3 sites
data_sites = []
for k in range(K):
    # Generate feature matrix following standard normal distribution
    X_k = np.random.randn(n, p)
    # Generate response variable y (Logistic GLM)
    logits = X_k @ true_beta
    probs = 1 / (1 + np.exp(-logits))
    y_k = np.random.binomial(1, probs)
    data_sites.append((X_k, y_k))
# ==========================================
# Supplement: Take a quick look at the generated data
# ==========================================
print("\n--- Inspecting data for [Site 1] ---")
X_0, y_0 = data_sites[0] # Extract data for the first site (Python indices start from 0)

print(f"Site 1 has {X_0.shape[0]} samples and {X_0.shape[1]} features.")
print("\nFeature matrix X for the first 5 individuals (first 5 rows only):")
print(np.round(X_0[:5], 2)) # Keep two decimal places for easier viewing

print("\nDiagnostic outcome labels y for the first 5 individuals (0 for negative, 1 for positive):")
print(y_0[:5])
print("-" * 40, "\n")
# ==========================================
# ==========================================
# 1. Algorithm 2 First Half: Local Site Task
# ==========================================
print("--- Algorithm 2: Local Site Support Set Extraction and Swap ---")
local_supports = []

# Step 1: Each site first uses GLM Lasso to screen out the initial support set S^(k)
for k in range(K):
    X_k, y_k = data_sites[k]
    # Use L1-regularized Logistic Regression with 5-fold cross-validation
    lasso = LogisticRegressionCV(cv=5, penalty='l1', solver='liblinear', max_iter=1000)
    lasso.fit(X_k, y_k)
    # Extract indices of non-zero features as the initial support set
    support_k = set(np.where(lasso.coef_[0] != 0)[0])
    local_supports.append(support_k)
    print(f"Site {k+1} initial Lasso support set: {support_k}")
# ================= Insert here: Teaching demonstration patch =================
print("\n(Note: To demonstrate the dimensionality reduction and zero-padding mechanism, an initial sparse support set is manually set here)")
print("\n(Note: To realistically demonstrate noise conflict and threshold-finding mechanism across multiple sites, a divergent initial support set is set here)")
# Signal features [0, 1] are found by everyone
# Noise features [2, 3, 4] randomly appear at each site, causing directional conflict
local_supports = [
    {0, 1, 2},    # Site 1 mistakenly selected feature 2
    {0, 1, 3},    # Site 2 mistakenly selected feature 3
    {0, 1, 4}     # Site 3 mistakenly selected feature 4
]
# =======================================================
beta_hat_lower = np.zeros((K, p))
# Step 2: 'Swap' support sets between sites and compute debiased estimates after dimensionality reduction
# (Corresponds to the between-site swap strategy in the paper)
beta_hat_lower = np.zeros((K, p))

for k in range(K):
    X_k, y_k = data_sites[k]
    
    # Construct S^(-k): the union of support sets of all sites except the current site k
    S_minus_k = set()
    for j in range(K):
        if j != k:
            S_minus_k = S_minus_k.union(local_supports[j])
    S_minus_k = list(S_minus_k)
    print(f"Site {k+1} borrowed union from other sites S^(-k): {S_minus_k}")
    # Use only features in S^(-k) to refit the unpenalized model to obtain debiased estimates
    # In a real high-dimensional scenario, use Xia et al. (2023)'s method; here in low dimension, MLE is directly used
    if len(S_minus_k) > 0:
        mle = LogisticRegression(penalty=None, max_iter=1000)
        mle.fit(X_k[:, S_minus_k], y_k)
        # Fill the estimated values into the corresponding dimensions; unselected dimensions remain 0
        beta_hat_lower[k, S_minus_k] = mle.coef_[0]
        print(f"    ↳ Local calculation result: Debiased coefficients for features {S_minus_k} are {np.round(mle.coef_[0], 3)}")
        print(f"    ↳ Actual transmitted data: 5-dimensional vector after zero-padding sent to center is {np.round(beta_hat_lower[k], 3)}\n")

print("\nDebiased coefficient matrix transmitted to the central site (each row represents a site):")
print(np.round(beta_hat_lower, 3))


# ==========================================
# 2. Algorithm 1: Central Site Task
# ==========================================
print("\n--- Algorithm 1: Central site FDR control based on generalized mirror statistics ---")
I_hat = np.zeros(p)
pair_count = 0

# Iterate through all pairwise combinations of sites (s, t)
for s, t in combinations(range(K), 2):
    pair_count += 1
    
    # Calculate generalized mirror statistic M_j^(st)
    M_st = np.zeros(p)
    for j in range(p):
        M_st[j] = np.sign(beta_hat_lower[s, j] * beta_hat_lower[t, j]) * \
                  (abs(beta_hat_lower[s, j]) + abs(beta_hat_lower[t, j]))
    
    # Calculate data-driven threshold tau
    candidate_taus = np.sort(np.abs(M_st[M_st != 0]))
    tau_alpha = np.inf
    selected_S_st = []
    
    for tau in candidate_taus:
        num_neg = np.sum(M_st < -tau)
        num_pos = max(np.sum(M_st > tau), 1)
        fdp_est = num_neg / num_pos
        
        if fdp_est <= alpha:
            tau_alpha = tau
            selected_S_st = np.where(M_st > tau)[0]
            break
            
    print(f"Combination (Site {s+1}, Site {t+1}) -> M statistic: {np.round(M_st, 3)}")
    print(f"    ↳ [Threshold Logic] Trying threshold tau = {tau_alpha:.3f}. Features with M < -{tau_alpha:.3f} are {np.sum(M_st < -tau_alpha)} (suspected false positives), features with M > {tau_alpha:.3f} are {np.sum(M_st > tau_alpha)} (candidate pool).")
    print(f"    ↳ [FDP Estimation] {np.sum(M_st < -tau_alpha)} / max({np.sum(M_st > tau_alpha)}, 1) = {np.sum(M_st < -tau_alpha)/max(np.sum(M_st > tau_alpha), 1):.3f} <= alpha({alpha}), passing line established!")
    print(f"    ↳ [Selected this round] Selected feature set with M > {tau_alpha:.3f}: {selected_S_st}")
    # The candidate threshold \tau is taken as the absolute values of all non-zero M statistics and sorted in ascending order.
    # For each candidate threshold \tau, we calculate the number of features satisfying M_j < -\tau (suspected false positives) and M_j > \tau (candidate pool).
    # Then estimate FDP = (number of false positives) / max(candidate pool size, 1). We choose the first \tau that makes FDP <= alpha as the final threshold and select the corresponding feature set.
    
    # Accumulate score
    num_selected = len(selected_S_st)
    if num_selected > 0:
        score = 1.0 / num_selected
        print(f"    ↳ [Score Allocation] {num_selected} features selected this round, each gets 1/{num_selected} = {score:.3f} points.\n")
        for j in selected_S_st:
            I_hat[j] += score
    else:
        print(f"    ↳ [Score Allocation] No features selected this round.\n")

print("--------------------------------------------------")
print("--- Algorithm 1 Final Step: Truncation based on inclusion rate (I_hat) sorting ---")

# Standardize inclusion rate I_hat (Bug fix: ensure division by pair_count only happens once globally!)
I_hat = I_hat / pair_count
print(f"Final inclusion rate I_hat for each feature: {np.round(I_hat, 3)}")

# Corresponds to Algorithm 1 Step 8 in the paper: Ascending sort (from smallest to largest)
sorted_indices = np.argsort(I_hat)
sorted_I = I_hat[sorted_indices]
print(f"1. Ascending sort of inclusion rates: {np.round(sorted_I, 3)} (Corresponding original feature indices: {sorted_indices})")

# Corresponds to Algorithm 1 Step 9 in the paper: Find the largest m such that the sum of the first m items <= alpha
m = 0
current_sum = 0
print(f"2. Finding maximum m so that the sum of the first m minimal inclusion rates <= alpha ({alpha}):")
for i, val in enumerate(sorted_I):
    if current_sum + val <= alpha:
        current_sum += val
        m += 1
        print(f"   - Accumulating the {i+1}-th value ({val:.3f}), current sum = {current_sum:.3f} (not exceeding {alpha}), advancing m to {m}")
    else:
        print(f"   - Trying to accumulate the {i+1}-th value ({val:.3f}), sum will become {(current_sum + val):.3f} (exceeds {alpha}!). Stopping accumulation, final m = {m}")
        break

# Corresponds to Algorithm 1 Step 10 in the paper: Obtain the final support set
threshold_I = sorted_I[m-1] if m > 0 else -1
final_S = np.where(I_hat > threshold_I)[0]

print(f"3. Determine the truncation line: the {m}-th value after sorting = {threshold_I:.3f}")
print(f"4. Strictly eliminate features with inclusion rate <= {threshold_I:.3f}, final feature support set: {final_S}")
print(f"The true feature set should be: {final_S}")