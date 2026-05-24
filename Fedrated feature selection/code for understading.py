import numpy as np
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')
# ==========================================
# 0. 数据集生成 (对应论文 Section 5 Simulation)
# ==========================================
np.random.seed(42)

K = 3          # 站点数量
p = 5          # 特征维度
n = 200        # 每个站点的样本量
alpha = 0.2    # 目标 FDR 控制水平 (玩具模型中为了展示效果适度放宽)

# 设定真实的系数：前两个特征为真实信号，后三个为纯噪声 (Null features)
true_beta = np.array([2.0, -1.5, 0.0, 0.0, 0.0])

# 生成分布在 3 个站点的本地数据
data_sites = []
for k in range(K):
    # 生成服从标准正态分布的特征矩阵
    X_k = np.random.randn(n, p)
    # 生成响应变量 y (Logistic GLM)
    logits = X_k @ true_beta
    probs = 1 / (1 + np.exp(-logits))
    y_k = np.random.binomial(1, probs)
    data_sites.append((X_k, y_k))
# ==========================================
# 补充：偷偷看一眼生成的数据长什么样
# ==========================================
print("\n--- 正在查看 [站点 1] 的数据 ---")
X_0, y_0 = data_sites[0] # 取出第一个站点的数据 (Python 索引从0开始)

print(f"站点 1 共有 {X_0.shape[0]} 个样本，{X_0.shape[1]} 个特征。")
print("\n前 5 个人的特征矩阵 X (只看前5行):")
print(np.round(X_0[:5], 2)) # 保留两位小数，方便查看

print("\n前 5 个人的诊断结果标签 y (0代表阴性，1代表阳性):")
print(y_0[:5])
print("-" * 40, "\n")
# ==========================================
# ==========================================
# 1. 算法 2 前半部分：本地站点的任务 (Local Site Task)
# ==========================================
print("--- Algorithm 2: 本地站点支持集提取与交叉 ---")
local_supports = []

# Step 1: 每个站点先用 GLM Lasso 筛选出初始支持集 S^(k)
for k in range(K):
    X_k, y_k = data_sites[k]
    # 使用带 5 折交叉验证的 L1 正则化 Logistic 回归
    lasso = LogisticRegressionCV(cv=5, penalty='l1', solver='liblinear', max_iter=1000)
    lasso.fit(X_k, y_k)
    # 提取非零特征的索引作为初始支持集
    support_k = set(np.where(lasso.coef_[0] != 0)[0])
    local_supports.append(support_k)
    print(f"站点 {k+1} 的初始 Lasso 支持集: {support_k}")
# ================= 填入这里：教学演示补丁 =================
print("\n(注：为了演示降维和补零机制，此处人为设定了初始稀疏支持集)")
print("\n(注：为了真实演示多站点下的噪声冲突与寻阈机制，此处设定了具有分歧的初始支持集)")
# 信号特征 [0, 1] 大家都找到了
# 噪声特征 [2, 3, 4] 在各个站点随机出现，方向冲突
local_supports = [
    {0, 1, 2},    # 站点 1 误选了特征 2
    {0, 1, 3},    # 站点 2 误选了特征 3
    {0, 1, 4}     # 站点 3 误选了特征 4
]
# =======================================================
beta_hat_lower = np.zeros((K, p))
# Step 2: 站点之间“交换”支持集，并计算降维后的去偏估计
# (对应论文中的 between-site swap strategy)
beta_hat_lower = np.zeros((K, p))

for k in range(K):
    X_k, y_k = data_sites[k]
    
    # 构建 S^(-k)：除当前站点 k 以外，其余所有站点支持集的并集
    S_minus_k = set()
    for j in range(K):
        if j != k:
            S_minus_k = S_minus_k.union(local_supports[j])
    S_minus_k = list(S_minus_k)
    print(f"站点 {k+1} 借用其他站点的并集 S^(-k): {S_minus_k}")
    # 仅使用 S^(-k) 中的特征，重新拟合无惩罚模型以获得去偏估计
    # 真实场景高维下使用 Xia et al. (2023) 的方法，此处低维直接用 MLE
    if len(S_minus_k) > 0:
        mle = LogisticRegression(penalty=None, max_iter=1000)
        mle.fit(X_k[:, S_minus_k], y_k)
        # 将估计值填入对应的维度，未选中的维度保持为 0
        beta_hat_lower[k, S_minus_k] = mle.coef_[0]
        print(f"    ↳ 本地计算结果: 针对特征 {S_minus_k} 算出的去偏系数为 {np.round(mle.coef_[0], 3)}")
        print(f"    ↳ 实际传输数据: 补零后传给中心的 5 维向量为 {np.round(beta_hat_lower[k], 3)}\n")

print("\n传输给中心站点的去偏系数矩阵 (每行代表一个站点):")
print(np.round(beta_hat_lower, 3))


# ==========================================
# 2. 算法 1：中心站点的任务 (Central Site Task)
# ==========================================
print("\n--- Algorithm 1: 中心站点基于广义镜像统计量的 FDR 控制 ---")
I_hat = np.zeros(p)
pair_count = 0

# 遍历所有站点的两两组合 (s, t)
for s, t in combinations(range(K), 2):
    pair_count += 1
    
    # 计算广义镜像统计量 M_j^(st)
    M_st = np.zeros(p)
    for j in range(p):
        M_st[j] = np.sign(beta_hat_lower[s, j] * beta_hat_lower[t, j]) * \
                  (abs(beta_hat_lower[s, j]) + abs(beta_hat_lower[t, j]))
    
    # 计算 data-driven 阈值 tau
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
            
    print(f"组合 (站点 {s+1}, 站点 {t+1}) -> M统计量: {np.round(M_st, 3)}")
    print(f"    ↳ 【寻阈逻辑】尝试阈值 tau = {tau_alpha:.3f}。此时 M < -{tau_alpha:.3f} 的特征有 {np.sum(M_st < -tau_alpha)} 个(推测假阳性)，M > {tau_alpha:.3f} 的有 {np.sum(M_st > tau_alpha)} 个(候选池)。")
    print(f"    ↳ 【FDP估计】{np.sum(M_st < -tau_alpha)} / max({np.sum(M_st > tau_alpha)}, 1) = {np.sum(M_st < -tau_alpha)/max(np.sum(M_st > tau_alpha), 1):.3f} <= alpha({alpha})，及格线成立！")
    print(f"    ↳ 【本轮入选】选择 M > {tau_alpha:.3f} 的特征集合: {selected_S_st}")
    #候选的阈值 $\tau$ 是取所有非零 M 统计量的绝对值，并从小到大排序。
    # 对于每个候选阈值 $\tau$，我们计算满足 $M_j < -\tau$ 的特征数量（推测的假阳性）和满足 $M_j > \tau$ 的特征数量（候选池）。
    # 然后估计 FDP = (假阳性数量) / max(候选池数量, 1)。我们选择第一个使得 FDP <= alpha 的 $\tau$ 作为最终的阈值，并选出对应的特征集合。
    
    # 累加得分
    num_selected = len(selected_S_st)
    if num_selected > 0:
        score = 1.0 / num_selected
        print(f"    ↳ 【得分分配】本轮选出 {num_selected} 个特征，每个获得 1/{num_selected} = {score:.3f} 分。\n")
        for j in selected_S_st:
            I_hat[j] += score
    else:
        print(f"    ↳ 【得分分配】本轮无特征入选。\n")

print("--------------------------------------------------")
print("--- 算法 1 最终步：基于包含率 (I_hat) 排序截断 ---")

# 标准化包含率 I_hat (修复 Bug：保证全局只除一次 pair_count！)
I_hat = I_hat / pair_count
print(f"每个特征的最终包含率 I_hat: {np.round(I_hat, 3)}")

# 对应论文 Algorithm 1 第 8 步：升序排列 (从小到大)
sorted_indices = np.argsort(I_hat)
sorted_I = I_hat[sorted_indices]
print(f"1. 升序排列包含率: {np.round(sorted_I, 3)} (对应原始特征索引: {sorted_indices})")

# 对应论文 Algorithm 1 第 9 步：寻找最大的 m，使前 m 项之和 <= alpha
m = 0
current_sum = 0
print(f"2. 寻找最大 m，使前 m 个极小包含率之和 <= alpha ({alpha}):")
for i, val in enumerate(sorted_I):
    if current_sum + val <= alpha:
        current_sum += val
        m += 1
        print(f"   - 累加第 {i+1} 个值 ({val:.3f})，当前和 = {current_sum:.3f} (未超 {alpha})，m 推进为 {m}")
    else:
        print(f"   - 尝试累加第 {i+1} 个值 ({val:.3f})，和将变为 {(current_sum + val):.3f} (超出 {alpha}！)。停止累加，最终 m = {m}")
        break

# 对应论文 Algorithm 1 第 10 步：获取最终支持集
threshold_I = sorted_I[m-1] if m > 0 else -1
final_S = np.where(I_hat > threshold_I)[0]

print(f"3. 确定截断线：排序后的第 {m} 个值 = {threshold_I:.3f}")
print(f"4. 严格剔除包含率 <= {threshold_I:.3f} 的特征，最终特征支持集: {final_S}")
print(f"真实特征集应为:{final_S}")
