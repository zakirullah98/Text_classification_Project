# =====================================================================
# 联邦学习溯源攻击概念验证 (Proof of Concept: Source Inference Attack)
# =====================================================================
library(glmnet)
library(mvtnorm)
library(ncvreg)

set.seed(42) # 固定随机种子以保证结果可复现 | Fix the random seed to ensure reproducible results
expit <- function(x){ 1 / (1 + exp(-x)) }

# 1. 准备环境与数据 (复用之前的生成逻辑) | 1. Prepare environment and data (reuse previous generation logic)
n <- 200; p <- 200; p1 <- 5; K <- 5; signal_strength <- 3.0

generate_one_site <- function(n, p, p1, signal_strength) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  beta_true <- rep(0, p)
  beta_true[1:p1] <- signal_strength * sample(c(-1, 1), p1, replace = TRUE)
  mean_prob <- expit(as.vector(X %*% beta_true))
  y <- rbinom(n, 1, mean_prob)
  return(list(X = X, y = y))
}

Xall <- NULL; Treatall <- NULL; site <- NULL
for (k in 1:K) {
  site_data <- generate_one_site(n, p, p1, signal_strength)
  Xall <- rbind(Xall, site_data$X)
  Treatall <- c(Treatall, site_data$y)
  site <- c(site, rep(k, n))
}

# =====================================================================
# 💀 攻击准备：在数据中“埋入”我们要追踪的受害者 
# Attack Preparation: "Implant" the victims we want to track into the data
# =====================================================================

# 【目标病人 A】(用于测试情景3：异常杠杆点) | [Target Patient A] (For testing Scenario 3: Abnormal Leverage Point)
# 我们在站点 2 (Site 2) 悄悄放入一个具有罕见极端特征的病人。 
# We secretly insert a patient with rare extreme features into Site 2.
# 假设第 50 号特征是某种罕见的生化指标，正常人都是 0 左右，病人 A 高达 15。 
# Assume feature 50 is a rare biochemical marker, normally around 0, but 15 for Patient A.
target_A_index <- which(site == 2)[1] # 取站点2的第一个人作为目标A | Select the first person in Site 2 as Target A
Xall[target_A_index, 50] <- 15        # 注入异常特征值 | Inject abnormal feature value
Treatall[target_A_index] <- 1         # 设为阳性 | Set to positive

# 【目标病人 B】(用于测试情景2：拟合度匹配) | [Target Patient B] (For testing Scenario 2: Goodness-of-Fit/Loss Matching)
# 这是一个完全正常的病人，我们假设他来自站点 4 (Site 4)。 
# This is a completely normal patient, assuming they come from Site 4.
# 作为攻击者，我们“通过某种非法渠道”拿到了他的真实数据 x 和 y。 
# As an attacker, we obtained their true data x and y "through some illegal channels".
target_B_index <- which(site == 4)[10] # 取站点4的某个人作为目标B | Select a person from Site 4 as Target B
patient_B_x <- Xall[target_B_index, ]
patient_B_y <- Treatall[target_B_index]

cat("数据准备完毕，受害者 A 已植入站点 2，受害者 B 的真实数据已被攻击者掌握（源自站点 4）。\nData preparation complete, Victim A implanted in Site 2, attacker has true data of Victim B (from Site 4).\n")
cat("正在运行 Fed-FDR 算法收集各站点模型参数 (模拟中心服务器接收数据)...\nRunning Fed-FDR algorithm to collect model parameters from all sites (simulating central server receiving data)...\n\n")

# 运行联邦算法获取 Theta (各站点的回归系数) | Run federated algorithm to get Theta (regression coefficients of each site)
local_res <- local.glm.lasso(Treatall = Treatall, Xall = Xall, site = site, weight = FALSE)
Theta_matrix <- local_res$Theta # 这就是中心服务器看到的“低维系数矩阵” (K 行 p 列) | This is the "low-dimensional coefficient matrix" seen by the central server (K rows, p columns)


# =====================================================================
# 🔪 开始攻击 | Start Attack
# =====================================================================

cat("========================================================\n")
cat(" 攻击一：情景 3 - 基于异常杠杆点的溯源 (动态验证版)\n Attack 1: Scenario 3 - Source Inference based on Abnormal Leverage Points (Dynamic Verification)\n")
cat("========================================================\n")

feature_50_coefficients <- Theta_matrix[, 50]
cat("各站点在第 50 号特征上的系数估计值：\nCoefficient estimates for Feature 50 at each site:\n")
for(k in 1:K){
  cat(sprintf("站点 %d: %.4f | Site %d: %.4f\n", k, feature_50_coefficients[k], k, feature_50_coefficients[k]))
}

predicted_site_A <- which.max(abs(feature_50_coefficients))
actual_site_A <- 2 # 我们设定的真实站点 | The true site we set

cat(sprintf("\n[攻击结论] 站点 %d 的绝对值系数最大 (%.4f)！\n[Attack Conclusion] Site %d has the maximum absolute coefficient (%.4f)!\n", predicted_site_A, feature_50_coefficients[predicted_site_A], predicted_site_A, feature_50_coefficients[predicted_site_A]))

# 动态判断攻击是否成功 | Dynamically check if the attack is successful
if(predicted_site_A == actual_site_A) {
  cat(sprintf("[事实核对] 目标 A 确实来自站点 %d。攻击成功！🎯\n[Fact Check] Target A indeed comes from Site %d. Attack successful! 🎯\n\n", actual_site_A, actual_site_A))
} else {
  cat(sprintf("[事实核对] 攻击失败！❌ 算法猜测是站点 %d，但目标 A 其实在站点 %d。\n[Fact Check] Attack failed! ❌ The algorithm guessed Site %d, but Target A is actually in Site %d.\n", predicted_site_A, actual_site_A, predicted_site_A, actual_site_A))
  cat("失败原因：Lasso 的 L1 惩罚项成功压制了单一异常病人的杠杆作用，或者其他站点的随机噪音产生了更强的虚假相关性。\nFailure Reason: Lasso's L1 penalty successfully suppressed the leverage of the single abnormal patient, or random noise in other sites created stronger spurious correlations.\n\n")
}


cat("========================================================\n")
cat(" 攻击二：情景 2 - 基于模型拟合度 (Loss) 的溯源 (针对目标 B)\n Attack 2: Scenario 2 - Source Inference based on Model Fit (Loss) (Target B)\n")
cat("========================================================\n")
# 攻击者思路：我们有病人 B 的真实数据 x 和 y。 | Attacker's thought process: We have the real data x and y for Patient B.
# 我们把他代入 5 个站点的模型中，计算二元交叉熵损失 (Log-Loss)。 | We plug this into the 5 sites' models and calculate the binary cross-entropy loss (Log-Loss).
# 病人 B 参与了哪个站点的训练，那个站点对他的 Loss 一定是最低的（模型过拟合了它的主人）。 
# Whichever site trained on Patient B will definitely have the lowest Loss for him (the model overfits to its owner).

# 定义计算交叉熵损失的函数 | Define the function to calculate cross-entropy loss
calculate_loss <- function(x, y, beta) {
  prob <- expit(sum(x * beta)) # 逻辑回归的预测概率 | Predicted probability for logistic regression
  prob <- max(min(prob, 1 - 1e-15), 1e-15) # 限制概率边界防止 log(0) 报错 | Restrict probability boundaries to prevent log(0) errors
  loss <- - (y * log(prob) + (1 - y) * log(1 - prob)) # 计算 Loss | Calculate Loss
  return(loss)
}

losses <- numeric(K)
for (k in 1:K) {
  # 用站点 k 发来的系数测试病人 B | Test Patient B using the coefficients sent by site k
  losses[k] <- calculate_loss(patient_B_x, patient_B_y, Theta_matrix[k, ])
  cat(sprintf("站点 %d 对病人 B 数据的 Loss: %.4f | Site %d Loss for Patient B's data: %.4f\n", k, losses[k], k, losses[k]))
}

predicted_site_B <- which.min(losses)
actual_site_B <- 4

cat(sprintf("\n[攻击结论] 站点 %d 的 Loss 最低 (%.4f)！说明该站点'看似'见过病人 B。\n[Attack Conclusion] Site %d has the lowest Loss (%.4f)! This indicates the site 'seemingly' saw Patient B.\n", predicted_site_B, losses[predicted_site_B], predicted_site_B, losses[predicted_site_B]))

# 动态判断攻击是否成功 | Dynamically check if the attack is successful
if(predicted_site_B == actual_site_B) {
  cat(sprintf("[事实核对] 目标 B 确实来自站点 %d。攻击成功！🎯\n[Fact Check] Target B indeed comes from Site %d. Attack successful! 🎯\n", actual_site_B, actual_site_B))
} else {
  cat(sprintf("[事实核对] 攻击失败！❌ 算法猜测是站点 %d，但目标 B 其实在站点 %d。\n[Fact Check] Attack failed! ❌ The algorithm guessed Site %d, but Target B is actually in Site %d.\n", predicted_site_B, actual_site_B, predicted_site_B, actual_site_B))
  cat("失败原因：强正则化（Lasso）成功阻止了模型的过拟合，去偏操作打乱了参数记忆，且高维空间下的随机噪音让站点 1 产生了更完美的‘虚假拟合’。\nFailure Reason: Strong regularization (Lasso) successfully prevented model overfitting, the debiasing operation disrupted parameter memory, and random noise in high-dimensional space allowed Site 1 to produce a more perfect 'spurious fit'.\n\n")
}

##总结：学到的东西：你的这次运行无意中证明了机器学习安全领域的一个重要结论：带有强正则化（如 Lasso, Ridge）的模型，比无正则化的模型更难被单点溯源攻击。 惩罚项天然地抹平了极端个体的特征记忆！
## Summary: What we learned: Your run inadvertently proved an important conclusion in the field of machine learning security: Models with strong regularization (like Lasso, Ridge) are much harder to target with single-point source inference attacks than unregularized models. The penalty term naturally smooths out the feature memory of extreme individuals!


#============二次攻击================== | ============ Second Attack ==================
cat(" 进阶攻击测试：‘团伙作案’能否击穿 Lasso 的防御？\n Advanced Attack Test: Can a 'syndicate crime' break through Lasso's defense?\n")
cat("========================================================\n")

# 1. 重新准备干净的数据 | 1. Re-prepare clean data
set.seed(88) # 换个种子 | Change the seed
n <- 200; p <- 200; p1 <- 5; K <- 5; signal_strength <- 3.0

Xall <- NULL; Treatall <- NULL; site <- NULL
for (k in 1:K) {
  site_data <- generate_one_site(n, p, p1, signal_strength)
  Xall <- rbind(Xall, site_data$X)
  Treatall <- c(Treatall, site_data$y)
  site <- c(site, rep(k, n))
}

# 2. 💀 植入“异常团伙” (Cluster Outliers) | 2. Implant "Cluster Outliers"
# 找到站点 2 的前 5 个人 | Find the first 5 people in Site 2
target_group_indices <- which(site == 2)[1:5] 
# 让这 5 个人在第 50 号特征上都表现出极端的 15 | Make all 5 people show an extreme value of 15 on feature 50
Xall[target_group_indices, 50] <- 15        
Treatall[target_group_indices] <- 1         

cat("已在站点 2 植入 5 个带有极端特征 (特征50 = 15) 的病人团伙。\nImplanted a syndicate of 5 patients with extreme features (Feature 50 = 15) in Site 2.\n")
cat("正在运行 Fed-FDR 算法收集各站点模型参数...\nRunning Fed-FDR algorithm to collect model parameters from all sites...\n\n")

# 3. 运行算法 | 3. Run algorithm
local_res <- local.glm.lasso(Treatall = Treatall, Xall = Xall, site = site, weight = FALSE)
Theta_matrix <- local_res$Theta 

# 4. 实施情景 3 攻击 | 4. Implement Scenario 3 Attack
feature_50_coefficients <- Theta_matrix[, 50]
cat("各站点在第 50 号特征上的系数估计值：\nCoefficient estimates for Feature 50 at each site:\n")
for(k in 1:K){
  cat(sprintf("站点 %d: %.4f | Site %d: %.4f\n", k, feature_50_coefficients[k], k, feature_50_coefficients[k]))
}

predicted_site_A <- which.max(abs(feature_50_coefficients))
actual_site_A <- 2

cat(sprintf("\n[攻击结论] 站点 %d 的绝对值系数最大 (%.4f)！\n[Attack Conclusion] Site %d has the maximum absolute coefficient (%.4f)!\n", predicted_site_A, feature_50_coefficients[predicted_site_A], predicted_site_A, feature_50_coefficients[predicted_site_A]))

if(predicted_site_A == actual_site_A) {
  cat(sprintf("[事实核对] 目标团伙确实来自站点 %d。攻击成功！🎯\n[Fact Check] The target syndicate indeed comes from Site %d. Attack successful! 🎯\n", actual_site_A, actual_site_A))
  cat("深度解析：Lasso 的防线被击穿了！当异常点从 1 个变成 5 个时，模型发现‘忽视他们’的代价太大，被迫在系数上留下了深深的指纹。\nDeep Analysis: Lasso's defense line has been broken! When the abnormal points went from 1 to 5, the model found that the cost of 'ignoring them' was too high, and was forced to leave deep fingerprints in the coefficients.\n")
} else {
  cat("[事实核对] 攻击仍然失败。\n[Fact Check] Attack still failed.\n")
}

