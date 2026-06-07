# =====================================================================
# 联邦学习溯源攻击概念验证 (Proof of Concept: Source Inference Attack)
# =====================================================================
library(glmnet)
library(mvtnorm)
library(ncvreg)

set.seed(42) # 固定随机种子以保证结果可复现
expit <- function(x){ 1 / (1 + exp(-x)) }

# 1. 准备环境与数据 (复用之前的生成逻辑)
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
# =====================================================================

# 【目标病人 A】(用于测试情景3：异常杠杆点)
# 我们在站点 2 (Site 2) 悄悄放入一个具有罕见极端特征的病人。
# 假设第 50 号特征是某种罕见的生化指标，正常人都是 0 左右，病人 A 高达 15。
target_A_index <- which(site == 2)[1] # 取站点2的第一个人作为目标A
Xall[target_A_index, 50] <- 15        # 注入异常特征值
Treatall[target_A_index] <- 1         # 设为阳性

# 【目标病人 B】(用于测试情景2：拟合度匹配)
# 这是一个完全正常的病人，我们假设他来自站点 4 (Site 4)。
# 作为攻击者，我们“通过某种非法渠道”拿到了他的真实数据 x 和 y。
target_B_index <- which(site == 4)[10] # 取站点4的某个人作为目标B
patient_B_x <- Xall[target_B_index, ]
patient_B_y <- Treatall[target_B_index]

cat("数据准备完毕，受害者 A 已植入站点 2，受害者 B 的真实数据已被攻击者掌握（源自站点 4）。\n")
cat("正在运行 Fed-FDR 算法收集各站点模型参数 (模拟中心服务器接收数据)...\n\n")

# 运行联邦算法获取 Theta (各站点的回归系数)
local_res <- local.glm.lasso(Treatall = Treatall, Xall = Xall, site = site, weight = FALSE)
Theta_matrix <- local_res$Theta # 这就是中心服务器看到的“低维系数矩阵” (K 行 p 列)


# =====================================================================
# 🔪 开始攻击
# =====================================================================

cat("========================================================\n")
cat(" 攻击一：情景 3 - 基于异常杠杆点的溯源 (动态验证版)\n")
cat("========================================================\n")

feature_50_coefficients <- Theta_matrix[, 50]
cat("各站点在第 50 号特征上的系数估计值：\n")
for(k in 1:K){
  cat(sprintf("站点 %d: %.4f\n", k, feature_50_coefficients[k]))
}

predicted_site_A <- which.max(abs(feature_50_coefficients))
actual_site_A <- 2 # 我们设定的真实站点

cat(sprintf("\n[攻击结论] 站点 %d 的绝对值系数最大 (%.4f)！\n", predicted_site_A, feature_50_coefficients[predicted_site_A]))

# 动态判断攻击是否成功
if(predicted_site_A == actual_site_A) {
  cat(sprintf("[事实核对] 目标 A 确实来自站点 %d。攻击成功！🎯\n\n", actual_site_A))
} else {
  cat(sprintf("[事实核对] 攻击失败！❌ 算法猜测是站点 %d，但目标 A 其实在站点 %d。\n", predicted_site_A, actual_site_A))
  cat("失败原因：Lasso 的 L1 惩罚项成功压制了单一异常病人的杠杆作用，或者其他站点的随机噪音产生了更强的虚假相关性。\n\n")
}


cat("========================================================\n")
cat(" 攻击二：情景 2 - 基于模型拟合度 (Loss) 的溯源 (针对目标 B)\n")
cat("========================================================\n")
# 攻击者思路：我们有病人 B 的真实数据 x 和 y。
# 我们把他代入 5 个站点的模型中，计算二元交叉熵损失 (Log-Loss)。
# 病人 B 参与了哪个站点的训练，那个站点对他的 Loss 一定是最低的（模型过拟合了它的主人）。

# 定义计算交叉熵损失的函数
calculate_loss <- function(x, y, beta) {
  # 逻辑回归的预测概率
  prob <- expit(sum(x * beta))
  # 限制概率边界防止 log(0) 报错
  prob <- max(min(prob, 1 - 1e-15), 1e-15)
  # 计算 Loss
  loss <- - (y * log(prob) + (1 - y) * log(1 - prob))
  return(loss)
}

losses <- numeric(K)
for (k in 1:K) {
  # 用站点 k 发来的系数测试病人 B
  losses[k] <- calculate_loss(patient_B_x, patient_B_y, Theta_matrix[k, ])
  cat(sprintf("站点 %d 对病人 B 数据的 Loss: %.4f\n", k, losses[k]))
}

predicted_site_B <- which.min(losses)
actual_site_B <- 4

cat(sprintf("\n[攻击结论] 站点 %d 的 Loss 最低 (%.4f)！说明该站点'看似'见过病人 B。\n", predicted_site_B, losses[predicted_site_B]))

# 动态判断攻击是否成功
if(predicted_site_B == actual_site_B) {
  cat(sprintf("[事实核对] 目标 B 确实来自站点 %d。攻击成功！🎯\n", actual_site_B))
} else {
  cat(sprintf("[事实核对] 攻击失败！❌ 算法猜测是站点 %d，但目标 B 其实在站点 %d。\n", predicted_site_B, actual_site_B))
  cat("失败原因：强正则化（Lasso）成功阻止了模型的过拟合，去偏操作打乱了参数记忆，且高维空间下的随机噪音让站点 1 产生了更完美的‘虚假拟合’。\n\n")
}




##总结：学到的东西：你的这次运行无意中证明了机器学习安全领域的一个重要结论：带有强正则化（如 Lasso, Ridge）的模型，比无正则化的模型更难被单点溯源攻击。 惩罚项天然地抹平了极端个体的特征记忆！


#============二次攻击==================
cat(" 进阶攻击测试：‘团伙作案’能否击穿 Lasso 的防御？\n")
cat("========================================================\n")

# 1. 重新准备干净的数据
set.seed(88) # 换个种子
n <- 200; p <- 200; p1 <- 5; K <- 5; signal_strength <- 3.0

Xall <- NULL; Treatall <- NULL; site <- NULL
for (k in 1:K) {
  site_data <- generate_one_site(n, p, p1, signal_strength)
  Xall <- rbind(Xall, site_data$X)
  Treatall <- c(Treatall, site_data$y)
  site <- c(site, rep(k, n))
}

# 2. 💀 植入“异常团伙” (Cluster Outliers)
# 找到站点 2 的前 5 个人
target_group_indices <- which(site == 2)[1:5] 
# 让这 5 个人在第 50 号特征上都表现出极端的 15
Xall[target_group_indices, 50] <- 15        
Treatall[target_group_indices] <- 1         

cat("已在站点 2 植入 5 个带有极端特征 (特征50 = 15) 的病人团伙。\n")
cat("正在运行 Fed-FDR 算法收集各站点模型参数...\n\n")

# 3. 运行算法
local_res <- local.glm.lasso(Treatall = Treatall, Xall = Xall, site = site, weight = FALSE)
Theta_matrix <- local_res$Theta 

# 4. 实施情景 3 攻击
feature_50_coefficients <- Theta_matrix[, 50]
cat("各站点在第 50 号特征上的系数估计值：\n")
for(k in 1:K){
  cat(sprintf("站点 %d: %.4f\n", k, feature_50_coefficients[k]))
}

predicted_site_A <- which.max(abs(feature_50_coefficients))
actual_site_A <- 2

cat(sprintf("\n[攻击结论] 站点 %d 的绝对值系数最大 (%.4f)！\n", predicted_site_A, feature_50_coefficients[predicted_site_A]))

if(predicted_site_A == actual_site_A) {
  cat(sprintf("[事实核对] 目标团伙确实来自站点 %d。攻击成功！🎯\n", actual_site_A))
  cat("深度解析：Lasso 的防线被击穿了！当异常点从 1 个变成 5 个时，模型发现‘忽视他们’的代价太大，被迫在系数上留下了深深的指纹。\n")
} else {
  cat("[事实核对] 攻击仍然失败。\n")
}