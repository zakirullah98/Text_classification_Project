# =====================================================================
# 终极测试：跨站点共谋攻击 (Collusion Attack) 击穿防线
# Ultimate Test: Cross-Site Collusion Attack breaking through the defense
# =====================================================================
library(glmnet)
library(mvtnorm)
library(ncvreg)

set.seed(99) # 换个新种子 | Change to a new seed

# 1. 基础函数与干净的数据生成 | 1. Basic functions and clean data generation
expit <- function(x){ 1 / (1 + exp(-x)) }

generate_one_site <- function(n, p, p1, signal_strength) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  beta_true <- rep(0, p)
  beta_true[1:p1] <- signal_strength * sample(c(-1, 1), p1, replace = TRUE)
  mean_prob <- expit(as.vector(X %*% beta_true))
  y <- rbinom(n, 1, mean_prob)
  return(list(X = X, y = y))
}

n <- 200; p <- 200; p1 <- 5; K <- 5; signal_strength <- 3.0
Xall <- NULL; Treatall <- NULL; site <- NULL

for (k in 1:K) {
  site_data <- generate_one_site(n, p, p1, signal_strength)
  Xall <- rbind(Xall, site_data$X)
  Treatall <- c(Treatall, site_data$y)
  site <- c(site, rep(k, n))
}

# =====================================================================
# 💀 核心攻击逻辑：分摊异常病人，形成“共谋”
# Core attack logic: Distribute abnormal patients to form a "collusion"
# =====================================================================

# 找到站点 2 的前 3 个人 | Find the first 3 people in Site 2
site2_indices <- which(site == 2)[1:3] 
# 找到站点 3 的前 2 个人 | Find the first 2 people in Site 3
site3_indices <- which(site == 3)[1:2] 

# 把这 5 个人的索引合并成一个“作案团伙” | Combine the indices of these 5 people into a "criminal syndicate"
target_group_indices <- c(site2_indices, site3_indices)

# 让这 5 个人在第 50 号特征上都表现出极端的 15 | Make all 5 people show an extreme value of 15 on feature 50
Xall[target_group_indices, 50] <- 15        
Treatall[target_group_indices] <- 1         

cat("已成功将 5 个极端异常病人分摊：3 个在站点 2，2 个在站点 3。\nSuccessfully distributed 5 extreme abnormal patients: 3 in Site 2, 2 in Site 3.\n")
cat("正在运行 Fed-FDR 算法收集各站点模型参数 (请稍候)...\nRunning Fed-FDR algorithm to collect model parameters from all sites (please wait)...\n\n")

# 运行算法获取中心服务器的 Theta 矩阵 | Run the algorithm to get the Theta matrix at the central server
local_res <- local.glm.lasso(Treatall = Treatall, Xall = Xall, site = site, weight = FALSE)
Theta_matrix <- local_res$Theta 

# =====================================================================
# 🔪 检阅攻击成果 | Review Attack Results
# =====================================================================
cat("========================================================\n")
cat(" 终极攻击结果：防线是否被撕裂？\n Ultimate Attack Result: Is the defense torn apart?\n")
cat("========================================================\n")

feature_50_coefficients <- Theta_matrix[, 50]
cat("各站点在第 50 号特征上的系数估计值：\nCoefficient estimates for Feature 50 at each site:\n")
for(k in 1:K){
  cat(sprintf("站点 %d: %.4f | Site %d: %.4f\n", k, feature_50_coefficients[k], k, feature_50_coefficients[k]))
}

# 找出系数绝对值最大的两个站点（因为我们分摊到了两个站点） 
# Find the two sites with the largest absolute coefficients (because we distributed them into two sites)
# order(..., decreasing = TRUE) 会返回从大到小的索引排名 
# order(..., decreasing = TRUE) returns index rankings from largest to smallest
top_2_sites <- order(abs(feature_50_coefficients), decreasing = TRUE)[1:2]

cat(sprintf("\n[攻击结论] 算法受到共谋欺骗！系数绝对值最大的两个站点被精准锁定为：站点 %d 和 站点 %d。\n[Attack Conclusion] The algorithm is deceived by collusion! The two sites with the largest absolute coefficients are precisely locked as: Site %d and Site %d.\n", 
            top_2_sites[1], top_2_sites[2], top_2_sites[1], top_2_sites[2]))

# 核对事实：只要排名前二的站点刚好是 2 和 3，攻击就完美成功 
# Fact check: As long as the top two sites are exactly 2 and 3, the attack is perfectly successful
if(all(sort(top_2_sites) == c(2, 3))) {
  cat("[事实核对] 目标团伙确实分布在站点 2 和 3。攻击完美成功！🎯🎯🎯\n[Fact Check] The target syndicate is indeed distributed in Sites 2 and 3. Attack perfectly successful! 🎯🎯🎯\n\n")
  cat("深度解析：\nDeep Analysis:\n")
  cat("1. 当算法计算站点 2 时，排除了站点 2 的特征池。但站点 3 举手说：‘第 50 号特征很重要！’\n1. When the algorithm calculates Site 2, it excludes Site 2's feature pool. But Site 3 raises a hand and says: 'Feature 50 is very important!'\n")
  cat("2. 于是算法强迫站点 2 拟合第 50 号特征。站点 2 内部的 3 个异常病人瞬间把系数拉爆。\n2. Thus, the algorithm forces Site 2 to fit Feature 50. The 3 abnormal patients inside Site 2 instantly blow up the coefficient.\n")
  cat("3. 同理，当计算站点 3 时，站点 2 的‘共谋’也帮它越过了防火墙。\n3. Similarly, when calculating Site 3, Site 2's 'collusion' also helps it bypass the firewall.\n")
  cat("这证明了：基于交叉验证的统计防御，在面对多节点共谋（Sybil Attack）时，是不堪一击的！\nThis proves that: Statistical defense based on cross-validation is extremely vulnerable when facing multi-node collusion (Sybil Attack)!\n")
} else {
  cat("[事实核对] 攻击出现意外。请检查数据或参数。\n[Fact Check] The attack had an unexpected outcome. Please check data or parameters.\n")
}

