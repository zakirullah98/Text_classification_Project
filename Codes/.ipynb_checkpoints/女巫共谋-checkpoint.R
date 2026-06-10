# =====================================================================
# 终极测试：跨站点共谋攻击 (Collusion Attack) 击穿防线
# =====================================================================
library(glmnet)
library(mvtnorm)
library(ncvreg)

set.seed(99) # 换个新种子

# 1. 基础函数与干净的数据生成
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
# =====================================================================

# 找到站点 2 的前 3 个人
site2_indices <- which(site == 2)[1:3] 
# 找到站点 3 的前 2 个人
site3_indices <- which(site == 3)[1:2] 

# 把这 5 个人的索引合并成一个“作案团伙”
target_group_indices <- c(site2_indices, site3_indices)

# 让这 5 个人在第 50 号特征上都表现出极端的 15
Xall[target_group_indices, 50] <- 15        
Treatall[target_group_indices] <- 1         

cat("已成功将 5 个极端异常病人分摊：3 个在站点 2，2 个在站点 3。\n")
cat("正在运行 Fed-FDR 算法收集各站点模型参数 (请稍候)...\n\n")

# 运行算法获取中心服务器的 Theta 矩阵
local_res <- local.glm.lasso(Treatall = Treatall, Xall = Xall, site = site, weight = FALSE)
Theta_matrix <- local_res$Theta 

# =====================================================================
# 🔪 检阅攻击成果
# =====================================================================
cat("========================================================\n")
cat(" 终极攻击结果：防线是否被撕裂？\n")
cat("========================================================\n")

feature_50_coefficients <- Theta_matrix[, 50]
cat("各站点在第 50 号特征上的系数估计值：\n")
for(k in 1:K){
  cat(sprintf("站点 %d: %.4f\n", k, feature_50_coefficients[k]))
}

# 找出系数绝对值最大的两个站点（因为我们分摊到了两个站点）
# order(..., decreasing = TRUE) 会返回从大到小的索引排名
top_2_sites <- order(abs(feature_50_coefficients), decreasing = TRUE)[1:2]

cat(sprintf("\n[攻击结论] 算法受到共谋欺骗！系数绝对值最大的两个站点被精准锁定为：站点 %d 和 站点 %d。\n", 
            top_2_sites[1], top_2_sites[2]))

# 核对事实：只要排名前二的站点刚好是 2 和 3，攻击就完美成功
  if(all(sort(top_2_sites) == c(2, 3))) {
  cat("[事实核对] 目标团伙确实分布在站点 2 和 3。攻击完美成功！🎯🎯🎯\n\n")
  cat("深度解析：\n")
  cat("1. 当算法计算站点 2 时，排除了站点 2 的特征池。但站点 3 举手说：‘第 50 号特征很重要！’\n")
  cat("2. 于是算法强迫站点 2 拟合第 50 号特征。站点 2 内部的 3 个异常病人瞬间把系数拉爆。\n")
  cat("3. 同理，当计算站点 3 时，站点 2 的‘共谋’也帮它越过了防火墙。\n")
  cat("这证明了：基于交叉验证的统计防御，在面对多节点共谋（Sybil Attack）时，是不堪一击的！\n")
} else {
  cat("[事实核对] 攻击出现意外。请检查数据或参数。\n")
}