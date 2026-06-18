# =====================================================================
# 联邦学习溯源攻击 批量测试 | 随机种子版
# =====================================================================
library(glmnet)
library(mvtnorm)
library(ncvreg)
library(pbapply)
library(parallel)

# -------------------------- 🔧 参数配置区 --------------------------
total_runs <- 10000        # 攻击次数，先100次试水，没问题改10000
GLOBAL_SEED <- 20260614  # 全局主种子，固定则结果可复现
USE_PARALLEL <- TRUE     # 是否启用并行
NUM_CORES <- detectCores() # 自动检测CPU核心数

# 模型参数
n <-100; p <- 200; p1 <- 5; K <-2; signal_strength <- 3.0
# -------------------------------------------------------------------

# 生成随机不重复的种子池
set.seed(GLOBAL_SEED)
seeds <- sample.int(10000, total_runs) # 从100万范围内随机抽total_runs个不重复种子

# 1. 基础工具函数
expit <- function(x){ 1 / (1 + exp(-x)) }

generate_one_site <- function(n, p, p1, signal_strength) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  beta_true <- rep(0, p)
  beta_true[1:p1] <- signal_strength * sample(c(-1, 1), p1, replace = TRUE)
  mean_prob <- expit(as.vector(X %*% beta_true))
  y <- rbinom(n, 1, mean_prob)
  return(list(X = X, y = y))
}

calculate_loss <- function(x, y, beta) {
  prob <- expit(sum(x * beta))
  prob <- max(min(prob, 1 - 1e-15), 1e-15)
  loss <- - (y * log(prob) + (1 - y) * log(1 - prob))
  return(loss)
}

# =====================================================================
# 🎯 单次溯源攻击函数（种子直接传入，无全局依赖）
# =====================================================================
single_source_attack <- function(seed, n, p, p1, K, signal_strength) {
  
  tryCatch({
    # ============= 第一部分：单点攻击 + Loss攻击 =============
    set.seed(seed)
    
    # 生成联邦数据
    Xall <- NULL; Treatall <- NULL; site <- NULL
    for (k in 1:K) {
      site_data <- generate_one_site(n, p, p1, signal_strength)
      Xall <- rbind(Xall, site_data$X)
      Treatall <- c(Treatall, site_data$y)
      site <- c(site, rep(k, n))
    }
    
    # 植入目标A：站点2第1个样本，特征50设为15
    target_A_index <- which(site == 2)[1]
    Xall[target_A_index, 50] <- 15
    Treatall[target_A_index] <- 1
    
    # 提取目标B：站点4第10个样本
    target_B_index <- which(site == 1)[10]
    patient_B_x <- Xall[target_B_index, ]
    patient_B_y <- Treatall[target_B_index]
    
    # 运行Fed-FDR
    local_res <- local.glm.lasso(Treatall = Treatall, Xall = Xall, site = site, weight = FALSE)
    Theta_matrix <- local_res$Theta
    
    # 攻击1：异常杠杆点溯源（单点）
    feature_50_coef <- Theta_matrix[, 50]
    pred_site_A1 <- which.max(abs(feature_50_coef))
    attack1_success <- (pred_site_A1 == 2)
    
    # 攻击2：Loss拟合度溯源
    losses <- numeric(K)
    for (k in 1:K) {
      losses[k] <- calculate_loss(patient_B_x, patient_B_y, Theta_matrix[k, ])
    }
    pred_site_B <- which.min(losses)
    attack2_success <- (pred_site_B == 4)
    
    # ============= 第二部分：团伙攻击 =============
    set.seed(seed + 200000) # 固定偏移量，保证和第一部分独立且无规律
    
    # 重新生成干净数据
    Xall2 <- NULL; Treatall2 <- NULL; site2 <- NULL
    for (k in 1:K) {
      site_data <- generate_one_site(n, p, p1, signal_strength)
      Xall2 <- rbind(Xall2, site_data$X)
      Treatall2 <- c(Treatall2, site_data$y)
      site2 <- c(site2, rep(k, n))
    }
    
    # 植入5人团伙
    target_group_idx <- which(site2 == 2)[1:5]
    Xall2[target_group_idx, 50] <- 15
    Treatall2[target_group_idx] <- 1
    
    # 运行Fed-FDR
    local_res2 <- local.glm.lasso(Treatall = Treatall2, Xall = Xall2, site = site2, weight = FALSE)
    Theta_matrix2 <- local_res2$Theta
    
    # 攻击3：团伙异常杠杆点溯源
    feature_50_coef2 <- Theta_matrix2[, 50]
    pred_site_A2 <- which.max(abs(feature_50_coef2))
    attack3_success <- (pred_site_A2 == 2)
    
    # 返回本次试验结果
    return(data.frame(
      Seed = seed,
      Attack1_单点杠杆 = attack1_success,
      Attack2_Loss溯源 = attack2_success,
      Attack3_团伙杠杆 = attack3_success,
      Error = NA
    ))
    
  }, error = function(e) {
    return(data.frame(
      Seed = seed,
      Attack1_单点杠杆 = NA,
      Attack2_Loss溯源 = NA,
      Attack3_团伙杠杆 = NA,
      Error = as.character(e)
    ))
  })
}

# =====================================================================
# 🚀 批量执行攻击
# =====================================================================
cat("========================================================\n")
cat(" 溯源攻击批量测试开始（随机种子模式）\n")
cat("========================================================\n")
cat(sprintf("测试次数：%d 次\n", total_runs))
cat(sprintf("全局主种子：%d\n", GLOBAL_SEED))
cat(sprintf("并行计算：%s\n", ifelse(USE_PARALLEL, paste("启用，", NUM_CORES, "核心"), "禁用")))
cat("攻击方式：1.单点杠杆溯源  2.Loss拟合溯源  3.团伙杠杆溯源\n\n")

if (USE_PARALLEL) {
  cl <- makeCluster(NUM_CORES)
  
  clusterExport(cl, c(
    "expit", "generate_one_site", "calculate_loss",
    "local.glm.lasso", "REF_DS_inf", "single_source_attack"
  ))
  
  clusterEvalQ(cl, { library(glmnet); library(mvtnorm); library(ncvreg) })
  
  results_list <- pblapply(
    seeds, single_source_attack,
    n = n, p = p, p1 = p1, K = K,
    signal_strength = signal_strength,
    cl = cl
  )
  
  stopCluster(cl)
} else {
  results_list <- pblapply(
    seeds, single_source_attack,
    n = n, p = p, p1 = p1, K = K,
    signal_strength = signal_strength
  )
}

# 合并结果
results <- do.call(rbind, results_list)
valid <- results[!is.na(results$Attack1_单点杠杆), ]

# =====================================================================
# 📊 结果统计
# =====================================================================
cat("\n\n========================================================\n")
cat(" 攻击成功率统计\n")
cat("========================================================\n")

rate1 <- mean(valid$Attack1_单点杠杆) * 100
rate2 <- mean(valid$Attack2_Loss溯源) * 100
rate3 <- mean(valid$Attack3_团伙杠杆) * 100

cat(sprintf("有效试验次数：%d / %d\n\n", nrow(valid), total_runs))

cat("1. 单点异常杠杆溯源（1个目标）：\n")
cat(sprintf("   成功率：%.2f%%\n\n", rate1))

cat("2. Loss拟合度溯源（正常样本）：\n")
cat(sprintf("   成功率：%.2f%%\n\n", rate2))

cat("3. 团伙异常杠杆溯源（5个目标）：\n")
cat(sprintf("   成功率：%.2f%%\n\n", rate3))

cat("========================================================\n")
cat(" 结果解读\n")
cat("========================================================\n")
cat("• 单点杠杆攻击成功率低 → Lasso正则化有效压制了单个异常点的指纹\n")
cat("• 团伙攻击成功率显著上升 → 异常样本越多，模型留下的痕迹越明显\n")
cat("• Loss溯源成功率接近随机 → 模型过拟合程度低，个体记忆弱\n\n")

# =====================================================================
# 💾 保存结果
# =====================================================================
write.csv(results, "source_attack_results.csv", row.names = FALSE)
saveRDS(results, "source_attack_results.rds")

cat("结果已保存：source_attack_results.csv / source_attack_results.rds\n")
cat("验证通过后，将 total_runs 改为 10000 即可运行完整测试\n")

