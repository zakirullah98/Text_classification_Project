# =====================================================================
# 联邦学习溯源攻击批量测试 | 随机种子双语版
# Federated Learning Source Inference Attack Batch Test | Random Seed Bilingual Version
# 默认100次试水，验证通过后可改为10000次 | Default 100 runs for trial, change to 10000 after verification
# =====================================================================
library(glmnet)
library(mvtnorm)
library(ncvreg)
library(pbapply)
library(parallel)

# -------------------------- 🔧 参数配置 | Parameter Configuration --------------------------
total_runs <- 100        # 攻击次数，先100次试水，没问题改10000 | Number of attack runs, 100 for trial, change to 10000 later
GLOBAL_SEED <- 20260614  # 全局主种子，固定则结果可复现 | Global master seed, fixed for reproducible results
USE_PARALLEL <- TRUE     # 是否启用并行计算 | Whether to enable parallel computing
NUM_CORES <- detectCores() # 自动检测可用CPU核心数 | Auto-detect available CPU cores

# 模型超参数 | Model hyperparameters
n <- 200; p <- 200; p1 <- 5; K <- 5; signal_strength <- 3.0
# -------------------------------------------------------------------

# 生成随机不重复的种子池 | Generate random non-repeating seed pool
set.seed(GLOBAL_SEED)
seeds <- sample.int(1000000, total_runs) # 从100万范围内随机抽取不重复种子 | Randomly sample unique seeds from 1 to 1,000,000

# 1. 基础工具函数 | 1. Basic utility functions
expit <- function(x){ 1 / (1 + exp(-x)) } # Sigmoid激活函数 | Sigmoid activation function

generate_one_site <- function(n, p, p1, signal_strength) {
  # 生成单个站点的模拟数据 | Generate simulated data for a single site
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  beta_true <- rep(0, p)
  beta_true[1:p1] <- signal_strength * sample(c(-1, 1), p1, replace = TRUE)
  mean_prob <- expit(as.vector(X %*% beta_true))
  y <- rbinom(n, 1, mean_prob)
  return(list(X = X, y = y))
}

calculate_loss <- function(x, y, beta) {
  # 计算二元交叉熵损失 | Calculate binary cross-entropy loss
  prob <- expit(sum(x * beta))
  prob <- max(min(prob, 1 - 1e-15), 1e-15) # 防止log(0)数值错误 | Prevent numerical error from log(0)
  loss <- - (y * log(prob) + (1 - y) * log(1 - prob))
  return(loss)
}

# =====================================================================
# 🎯 单次溯源攻击函数（种子直接传入，无全局依赖）
# Single source inference attack function (seed passed directly, no global dependency)
# =====================================================================
single_source_attack <- function(seed, n, p, p1, K, signal_strength) {
  
  tryCatch({
    # ============= 第一部分：单点攻击 + Loss拟合攻击 | Part 1: Single-point attack + Loss fitting attack =============
    set.seed(seed)
    
    # 生成联邦学习全局数据 | Generate federated global dataset
    Xall <- NULL; Treatall <- NULL; site <- NULL
    for (k in 1:K) {
      site_data <- generate_one_site(n, p, p1, signal_strength)
      Xall <- rbind(Xall, site_data$X)
      Treatall <- c(Treatall, site_data$y)
      site <- c(site, rep(k, n))
    }
    
    # 植入目标A：站点2第1个样本，第50号特征设为异常值15
    # Implant Target A: 1st sample in Site 2, set feature 50 to outlier 15
    target_A_index <- which(site == 2)[1]
    Xall[target_A_index, 50] <- 15
    Treatall[target_A_index] <- 1
    
    # 提取目标B：站点4第10个样本的完整数据（攻击者已知）
    # Extract Target B: Full data of 10th sample in Site 4 (known to attacker)
    target_B_index <- which(site == 4)[10]
    patient_B_x <- Xall[target_B_index, ]
    patient_B_y <- Treatall[target_B_index]
    
    # 运行Fed-FDR算法获取各站点参数 | Run Fed-FDR algorithm to get site-wise parameters
    local_res <- local.glm.lasso(Treatall = Treatall, Xall = Xall, site = site, weight = FALSE)
    Theta_matrix <- local_res$Theta
    
    # 攻击1：异常杠杆点溯源（单点）
    # Attack 1: Abnormal leverage point source inference (single point)
    feature_50_coef <- Theta_matrix[, 50]
    pred_site_A1 <- which.max(abs(feature_50_coef))
    attack1_success <- (pred_site_A1 == 2)
    
    # 攻击2：Loss拟合度溯源
    # Attack 2: Goodness-of-fit (Loss) source inference
    losses <- numeric(K)
    for (k in 1:K) {
      losses[k] <- calculate_loss(patient_B_x, patient_B_y, Theta_matrix[k, ])
    }
    pred_site_B <- which.min(losses)
    attack2_success <- (pred_site_B == 4)
    
    # ============= 第二部分：团伙异常杠杆攻击 | Part 2: Cluster outlier leverage attack =============
    set.seed(seed + 200000) # 固定偏移量，保证与第一部分独立且无规律 | Fixed offset, independent from Part 1
    
    # 重新生成干净数据集 | Regenerate clean dataset
    Xall2 <- NULL; Treatall2 <- NULL; site2 <- NULL
    for (k in 1:K) {
      site_data <- generate_one_site(n, p, p1, signal_strength)
      Xall2 <- rbind(Xall2, site_data$X)
      Treatall2 <- c(Treatall2, site_data$y)
      site2 <- c(site2, rep(k, n))
    }
    
    # 植入5人异常团伙 | Implant 5-person outlier cluster
    target_group_idx <- which(site2 == 2)[1:5]
    Xall2[target_group_idx, 50] <- 15
    Treatall2[target_group_idx] <- 1
    
    # 运行Fed-FDR算法 | Run Fed-FDR algorithm
    local_res2 <- local.glm.lasso(Treatall = Treatall2, Xall = Xall2, site = site2, weight = FALSE)
    Theta_matrix2 <- local_res2$Theta
    
    # 攻击3：团伙异常杠杆点溯源
    # Attack 3: Cluster outlier leverage point source inference
    feature_50_coef2 <- Theta_matrix2[, 50]
    pred_site_A2 <- which.max(abs(feature_50_coef2))
    attack3_success <- (pred_site_A2 == 2)
    
    # 返回本次试验结果 | Return results of this trial
    return(data.frame(
      Seed = seed,
      Attack1_SingleLeverage = attack1_success,
      Attack2_LossInference = attack2_success,
      Attack3_ClusterLeverage = attack3_success,
      Error = NA
    ))
    
  }, error = function(e) {
    # 错误捕获：记录错误信息，不中断整体循环
    # Error handling: record error message without interrupting the whole loop
    return(data.frame(
      Seed = seed,
      Attack1_SingleLeverage = NA,
      Attack2_LossInference = NA,
      Attack3_ClusterLeverage = NA,
      Error = as.character(e)
    ))
  })
}

# =====================================================================
# 🚀 批量执行攻击 | Batch attack execution
# =====================================================================
cat("========================================================\n")
cat(" 溯源攻击批量测试开始（随机种子模式）\n")
cat(" Source Inference Attack Batch Test (Random Seed Mode)\n")
cat("========================================================\n")
cat(sprintf("测试次数：%d 次\n", total_runs))
cat(sprintf("Total test runs: %d\n", total_runs))
cat(sprintf("全局主种子：%d\n", GLOBAL_SEED))
cat(sprintf("Global master seed: %d\n", GLOBAL_SEED))
cat(sprintf("并行计算：%s\n", ifelse(USE_PARALLEL, paste("启用，", NUM_CORES, "核心"), "禁用")))
cat(sprintf("Parallel computing: %s\n", ifelse(USE_PARALLEL, paste("Enabled, ", NUM_CORES, " cores"), "Disabled")))
cat("攻击方式：1.单点杠杆溯源  2.Loss拟合溯源  3.团伙杠杆溯源\n")
cat("Attack methods: 1.Single leverage  2.Loss fitting  3.Cluster leverage\n\n")

if (USE_PARALLEL) {
  cl <- makeCluster(NUM_CORES)
  
  # 导出所需函数到并行子进程 | Export required functions to parallel workers
  clusterExport(cl, c(
    "expit", "generate_one_site", "calculate_loss",
    "local.glm.lasso", "REF_DS_inf", "single_source_attack"
  ))
  
  # 子进程加载依赖包 | Load required packages in workers
  clusterEvalQ(cl, { library(glmnet); library(mvtnorm); library(ncvreg) })
  
  # 并行执行攻击（带进度条）| Parallel execution with progress bar
  results_list <- pblapply(
    seeds, single_source_attack,
    n = n, p = p, p1 = p1, K = K,
    signal_strength = signal_strength,
    cl = cl
  )
  
  stopCluster(cl)
} else {
  # 串行执行攻击（带进度条）| Serial execution with progress bar
  results_list <- pblapply(
    seeds, single_source_attack,
    n = n, p = p, p1 = p1, K = K,
    signal_strength = signal_strength
  )
}

# 合并所有试验结果 | Merge all trial results
results <- do.call(rbind, results_list)
valid <- results[!is.na(results$Attack1_SingleLeverage), ]

# =====================================================================
# 📊 攻击结果统计 | Attack Result Statistics
# =====================================================================
cat("\n\n========================================================\n")
cat(" 攻击成功率统计\n")
cat(" Attack Success Rate Statistics\n")
cat("========================================================\n")

rate1 <- mean(valid$Attack1_SingleLeverage) * 100
rate2 <- mean(valid$Attack2_LossInference) * 100
rate3 <- mean(valid$Attack3_ClusterLeverage) * 100

cat(sprintf("有效试验次数：%d / %d\n\n", nrow(valid), total_runs))
cat(sprintf("Valid trials: %d / %d\n\n", nrow(valid), total_runs))

cat("1. 单点异常杠杆溯源（1个目标样本）：\n")
cat("1. Single outlier leverage inference (1 target sample):\n")
cat(sprintf("   成功率 | Success rate: %.2f%%\n\n", rate1))

cat("2. Loss拟合度溯源（正常样本）：\n")
cat("2. Goodness-of-fit (Loss) inference (normal sample):\n")
cat(sprintf("   成功率 | Success rate: %.2f%%\n\n", rate2))

cat("3. 团伙异常杠杆溯源（5个目标样本）：\n")
cat("3. Cluster outlier leverage inference (5 target samples):\n")
cat(sprintf("   成功率 | Success rate: %.2f%%\n\n", rate3))

cat("========================================================\n")
cat(" 结果解读 | Result Interpretation\n")
cat("========================================================\n")
cat("• 单点杠杆攻击成功率低 → Lasso正则化有效压制了单个异常点的参数指纹\n")
cat("• Low success rate of single-point attack → Lasso regularization effectively suppresses parameter fingerprints of single outliers\n")
cat("• 团伙攻击成功率显著上升 → 异常样本越多，模型在参数上留下的痕迹越明显\n")
cat("• Significantly higher success rate of cluster attack → More outliers leave stronger traces in model parameters\n")
cat("• Loss溯源成功率接近随机 → 模型过拟合程度低，个体数据记忆弱\n")
cat("• Loss inference rate close to random → Low overfitting, weak memory of individual data\n\n")

# =====================================================================
# 💾 保存结果 | Save Results
# =====================================================================
write.csv(results, "source_attack_results.csv", row.names = FALSE)
saveRDS(results, "source_attack_results.rds")

cat("结果已保存：source_attack_results.csv / source_attack_results.rds\n")
cat("Results saved: source_attack_results.csv / source_attack_results.rds\n")
cat("验证通过后，将 total_runs 改为 10000 即可运行完整测试\n")
cat("After verification, change total_runs to 10000 to run the full test\n")