# =====================================================================
# 成员推断攻击 (MIA) 批量测试 | 随机种子版
# =====================================================================
library(glmnet)
library(mvtnorm)
library(ncvreg)
library(pbapply)
library(parallel)
library(ggplot2)

# -------------------------- 🔧 参数配置区 --------------------------
total_runs <- 10000        # 攻击次数，先100次试水，没问题改10000
GLOBAL_SEED <- 20260614  # 全局主种子，固定则结果可复现
USE_PARALLEL <- TRUE     # 是否启用并行
NUM_CORES <- detectCores() # 自动检测CPU核心数

# 模型参数
n <- 200; p <- 200; p1 <- 5; K <- 5; signal_strength <- 3.0
# -------------------------------------------------------------------

# 生成随机不重复的种子池
set.seed(GLOBAL_SEED)
seeds <- sample.int(1000000, total_runs)

# 基础函数
expit <- function(x) { 1 / (1 + exp(-x)) }

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
  - (y * log(prob) + (1 - y) * log(1 - prob))
}

# =====================================================================
# 🎯 单次攻击函数（种子直接传入）
# =====================================================================
single_attack <- function(seed, n, p, p1, K, signal_strength) {
  
  tryCatch({
    set.seed(seed)
    
    # 生成联邦数据
    Xall <- NULL; Treatall <- NULL; site <- NULL
    for (k in 1:K) {
      site_data <- generate_one_site(n, p, p1, signal_strength)
      Xall <- rbind(Xall, site_data$X)
      Treatall <- c(Treatall, site_data$y)
      site <- c(site, rep(k, n))
    }
    
    # 提取目标样本
    member_idx <- which(site == 2)[42]
    member_x <- Xall[member_idx, ]; member_y <- Treatall[member_idx]
    non_member_data <- generate_one_site(1, p, p1, signal_strength)
    non_member_x <- non_member_data$X[1, ]; non_member_y <- non_member_data$y[1]
    
    # 训练模型
    local_res <- local.glm.lasso(Treatall = Treatall, Xall = Xall, site = site, weight = FALSE)
    Theta_matrix <- local_res$Theta
    
    # 计算损失
    member_losses <- numeric(K)
    non_member_losses <- numeric(K)
    for (k in 1:K) {
      member_losses[k] <- calculate_loss(member_x, member_y, Theta_matrix[k, ])
      non_member_losses[k] <- calculate_loss(non_member_x, non_member_y, Theta_matrix[k, ])
    }
    
    min_loss_member <- min(member_losses)
    min_loss_non_member <- min(non_member_losses)
    
    return(data.frame(
      Seed = seed,
      MinLossMember = min_loss_member,
      MinLossNonMember = min_loss_non_member,
      Success = (min_loss_member < min_loss_non_member),
      Error = NA
    ))
  }, error = function(e) {
    return(data.frame(
      Seed = seed,
      MinLossMember = NA,
      MinLossNonMember = NA,
      Success = NA,
      Error = as.character(e)
    ))
  })
}

# =====================================================================
# 🚀 执行攻击
# =====================================================================
cat(sprintf("开始执行 %d 次MIA攻击（随机种子模式）...\n", total_runs))
cat(sprintf("全局主种子：%d\n\n", GLOBAL_SEED))

if (USE_PARALLEL) {
  cl <- makeCluster(NUM_CORES)
  
  clusterExport(cl, c(
    "expit", "generate_one_site", "calculate_loss",
    "local.glm.lasso", "REF_DS_inf", "single_attack"
  ))
  
  clusterEvalQ(cl, { library(glmnet); library(mvtnorm); library(ncvreg) })
  
  results_list <- pblapply(
    seeds, single_attack,
    n = n, p = p, p1 = p1, K = K,
    signal_strength = signal_strength,
    cl = cl
  )
  
  stopCluster(cl)
} else {
  results_list <- pblapply(
    seeds, single_attack,
    n = n, p = p, p1 = p1, K = K,
    signal_strength = signal_strength
  )
}

# 合并结果
results <- do.call(rbind, results_list)
valid_results <- results[!is.na(results$Success), ]

# =====================================================================
# 📊 统计分析
# =====================================================================
cat("\n\n========================================================\n")
cat(" MIA攻击结果统计\n")
cat("========================================================\n")

success_count <- sum(valid_results$Success)
success_rate <- mean(valid_results$Success) * 100
diffs <- valid_results$MinLossNonMember - valid_results$MinLossMember

cat(sprintf("总攻击次数：%d\n", total_runs))
cat(sprintf("成功次数：%d\n", success_count))
cat(sprintf("攻击成功率：%.2f%%\n\n", success_rate))

cat(sprintf("平均损失差：%.4f (±%.4f)\n", mean(diffs), sd(diffs)))
cat(sprintf("损失差中位数：%.4f\n", median(diffs)))
cat(sprintf("损失差范围：[%.4f, %.4f]\n\n", min(diffs), max(diffs)))

# =====================================================================
# 💾 保存结果 + 绘图
# =====================================================================
write.csv(results, "mia_attack_results.csv", row.names = FALSE)
saveRDS(results, "mia_attack_results.rds")

# 绘制损失差分布直方图 | Plot loss difference distribution histogram
p <- ggplot(valid_results, aes(x = MinLossNonMember - MinLossMember)) +
  geom_histogram(bins = 30, fill="steelblue", color="black", alpha=0.7) +
  geom_vline(xintercept = 0, color="red", linetype="dashed", linewidth=1) +
  labs(title="Fed-FDR 成员推断攻击损失差分布 | Fed-FDR MIA Loss Difference Distribution", 
       subtitle=sprintf("攻击成功率 Success rate: %.2f%%，平均损失差 Mean diff: %.4f", success_rate, mean(diffs)),
       x="非成员损失 - 成员损失 | Non-member Loss - Member Loss", y="频数 | Frequency") +
  theme_minimal()

# 直接在绘图窗口显示图片 | Display plot directly in graphics window
print(p)

cat("结果已保存：mia_attack_results.csv / mia_attack_results.rds\n")
cat("Results saved: mia_attack_results.csv / mia_attack_results.rds\n")


























# =====================================================================
# 成员推断攻击 (MIA) 批量测试 | 随机种子双语版
# Membership Inference Attack (MIA) Batch Test | Random Seed Bilingual Version
# =====================================================================
library(glmnet)
library(mvtnorm)
library(ncvreg)
library(pbapply)
library(parallel)
library(ggplot2)

# -------------------------- 🔧 参数配置区 | Parameter Configuration --------------------------
total_runs <- 10000        # 攻击次数，先100次试水，没问题改10000 | Number of attack runs, 100 for trial, change to 10000 later
GLOBAL_SEED <- 20260614  # 全局主种子，固定则结果可复现 | Global master seed, fixed for reproducible results
USE_PARALLEL <- TRUE     # 是否启用并行计算 | Whether to enable parallel computing
NUM_CORES <- detectCores() # 自动检测可用CPU核心数 | Auto-detect available CPU cores

# 模型超参数 | Model hyperparameters
n <- 200; p <- 200; p1 <- 5; K <- 5; signal_strength <- 3.0
# -------------------------------------------------------------------

# 生成随机不重复的种子池 | Generate random non-repeating seed pool
set.seed(GLOBAL_SEED)
seeds <- sample.int(1000000, total_runs) # 从100万范围内随机抽取不重复种子 | Randomly sample unique seeds from 1 to 1,000,000

# 基础工具函数 | Basic utility functions
expit <- function(x) { 1 / (1 + exp(-x)) } # Sigmoid激活函数 | Sigmoid activation function

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
  - (y * log(prob) + (1 - y) * log(1 - prob))
}

# =====================================================================
# 🎯 单次攻击函数（种子直接传入，无全局依赖）
# Single attack function (seed passed directly, no global dependency)
# =====================================================================
single_attack <- function(seed, n, p, p1, K, signal_strength) {
  
  tryCatch({
    set.seed(seed)
    
    # 生成联邦学习全局数据 | Generate federated global dataset
    Xall <- NULL; Treatall <- NULL; site <- NULL
    for (k in 1:K) {
      site_data <- generate_one_site(n, p, p1, signal_strength)
      Xall <- rbind(Xall, site_data$X)
      Treatall <- c(Treatall, site_data$y)
      site <- c(site, rep(k, n))
    }
    
    # 提取目标样本 | Extract target samples
    member_idx <- which(site == 2)[42]
    member_x <- Xall[member_idx, ]; member_y <- Treatall[member_idx]
    non_member_data <- generate_one_site(1, p, p1, signal_strength)
    non_member_x <- non_member_data$X[1, ]; non_member_y <- non_member_data$y[1]
    
    # 训练Fed-FDR模型 | Train Fed-FDR model
    local_res <- local.glm.lasso(Treatall = Treatall, Xall = Xall, site = site, weight = FALSE)
    Theta_matrix <- local_res$Theta
    
    # 计算各站点损失 | Calculate site-wise losses
    member_losses <- numeric(K)
    non_member_losses <- numeric(K)
    for (k in 1:K) {
      member_losses[k] <- calculate_loss(member_x, member_y, Theta_matrix[k, ])
      non_member_losses[k] <- calculate_loss(non_member_x, non_member_y, Theta_matrix[k, ])
    }
    
    min_loss_member <- min(member_losses)
    min_loss_non_member <- min(non_member_losses)
    
    return(data.frame(
      Seed = seed,
      MinLossMember = min_loss_member,
      MinLossNonMember = min_loss_non_member,
      Success = (min_loss_member < min_loss_non_member),
      Error = NA
    ))
  }, error = function(e) {
    return(data.frame(
      Seed = seed,
      MinLossMember = NA,
      MinLossNonMember = NA,
      Success = NA,
      Error = as.character(e)
    ))
  })
}

# =====================================================================
# 🚀 批量执行攻击 | Batch attack execution
# =====================================================================
cat(sprintf("开始执行 %d 次MIA攻击（随机种子模式）...\n", total_runs))
cat(sprintf("Starting %d MIA attacks (random seed mode)...\n", total_runs))
cat(sprintf("全局主种子：%d\n\n", GLOBAL_SEED))
cat(sprintf("Global master seed: %d\n\n", GLOBAL_SEED))

if (USE_PARALLEL) {
  cl <- makeCluster(NUM_CORES)
  
  # 导出所需函数到并行子进程 | Export required functions to parallel workers
  clusterExport(cl, c(
    "expit", "generate_one_site", "calculate_loss",
    "local.glm.lasso", "REF_DS_inf", "single_attack"
  ))
  
  # 子进程加载依赖包 | Load required packages in workers
  clusterEvalQ(cl, { library(glmnet); library(mvtnorm); library(ncvreg) })
  
  # 并行执行攻击（带进度条）| Parallel execution with progress bar
  results_list <- pblapply(
    seeds, single_attack,
    n = n, p = p, p1 = p1, K = K,
    signal_strength = signal_strength,
    cl = cl
  )
  
  stopCluster(cl)
} else {
  # 串行执行攻击（带进度条）| Serial execution with progress bar
  results_list <- pblapply(
    seeds, single_attack,
    n = n, p = p, p1 = p1, K = K,
    signal_strength = signal_strength
  )
}

# 合并所有试验结果 | Merge all trial results
results <- do.call(rbind, results_list)
valid_results <- results[!is.na(results$Success), ]

# =====================================================================
# 📊 攻击结果统计分析 | Attack Result Statistical Analysis
# =====================================================================
cat("\n\n========================================================\n")
cat(" MIA攻击结果统计\n")
cat(" MIA Attack Result Statistics\n")
cat("========================================================\n")

success_count <- sum(valid_results$Success)
success_rate <- mean(valid_results$Success) * 100
diffs <- valid_results$MinLossNonMember - valid_results$MinLossMember

cat(sprintf("总攻击次数：%d\n", total_runs))
cat(sprintf("Total attacks: %d\n", total_runs))
cat(sprintf("成功次数：%d\n", success_count))
cat(sprintf("Successful attacks: %d\n", success_count))
cat(sprintf("攻击成功率：%.2f%%\n\n", success_rate))
cat(sprintf("Attack success rate: %.2f%%\n\n", success_rate))

cat(sprintf("平均损失差：%.4f (±%.4f)\n", mean(diffs), sd(diffs)))
cat(sprintf("Mean loss difference: %.4f (±%.4f)\n", mean(diffs), sd(diffs)))
cat(sprintf("损失差中位数：%.4f\n", median(diffs)))
cat(sprintf("Median loss difference: %.4f\n", median(diffs)))
cat(sprintf("损失差范围：[%.4f, %.4f]\n\n", min(diffs), max(diffs)))
cat(sprintf("Loss difference range: [%.4f, %.4f]\n\n", min(diffs), max(diffs)))

# =====================================================================
# 💾 保存结果 + 绘制直方图 | Save Results + Plot Histogram
# =====================================================================
write.csv(results, "mia_attack_results.csv", row.names = FALSE)
saveRDS(results, "mia_attack_results.rds")

# 绘制损失差分布直方图 | Plot loss difference distribution histogram
p <- ggplot(valid_results, aes(x = MinLossNonMember - MinLossMember)) +
  geom_histogram(bins = 30, fill="steelblue", color="black", alpha=0.7) +
  geom_vline(xintercept = 0, color="red", linetype="dashed", linewidth=1) +
  labs(title="Fed-FDR 成员推断攻击损失差分布 | Fed-FDR MIA Loss Difference Distribution", 
       subtitle=sprintf("攻击成功率 Success rate: %.2f%%，平均损失差 Mean diff: %.4f", success_rate, mean(diffs)),
       x="非成员损失 - 成员损失 | Non-member Loss - Member Loss", y="频数 | Frequency") +
  theme_minimal()

# 直接在绘图窗口显示图片 | Display plot directly in graphics window
print(p)

cat("结果已保存：mia_attack_results.csv / mia_attack_results.rds\n")
cat("Results saved: mia_attack_results.csv / mia_attack_results.rds\n")







# =====================================================================
# 成员推断攻击 (MIA) 批量测试 | 随机种子双语版
# Membership Inference Attack (MIA) Batch Test | Random Seed Bilingual Version
# =====================================================================
library(glmnet)
library(mvtnorm)
library(ncvreg)
library(pbapply)
library(parallel)
library(ggplot2)

# -------------------------- 🔧 参数配置区 | Parameter Configuration --------------------------
total_runs <- 10000        # 攻击次数，先100次试水，没问题改10000 | Number of attack runs, 100 for trial, change to 10000 later
GLOBAL_SEED <- 20260614  # 全局主种子，固定则结果可复现 | Global master seed, fixed for reproducible results
USE_PARALLEL <- TRUE     # 是否启用并行计算 | Whether to enable parallel computing
NUM_CORES <- detectCores() # 自动检测可用CPU核心数 | Auto-detect available CPU cores

# 模型超参数 | Model hyperparameters
n <- 100; p <- 500; p1 <- 5; K <- 3; signal_strength <- 3.0
# -------------------------------------------------------------------

# 生成随机不重复的种子池 | Generate random non-repeating seed pool
set.seed(GLOBAL_SEED)
seeds <- sample.int(1000000, total_runs) # 从100万范围内随机抽取不重复种子 | Randomly sample unique seeds from 1 to 1,000,000

# 基础工具函数 | Basic utility functions
expit <- function(x) { 1 / (1 + exp(-x)) } # Sigmoid激活函数 | Sigmoid activation function

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
  - (y * log(prob) + (1 - y) * log(1 - prob))
}

# =====================================================================
# 🎯 单次攻击函数（种子直接传入，无全局依赖）
# Single attack function (seed passed directly, no global dependency)
# =====================================================================
single_attack <- function(seed, n, p, p1, K, signal_strength) {
  
  tryCatch({
    set.seed(seed)
    
    # 生成联邦学习全局数据 | Generate federated global dataset
    Xall <- NULL; Treatall <- NULL; site <- NULL
    for (k in 1:K) {
      site_data <- generate_one_site(n, p, p1, signal_strength)
      Xall <- rbind(Xall, site_data$X)
      Treatall <- c(Treatall, site_data$y)
      site <- c(site, rep(k, n))
    }
    
    # 提取目标样本 | Extract target samples
    member_idx <- which(site == 2)[42]
    member_x <- Xall[member_idx, ]; member_y <- Treatall[member_idx]
    non_member_data <- generate_one_site(1, p, p1, signal_strength)
    non_member_x <- non_member_data$X[1, ]; non_member_y <- non_member_data$y[1]
    
    # 训练Fed-FDR模型 | Train Fed-FDR model
    local_res <- local.glm.lasso(Treatall = Treatall, Xall = Xall, site = site, weight = FALSE)
    Theta_matrix <- local_res$Theta
    
    # 计算各站点损失 | Calculate site-wise losses
    member_losses <- numeric(K)
    non_member_losses <- numeric(K)
    for (k in 1:K) {
      member_losses[k] <- calculate_loss(member_x, member_y, Theta_matrix[k, ])
      non_member_losses[k] <- calculate_loss(non_member_x, non_member_y, Theta_matrix[k, ])
    }
    
    min_loss_member <- min(member_losses)
    min_loss_non_member <- min(non_member_losses)
    
    return(data.frame(
      Seed = seed,
      MinLossMember = min_loss_member,
      MinLossNonMember = min_loss_non_member,
      Success = (min_loss_member < min_loss_non_member),
      Error = NA
    ))
  }, error = function(e) {
    return(data.frame(
      Seed = seed,
      MinLossMember = NA,
      MinLossNonMember = NA,
      Success = NA,
      Error = as.character(e)
    ))
  })
}

# =====================================================================
# 🚀 批量执行攻击 | Batch attack execution
# =====================================================================
cat(sprintf("开始执行 %d 次MIA攻击（随机种子模式）...\n", total_runs))
cat(sprintf("Starting %d MIA attacks (random seed mode)...\n", total_runs))
cat(sprintf("全局主种子：%d\n\n", GLOBAL_SEED))
cat(sprintf("Global master seed: %d\n\n", GLOBAL_SEED))

if (USE_PARALLEL) {
  cl <- makeCluster(NUM_CORES)
  
  # 导出所需函数到并行子进程 | Export required functions to parallel workers
  clusterExport(cl, c(
    "expit", "generate_one_site", "calculate_loss",
    "local.glm.lasso", "REF_DS_inf", "single_attack"
  ))
  
  # 子进程加载依赖包 | Load required packages in workers
  clusterEvalQ(cl, { library(glmnet); library(mvtnorm); library(ncvreg) })
  
  # 并行执行攻击（带进度条）| Parallel execution with progress bar
  results_list <- pblapply(
    seeds, single_attack,
    n = n, p = p, p1 = p1, K = K,
    signal_strength = signal_strength,
    cl = cl
  )
  
  stopCluster(cl)
} else {
  # 串行执行攻击（带进度条）| Serial execution with progress bar
  results_list <- pblapply(
    seeds, single_attack,
    n = n, p = p, p1 = p1, K = K,
    signal_strength = signal_strength
  )
}

# 合并所有试验结果 | Merge all trial results
results <- do.call(rbind, results_list)
valid_results <- results[!is.na(results$Success), ]

# =====================================================================
# 📊 攻击结果统计分析 | Attack Result Statistical Analysis
# =====================================================================
cat("\n\n========================================================\n")
cat(" MIA攻击结果统计\n")
cat(" MIA Attack Result Statistics\n")
cat("========================================================\n")

success_count <- sum(valid_results$Success)
success_rate <- mean(valid_results$Success) * 100
diffs <- valid_results$MinLossNonMember - valid_results$MinLossMember

cat(sprintf("总攻击次数：%d\n", total_runs))
cat(sprintf("Total attacks: %d\n", total_runs))
cat(sprintf("成功次数：%d\n", success_count))
cat(sprintf("Successful attacks: %d\n", success_count))
cat(sprintf("攻击成功率：%.2f%%\n\n", success_rate))
cat(sprintf("Attack success rate: %.2f%%\n\n", success_rate))

cat(sprintf("平均损失差：%.4f (±%.4f)\n", mean(diffs), sd(diffs)))
cat(sprintf("Mean loss difference: %.4f (±%.4f)\n", mean(diffs), sd(diffs)))
cat(sprintf("损失差中位数：%.4f\n", median(diffs)))
cat(sprintf("Median loss difference: %.4f\n", median(diffs)))
cat(sprintf("损失差范围：[%.4f, %.4f]\n\n", min(diffs), max(diffs)))
cat(sprintf("Loss difference range: [%.4f, %.4f]\n\n", min(diffs), max(diffs)))

# =====================================================================
# 💾 保存结果 + 绘制直方图 | Save Results + Plot Histogram
# =====================================================================
write.csv(results, "mia_attack_results.csv", row.names = FALSE)
saveRDS(results, "mia_attack_results.rds")

# =====================================================================
# 💾 绘制并显示直方图 | Plot and display histogram
# =====================================================================
write.csv(results, "mia_attack_results.csv", row.names = FALSE)
saveRDS(results, "mia_attack_results.rds")

# 绘制损失差分布直方图 | Plot loss difference distribution histogram
p <- ggplot(valid_results, aes(x = MinLossNonMember - MinLossMember)) +
  geom_histogram(bins = 30, fill="steelblue", color="black", alpha=0.7) +
  geom_vline(xintercept = 0, color="red", linetype="dashed", linewidth=1) +
  labs(title="Fed-FDR 成员推断攻击损失差分布 | Fed-FDR MIA Loss Difference Distribution", 
       subtitle=sprintf("攻击成功率 Success rate: %.2f%%，平均损失差 Mean diff: %.4f", success_rate, mean(diffs)),
       x="非成员损失 - 成员损失 | Non-member Loss - Member Loss", y="频数 | Frequency") +
  theme_minimal()

# 直接在绘图窗口显示图片 | Display plot directly in graphics window
print(p)

cat("结果已保存：mia_attack_results.csv / mia_attack_results.rds\n")
cat("Results saved: mia_attack_results.csv / mia_attack_results.rds\n")







# =====================================================================
# 成员推断攻击 (MIA) 批量测试 | 随机种子双语版
# Membership Inference Attack (MIA) Batch Test | Random Seed Bilingual Version
# =====================================================================
library(glmnet)
library(mvtnorm)
library(ncvreg)
library(pbapply)
library(parallel)
library(ggplot2)

# -------------------------- 🔧 参数配置区 | Parameter Configuration --------------------------
total_runs <- 10000        # 攻击次数，先100次试水，没问题改10000 | Number of attack runs, 100 for trial, change to 10000 later
GLOBAL_SEED <- 20260614  # 全局主种子，固定则结果可复现 | Global master seed, fixed for reproducible results
USE_PARALLEL <- TRUE     # 是否启用并行计算 | Whether to enable parallel computing
NUM_CORES <- detectCores() # 自动检测可用CPU核心数 | Auto-detect available CPU cores

# 模型超参数 | Model hyperparameters
n <- 100; p <- 500; p1 <- 10; K <- 2; signal_strength <- 3.0
# -------------------------------------------------------------------

# 生成随机不重复的种子池 | Generate random non-repeating seed pool
set.seed(GLOBAL_SEED)
seeds <- sample.int(1000000, total_runs) # 从100万范围内随机抽取不重复种子 | Randomly sample unique seeds from 1 to 1,000,000

# 基础工具函数 | Basic utility functions
expit <- function(x) { 1 / (1 + exp(-x)) } # Sigmoid激活函数 | Sigmoid activation function

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
  - (y * log(prob) + (1 - y) * log(1 - prob))
}

# =====================================================================
# 🎯 单次攻击函数（种子直接传入，无全局依赖）
# Single attack function (seed passed directly, no global dependency)
# =====================================================================
single_attack <- function(seed, n, p, p1, K, signal_strength) {
  
  tryCatch({
    set.seed(seed)
    
    # 生成联邦学习全局数据 | Generate federated global dataset
    Xall <- NULL; Treatall <- NULL; site <- NULL
    for (k in 1:K) {
      site_data <- generate_one_site(n, p, p1, signal_strength)
      Xall <- rbind(Xall, site_data$X)
      Treatall <- c(Treatall, site_data$y)
      site <- c(site, rep(k, n))
    }
    
    # 提取目标样本 | Extract target samples
    member_idx <- which(site == 2)[42]
    member_x <- Xall[member_idx, ]; member_y <- Treatall[member_idx]
    non_member_data <- generate_one_site(1, p, p1, signal_strength)
    non_member_x <- non_member_data$X[1, ]; non_member_y <- non_member_data$y[1]
    
    # 训练Fed-FDR模型 | Train Fed-FDR model
    local_res <- local.glm.lasso(Treatall = Treatall, Xall = Xall, site = site, weight = FALSE)
    Theta_matrix <- local_res$Theta
    
    # 计算各站点损失 | Calculate site-wise losses
    member_losses <- numeric(K)
    non_member_losses <- numeric(K)
    for (k in 1:K) {
      member_losses[k] <- calculate_loss(member_x, member_y, Theta_matrix[k, ])
      non_member_losses[k] <- calculate_loss(non_member_x, non_member_y, Theta_matrix[k, ])
    }
    
    min_loss_member <- min(member_losses)
    min_loss_non_member <- min(non_member_losses)
    
    return(data.frame(
      Seed = seed,
      MinLossMember = min_loss_member,
      MinLossNonMember = min_loss_non_member,
      Success = (min_loss_member < min_loss_non_member),
      Error = NA
    ))
  }, error = function(e) {
    return(data.frame(
      Seed = seed,
      MinLossMember = NA,
      MinLossNonMember = NA,
      Success = NA,
      Error = as.character(e)
    ))
  })
}

# =====================================================================
# 🚀 批量执行攻击 | Batch attack execution
# =====================================================================
cat(sprintf("开始执行 %d 次MIA攻击（随机种子模式）...\n", total_runs))
cat(sprintf("Starting %d MIA attacks (random seed mode)...\n", total_runs))
cat(sprintf("全局主种子：%d\n\n", GLOBAL_SEED))
cat(sprintf("Global master seed: %d\n\n", GLOBAL_SEED))

if (USE_PARALLEL) {
  cl <- makeCluster(NUM_CORES)
  
  # 导出所需函数到并行子进程 | Export required functions to parallel workers
  clusterExport(cl, c(
    "expit", "generate_one_site", "calculate_loss",
    "local.glm.lasso", "REF_DS_inf", "single_attack"
  ))
  
  # 子进程加载依赖包 | Load required packages in workers
  clusterEvalQ(cl, { library(glmnet); library(mvtnorm); library(ncvreg) })
  
  # 并行执行攻击（带进度条）| Parallel execution with progress bar
  results_list <- pblapply(
    seeds, single_attack,
    n = n, p = p, p1 = p1, K = K,
    signal_strength = signal_strength,
    cl = cl
  )
  
  stopCluster(cl)
} else {
  # 串行执行攻击（带进度条）| Serial execution with progress bar
  results_list <- pblapply(
    seeds, single_attack,
    n = n, p = p, p1 = p1, K = K,
    signal_strength = signal_strength
  )
}

# 合并所有试验结果 | Merge all trial results
results <- do.call(rbind, results_list)
valid_results <- results[!is.na(results$Success), ]

# =====================================================================
# 📊 攻击结果统计分析 | Attack Result Statistical Analysis
# =====================================================================
cat("\n\n========================================================\n")
cat(" MIA攻击结果统计\n")
cat(" MIA Attack Result Statistics\n")
cat("========================================================\n")

success_count <- sum(valid_results$Success)
success_rate <- mean(valid_results$Success) * 100
diffs <- valid_results$MinLossNonMember - valid_results$MinLossMember

cat(sprintf("总攻击次数：%d\n", total_runs))
cat(sprintf("Total attacks: %d\n", total_runs))
cat(sprintf("成功次数：%d\n", success_count))
cat(sprintf("Successful attacks: %d\n", success_count))
cat(sprintf("攻击成功率：%.2f%%\n\n", success_rate))
cat(sprintf("Attack success rate: %.2f%%\n\n", success_rate))

cat(sprintf("平均损失差：%.4f (±%.4f)\n", mean(diffs), sd(diffs)))
cat(sprintf("Mean loss difference: %.4f (±%.4f)\n", mean(diffs), sd(diffs)))
cat(sprintf("损失差中位数：%.4f\n", median(diffs)))
cat(sprintf("Median loss difference: %.4f\n", median(diffs)))
cat(sprintf("损失差范围：[%.4f, %.4f]\n\n", min(diffs), max(diffs)))
cat(sprintf("Loss difference range: [%.4f, %.4f]\n\n", min(diffs), max(diffs)))

# =====================================================================
# 💾 保存结果 + 绘制直方图 | Save Results + Plot Histogram
# =====================================================================
write.csv(results, "mia_attack_results.csv", row.names = FALSE)
saveRDS(results, "mia_attack_results.rds")

# =====================================================================
# 💾 绘制并显示直方图 | Plot and display histogram
# =====================================================================
write.csv(results, "mia_attack_results.csv", row.names = FALSE)
saveRDS(results, "mia_attack_results.rds")

# 绘制损失差分布直方图 | Plot loss difference distribution histogram
p <- ggplot(valid_results, aes(x = MinLossNonMember - MinLossMember)) +
  geom_histogram(bins = 30, fill="steelblue", color="black", alpha=0.7) +
  geom_vline(xintercept = 0, color="red", linetype="dashed", linewidth=1) +
  labs(title="Fed-FDR 成员推断攻击损失差分布 | Fed-FDR MIA Loss Difference Distribution", 
       subtitle=sprintf("攻击成功率 Success rate: %.2f%%，平均损失差 Mean diff: %.4f", success_rate, mean(diffs)),
       x="非成员损失 - 成员损失 | Non-member Loss - Member Loss", y="频数 | Frequency") +
  theme_minimal()

# 直接在绘图窗口显示图片 | Display plot directly in graphics window
print(p)

cat("结果已保存：mia_attack_results.csv / mia_attack_results.rds\n")
cat("Results saved: mia_attack_results.csv / mia_attack_results.rds\n")