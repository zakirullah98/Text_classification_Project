# 改进版：结合你的优点 + 我的优点
library(ggplot2)
library(pbapply)
library(parallel)

# ---------- 🔧 集中参数配置区 ----------
total_runs <- 10000        # 攻击次数
start_seed <- 2026         # 起始随机种子
seed_step <- 1             # 种子步长
USE_PARALLEL <- TRUE       # 是否启用并行
NUM_CORES <- detectCores() # 使用的CPU核心数

# 模型参数
n <- 200; p <- 200; p1 <- 5; K <- 5; signal_strength <- 3.0

# ---------- 基础函数 ----------
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

# ---------- 单次攻击函数（所有参数显式传入，不依赖全局变量） ----------
single_attack <- function(i, start_seed, seed_step, n, p, p1, K, signal_strength) {
  seed_i <- start_seed + (i-1) * seed_step
  set.seed(seed_i)
  
  tryCatch({
    # 生成 K 个站点数据
    Xall <- NULL; Treatall <- NULL; site <- NULL
    for (k in 1:K) {
      site_data <- generate_one_site(n, p, p1, signal_strength)
      Xall <- rbind(Xall, site_data$X)
      Treatall <- c(Treatall, site_data$y)
      site <- c(site, rep(k, n))
    }
    
    # 提取成员与非成员样本
    member_idx <- which(site == 2)[42]
    member_x <- Xall[member_idx, ]; member_y <- Treatall[member_idx]
    non_member_data <- generate_one_site(1, p, p1, signal_strength)
    non_member_x <- non_member_data$X[1, ]; non_member_y <- non_member_data$y[1]
    
    # 训练 Fed-FDR 模型
    local_res <- local.glm.lasso(Treatall = Treatall, Xall = Xall, site = site, weight = FALSE)
    Theta_matrix <- local_res$Theta
    
    # 计算各站点损失
    member_losses <- numeric(K)
    non_member_losses <- numeric(K)
    for (k in 1:K) {
      member_losses[k] <- calculate_loss(member_x, member_y, Theta_matrix[k, ])
      non_member_losses[k] <- calculate_loss(non_member_x, non_member_y, Theta_matrix[k, ])
    }
    
    min_loss_member <- min(member_losses)
    min_loss_non_member <- min(non_member_losses)
    
    return(data.frame(
      Seed = seed_i,
      MinLossMember = min_loss_member,
      MinLossNonMember = min_loss_non_member,
      Success = (min_loss_member < min_loss_non_member),
      Error = NA
    ))
  }, error = function(e) {
    # 出错时记录错误，不中断整个循环
    return(data.frame(
      Seed = seed_i,
      MinLossMember = NA,
      MinLossNonMember = NA,
      Success = NA,
      Error = as.character(e)
    ))
  })
}


# ---------- 执行攻击 ----------
cat(sprintf("开始执行 %d 次攻击...\n", total_runs))

if (USE_PARALLEL) {
  cl <- makeCluster(NUM_CORES)
  
  # ✅ 这里只需要导出函数，不需要导出任何参数变量了
  clusterExport(cl, c(
    "expit", "generate_one_site", "calculate_loss",
    "local.glm.lasso", "REF_DS_inf", "single_attack"
  ))
  
  clusterEvalQ(cl, { library(glmnet); library(mvtnorm); library(ncvreg) })
  
  # ✅ 把所有参数通过 pblapply 直接传给函数
  results_list <- pblapply(
    1:total_runs, single_attack,
    start_seed = start_seed,
    seed_step = seed_step,
    n = n, p = p, p1 = p1, K = K,
    signal_strength = signal_strength,
    cl = cl
  )
  
  stopCluster(cl)
} else {
  # 串行模式也同样显式传参
  results_list <- pblapply(
    1:total_runs, single_attack,
    start_seed = start_seed,
    seed_step = seed_step,
    n = n, p = p, p1 = p1, K = K,
    signal_strength = signal_strength
  )
}

# 合并所有结果
results <- do.call(rbind, results_list)

# ---------- 保存结果 ----------
write.csv(results, "attack_results_full.csv", row.names = FALSE)
saveRDS(results, "attack_results_full.rds")

# ---------- 统计分析 ----------
valid_results <- results[!is.na(results$Success), ]
success_count <- sum(valid_results$Success)
success_rate <- mean(valid_results$Success) * 100
diffs <- valid_results$MinLossNonMember - valid_results$MinLossMember

cat("\n\n========================================================\n")
cat(" 最终统计结果\n")
cat("========================================================\n")
cat(sprintf("总攻击次数：%d\n", total_runs))
cat(sprintf("成功次数：%d\n", success_count))
cat(sprintf("攻击成功率：%.2f%%\n\n", success_rate))

cat(sprintf("平均损失差：%.4f (±%.4f)\n", mean(diffs), sd(diffs)))
cat(sprintf("损失差中位数：%.4f\n", median(diffs)))
cat(sprintf("损失差范围：[%.4f, %.4f]\n\n", min(diffs), max(diffs)))

# ---------- 绘制直方图 ----------
p <- ggplot(valid_results, aes(x = MinLossNonMember - MinLossMember)) +
  geom_histogram(bins = 30, fill="steelblue", color="black", alpha=0.7) +
  geom_vline(xintercept = 0, color="red", linetype="dashed", linewidth=1) +
  labs(title="Fed-FDR 成员推断攻击损失差分布", 
       subtitle=sprintf("攻击成功率：%.2f%%，平均损失差：%.4f", success_rate, mean(diffs)),
       x="非成员损失 - 成员损失", y="频数") +
  theme_minimal()

ggsave("loss_difference_histogram.png", p, width=8, height=6, dpi=300)
cat("直方图已保存为 loss_difference_histogram.png\n")
# 列出当前工作目录下的所有文件
list.files()

# 显示文件所在的完整文件夹路径
getwd()
