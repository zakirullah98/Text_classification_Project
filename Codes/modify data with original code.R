# =======大数据版本=========（用生成的数据跑原文的方法）
# ======= Big Data Version ========= (Running the original method with generated data)
# ===================== 1. 生成 5 个站点的大规模数据 | 1. Generate large-scale data for 5 sites =====================
set.seed(42)

# 定义 expit 函数 | Define the expit function
expit <- function(x){ 1 / (1 + exp(-x)) }

# 数据生成函数 | Data generation function
generate_one_site <- function(n, p, p1, signal_strength) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  beta_true <- rep(0, p)
  beta_true[1:p1] <- signal_strength * sample(c(-1, 1), p1, replace = TRUE)
  mean_prob <- expit(as.vector(X %*% beta_true))
  y <- rbinom(n, 1, mean_prob)
  return(list(X = X, y = y, beta_true = beta_true))
}

n <- 200                # 单站点样本量 | Sample size per site
p <- 200                # 总特征数 | Total number of features
p1 <- 5                 # 真实特征数 | Number of true features
K <- 5                  # 站点数 | Number of sites
signal_strength <- 3.0  # 信号强度 | Signal strength

Xall <- NULL
Treatall <- NULL
site <- NULL

cat("正在生成 5 个站点的数据... | Generating data for 5 sites...\n")
for (k in 1:K) {
  site_data <- generate_one_site(n, p, p1, signal_strength)
  Xall <- rbind(Xall, site_data$X)
  Treatall <- c(Treatall, site_data$y)
  site <- c(site, rep(k, n))
}

cat(sprintf("数据生成完毕！总样本量：%d，特征数：%d，真实特征：{%s}\nData generation completed! Total sample size: %d, Features: %d, True features: {%s}\n\n", 
            n * K, p, paste(1:p1, collapse = ", "), 
            n * K, p, paste(1:p1, collapse = ", ")))

# ===================== 2. 直接运行 Fed-FDR (不要插播旧代码) | 2. Run Fed-FDR directly (Do not insert old code) =====================
cat("正在执行本地特征筛选与去偏估计 (可能需要十几秒)... \nExecuting local feature screening and debiased estimation (may take a few seconds)...\n")
# 确保已经加载了需要的包和源文件 | Ensure required packages and source files are loaded
# library(glmnet); library(mvtnorm); library(ncvreg)
source("C:/Users/Wangzijing/Desktop/Fed-FDR-master/simulation_result/Fed_FDR_functions.R")
local_res <- local.glm.lasso(Treatall = Treatall, 
                             Xall = Xall, 
                             site = site, 
                             weight = FALSE)

Theta_matrix <- local_res$Theta

q_fdr <- 0.1
cat(sprintf("正在中心服务器执行 FDR 控制 (q = %s)... \nExecuting FDR control on the central server (q = %s)...\n", q_fdr, q_fdr))
final_res <- MFDR(Theta = Theta_matrix, q = q_fdr)

# ===================== 3. 结果输出 | 3. Output Results =====================
selected_features <- which(final_res$FDR == 1)
true_features <- 1:p1

cat("\n=== Fed-FDR 特征选择结果 | Fed-FDR Feature Selection Results ===\n")
cat(sprintf("真实特征的索引 | Indices of true features: {%s}\n", paste(true_features, collapse = ", ")))
cat(sprintf("Fed-FDR 选出的特征 | Features selected by Fed-FDR: {%s}\n", paste(selected_features, collapse = ", ")))