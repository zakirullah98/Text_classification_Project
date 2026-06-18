# =======大数据版本=========（用生成的数据跑原文的方法）
# ===================== 1. 生成 5 个站点的大规模数据 =====================
set.seed(42)

# 定义 expit 函数
expit <- function(x){ 1 / (1 + exp(-x)) }

# 数据生成函数
generate_one_site <- function(n, p, p1, signal_strength) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  beta_true <- rep(0, p)
  beta_true[1:p1] <- signal_strength * sample(c(-1, 1), p1, replace = TRUE)
  mean_prob <- expit(as.vector(X %*% beta_true))
  y <- rbinom(n, 1, mean_prob)
  return(list(X = X, y = y, beta_true = beta_true))
}

n <- 200                # 单站点样本量
p <- 200                # 总特征数
p1 <- 5                 # 真实特征数
K <- 5                  # 站点数
signal_strength <- 3.0  # 信号强度

Xall <- NULL
Treatall <- NULL
site <- NULL

cat("正在生成 5 个站点的数据...\n")
for (k in 1:K) {
  site_data <- generate_one_site(n, p, p1, signal_strength)
  Xall <- rbind(Xall, site_data$X)
  Treatall <- c(Treatall, site_data$y)
  site <- c(site, rep(k, n))
}
cat(sprintf("数据生成完毕！总样本量：%d，特征数：%d，真实特征：{%s}\n\n", 
            n * K, p, paste(1:p1, collapse = ", ")))

# ===================== 2. 直接运行 Fed-FDR (不要插播旧代码) =====================
cat("正在执行本地特征筛选与去偏估计 (可能需要十几秒)...\n")
# 确保已经加载了需要的包和源文件
# library(glmnet); library(mvtnorm); library(ncvreg)
source("C:/Users/Wangzijing/Desktop/Fed-FDR-master/simulation_result/Fed_FDR_functions.R")
local_res <- local.glm.lasso(Treatall = Treatall, 
                             Xall = Xall, 
                             site = site, 
                             weight = FALSE)

Theta_matrix <- local_res$Theta

q_fdr <- 0.1
cat(sprintf("正在中心服务器执行 FDR 控制 (q = %s)...\n", q_fdr))
final_res <- MFDR(Theta = Theta_matrix, q = q_fdr)

# ===================== 3. 结果输出 =====================
selected_features <- which(final_res$FDR == 1)
true_features <- 1:p1

cat("\n=== Fed-FDR 特征选择结果 ===\n")
cat(sprintf("真实特征的索引: {%s}\n", paste(true_features, collapse = ", ")))
cat(sprintf("Fed-FDR 选出的特征: {%s}\n", paste(selected_features, collapse = ", ")))

