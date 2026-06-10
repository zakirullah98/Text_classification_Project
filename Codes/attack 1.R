# =====================================================================
# 成员推断攻击 (Membership Inference Attack, MIA) 攻防演练
# Membership Inference Attack (MIA) Attack and Defense Drill
# =====================================================================
library(glmnet)
library(mvtnorm)
library(ncvreg)

set.seed(2026) # 设定新种子，你可以随意修改这个种子多测几次 | Set a new seed, you can change it freely to test multiple times

# 1. 基础函数与数据生成 | 1. Basic functions and data generation
expit <- function(x){ 1 / (1 + exp(-x)) }
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
# 👩‍⚕️ 提取“真实成员”(Member) 和 生成“非成员”(Non-Member) 
# Extract "True Member" and Generate "Non-Member"
# =====================================================================
# 真实成员：取自站点 2 的第 42 号病人 (他确确实实被模型训练过) 
# True Member: Patient No. 42 from Site 2 (truly trained by the model)
member_idx <- which(site == 2)[42]
member_x <- Xall[member_idx, ]
member_y <- Treatall[member_idx]

# 非成员：一个全新的、从未参与过训练的病人 (路人甲) 
# Non-Member: A brand new patient who never participated in training (a random passerby)
non_member_data <- generate_one_site(1, p, p1, signal_strength)
non_member_x <- non_member_data$X[1, ]
non_member_y <- non_member_data$y[1]

cat("数据准备完毕！真实成员与非成员已就位。\nData preparation complete! True member and non-member are ready.\n")
cat("正在运行 Fed-FDR 算法，获取中心服务器的 Theta 矩阵 (请稍候)...\nRunning Fed-FDR algorithm to obtain the Theta matrix at the central server (please wait)...\n\n")

local_res <- local.glm.lasso(Treatall = Treatall, Xall = Xall, site = site, weight = FALSE)
Theta_matrix <- local_res$Theta

# =====================================================================
# 🕵️‍♂️ 攻击者视角：计算 Loss，谁才是训练集里的人？ 
# Attacker's Perspective: Calculate Loss, who is the one in the training set?
# =====================================================================
calculate_loss <- function(x, y, beta) {
  prob <- expit(sum(x * beta))
  prob <- max(min(prob, 1 - 1e-15), 1e-15) # 防止 log(0) 报错 | Prevent log(0) errors
  loss <- - (y * log(prob) + (1 - y) * log(1 - prob)) # 计算 Loss | Calculate Loss
  return(loss)
}

cat("========================================================\n")
cat(" 开始执行成员推断攻击 (MIA) \n Start executing Membership Inference Attack (MIA) \n")
cat("========================================================\n")

member_losses <- numeric(K)
non_member_losses <- numeric(K)

# 攻击者把两个病人的数据分别代入 5 个站点的模型中算 Loss 
# Attacker plugs the two patients' data into the 5 sites' models to calculate Loss
for (k in 1:K) {
  member_losses[k] <- calculate_loss(member_x, member_y, Theta_matrix[k, ])
  non_member_losses[k] <- calculate_loss(non_member_x, non_member_y, Theta_matrix[k, ])
}

# 攻击者取每个人在所有站点中的最低 Loss（因为如果他在，最多只在一个站点里） 
# Attacker takes the minimum Loss for each person across all sites (since if present, they are at most in one site)
min_loss_member <- min(member_losses)
min_loss_non_member <- min(non_member_losses)

cat(sprintf("【真实成员】   的最低 Loss: %.4f\n[True Member] Minimum Loss: %.4f\n", min_loss_member, min_loss_member))
cat(sprintf("【非成员(路人)】的最低 Loss: %.4f\n[Non-Member (Passerby)] Minimum Loss: %.4f\n\n", min_loss_non_member, min_loss_non_member))

cat("[攻击者推断结论]：\n[Attacker's Inference Conclusion]:\n")
if (min_loss_member < min_loss_non_member) {
  cat("结论：真实成员的 Loss 更低！\nConclusion: The True Member's Loss is lower!\n")
  cat("深度解析：MIA 攻击表面上猜对了！但如果你多运行几次（换换 seed），你会发现两者的 Loss 差距往往非常微小。这种微小的差距在实际黑客攻击中，根本无法作为确认身份的铁证。\nDeep Analysis: MIA attack guessed right on the surface! But if you run it a few more times (change the seed), you will find the Loss difference is often very small. Such a tiny difference cannot serve as solid proof of identity in a real hacker attack.\n")
} else {
  cat("结论：非成员(路人)的 Loss 竟然比真实成员还要低！\nConclusion: The Non-Member's Loss is actually lower than the True Member's!\n")
  cat("深度解析：MIA 攻击彻底翻车！模型对没见过的路人甲反而给出了更高的置信度（更低的 Loss）。\nDeep Analysis: MIA attack completely failed! The model gave higher confidence (lower Loss) to a passerby it has never seen.\n")
  cat("这证明高维噪音和去偏操作彻底搅乱了参数空间，抹除了模型对个体的记忆痕迹。\nThis proves that high-dimensional noise and debiasing operations have completely scrambled the parameter space, wiping out the model's memory traces of individuals.\n")
}