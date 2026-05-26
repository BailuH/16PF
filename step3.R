# 导入因子分析相关的R包
library(psych)
library(polycor)
library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================
# 1. 读取数据 & 提取人格题项
# ============================================================

# 根据你的 Python 导出格式调整 sep：逗号用 ",", 制表符用 "\t"
df <- read.csv("data_cleaned.csv", sep = "\t", check.names = FALSE)

# 排除人口统计学变量（根据你的实际列名调整！）
exclude_cols <- c("age", "gender", "accuracy", "country", "source", "elapsed")
personality_items <- df[, !(names(df) %in% exclude_cols)]

# 确保全为数值型
personality_items <- as.data.frame(lapply(personality_items, as.numeric))

n_obs  <- nrow(personality_items)
n_vars <- ncol(personality_items)

cat("数据维度:", n_obs, "行 x", n_vars, "列\n")
cat("题项示例:", head(names(personality_items), 5), "...\n")

# ============================================================
# 2. 相关矩阵计算
# ============================================================

# 李克特 5 点量表理论上应该用 Polychoric，但 n=50000, p=163 时计算较慢。
# 如果内存不足或等待时间过长，可将下面设为 FALSE，改用 Pearson。
use_polychoric <- TRUE

if (use_polychoric) {
  cat("\n正在计算 Polychoric 相关矩阵，请耐心等待（大矩阵可能需要 2-5 分钟）...\n")
  # ML=FALSE 和 std.err=FALSE 可显著加速
  cor_obj <- hetcor(personality_items, ML = FALSE, std.err = FALSE)
  cor_matrix <- cor_obj$correlations
  saveRDS(cor_matrix, "polychoric_cor_matrix.rds")
  cat("Polychoric 相关矩阵已保存到 polychoric_cor_matrix.rds\n")
} else {
  cat("\n使用 Pearson 相关矩阵（快速模式）\n")
  cor_matrix <- cor(personality_items, use = "pairwise.complete.obs")
}

# ============================================================
# 3. 碎石图 (Scree Plot) + Kaiser-Guttman 准则
# ============================================================

cat("\n===== 3. 碎石图 & Kaiser-Guttman =====\n")

ev <- eigen(cor_matrix)$values
n_kaiser <- sum(ev > 1)

cat("特征值 > 1 的数量 (Kaiser-Guttman):", n_kaiser, "\n")

png("03_scree_plot.png", width = 900, height = 600, res = 120)
scree(cor_matrix, factors = TRUE, pc = TRUE,
      main = "Scree Plot with Kaiser-Guttman Criterion")
abline(h = 1, col = "red", lty = 2, lwd = 2)
text(x = n_vars * 0.7, y = 1.3, labels = "Eigenvalue = 1", col = "red")
dev.off()
cat("碎石图已保存: 03_scree_plot.png\n")

# ============================================================
# 4. 平行分析 (Parallel Analysis) —— 核心方法
# ============================================================

cat("\n===== 4. 平行分析 (Parallel Analysis) =====\n")

# 建议：测试时用 n.iter = 100，正式报告用 n.iter = 1000
n_iter <- 100
cat("迭代次数:", n_iter, "（建议正式分析改为 1000）\n")

cat("\n[4.1] PCA 平行分析...\n")
pa_pc <- fa.parallel(cor_matrix, n.obs = n_obs, fa = "pc",
                     n.iter = n_iter, plot = FALSE)
cat("  PCA 平行分析建议因子数:", pa_pc$ncomp, "\n")

cat("\n[4.2] EFA 平行分析...\n")
pa_fa <- fa.parallel(cor_matrix, n.obs = n_obs, fa = "fa",
                     n.iter = n_iter, plot = FALSE)
cat("  EFA 平行分析建议因子数:", pa_fa$nfact, "\n")

# 保存平行分析图
png("03_parallel_analysis.png", width = 1000, height = 700, res = 120)
fa.parallel(cor_matrix, n.obs = n_obs, fa = "fa",
            n.iter = n_iter, main = "Parallel Analysis (EFA)")
dev.off()
cat("平行分析图已保存: 03_parallel_analysis.png\n")

# ============================================================
# 5. VSS (Very Simple Structure) + MAP
# ============================================================

cat("\n===== 5. VSS & MAP =====\n")

vss_max <- min(20, n_vars - 1)
cat("测试因子数范围: 1 -", vss_max, "\n")

vss_result <- VSS(cor_matrix, n = vss_max, n.obs = n_obs, rotate = "promax")

# VSS 复杂度最低 —— 安全检查
vss_complexity <- vss_result$vss.stats$complexity
n_vss <- which.min(vss_complexity)

# 修复：which.min 可能返回 integer(0)
if (length(n_vss) == 0 || is.infinite(min(vss_complexity, na.rm = TRUE))) {
  n_vss <- NA_integer_
  cat("VSS 复杂度全部为 Inf/NA，无法确定最优因子数\n")
} else {
  cat("VSS 复杂度最低时的因子数:", n_vss, 
      "(complexity =", min(vss_complexity, na.rm = TRUE), ")\n")
}

# MAP (Minimum Average Partial)
map_vals <- vss_result$map
valid_map <- which(!is.na(map_vals))
n_map <- valid_map[which.min(map_vals[valid_map])]

# 修复：同样检查 MAP 结果
if (length(n_map) == 0) n_map <- NA_integer_

cat("MAP 最小值对应的因子数:", n_map, 
    "(MAP =", min(map_vals, na.rm = TRUE), ")\n")

# 保存 VSS 图
png("03_vss_plot.png", width = 800, height = 600, res = 120)
VSS.plot(vss_result, title = "VSS Complexity by Factor Number")
dev.off()
cat("VSS 图已保存: 03_vss_plot.png\n")

# ============================================================
# 6. Hull Method（基于拟合指数 vs 简约性）
# ============================================================

cat("\n===== 6. Hull Method =====\n")
cat("正在拟合 1-25 因子模型以计算 RMSEA / BIC / CFI（约需 1-3 分钟）...\n")

hull_n <- 1:25
hull_results <- data.frame(
  nfac  = hull_n,
  df    = NA_real_,
  chisq = NA_real_,
  rmsea = NA_real_,
  bic   = NA_real_,
  cfi   = NA_real_
)

for (i in seq_along(hull_n)) {
  k <- hull_n[i]
  fit <- tryCatch(
    fa(cor_matrix, nfactors = k, n.obs = n_obs, fm = "ml",
       rotate = "oblimin", max.iter = 100, warnings = FALSE),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  if (!is.null(fit) && !is.null(fit$RMSEA)) {
    hull_results$df[i]    <- fit$dof
    hull_results$chisq[i] <- fit$STATISTIC
    hull_results$rmsea[i] <- fit$RMSEA[1]
    hull_results$bic[i]   <- fit$BIC
    hull_results$cfi[i]   <- ifelse(is.null(fit$CFI), NA, fit$CFI)
  }
}

hull_results <- na.omit(hull_results)

# 绘制 Hull 图（以 RMSEA 为例）
hull_plot <- ggplot(hull_results, aes(x = nfac, y = rmsea)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "red", size = 3) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "darkgreen") +
  annotate("text", x = max(hull_results$nfac) * 0.8, y = 0.055,
           label = "RMSEA = 0.05", color = "darkgreen") +
  labs(title = "Hull Method: RMSEA vs. Number of Factors",
       x = "Number of Factors", y = "RMSEA") +
  theme_minimal(base_size = 14)

ggsave("03_hull_rmsea.png", hull_plot, width = 8, height = 6, dpi = 150)
cat("Hull 图已保存: 03_hull_rmsea.png\n")

# 找 RMSEA < 0.05 的最小因子数
rmsea_candidates <- hull_results$nfac[hull_results$rmsea < 0.05]
n_rmsea <- ifelse(length(rmsea_candidates) > 0, min(rmsea_candidates), NA)
cat("Hull/RMSEA < 0.05 的最小因子数:", n_rmsea, "\n")

# ============================================================
# 7. 样本量敏感性分析
# ============================================================

cat("\n===== 7. 样本量敏感性分析 =====\n")
cat("从总体中抽取 25%, 50%, 75%, 100% 子样本，重复平行分析...\n")

set.seed(42)
proportions <- c(0.25, 0.50, 0.75, 1.0)
sensitivity <- data.frame(
  proportion   = proportions,
  n_sample     = NA_integer_,
  pa_suggested = NA_integer_
)

for (i in seq_along(proportions)) {
  prop <- proportions[i]
  n_samp <- round(n_obs * prop)
  sensitivity$n_sample[i] <- n_samp

  idx <- sample(1:n_obs, n_samp)
  sub_data <- personality_items[idx, ]

  # 子样本用 Pearson 相关加速
  sub_cor <- cor(sub_data, use = "pairwise.complete.obs")

  # 子样本迭代次数设为 50 加速（正式可改 100-200）
  pa_sub <- fa.parallel(sub_cor, n.obs = n_samp, fa = "fa",
                        n.iter = 50, plot = FALSE)
  sensitivity$pa_suggested[i] <- pa_sub$nfact

  cat("  样本比例", sprintf("%.0f%%", prop * 100),
      "(n =", n_samp, ") -> 建议因子数:", pa_sub$nfact, "\n")
}

# ============================================================
# 8. 结果汇总与决策表
# ============================================================

cat("\n===== 8. 因子保留决策汇总 =====\n")
cat("==============================================\n")

# 确保所有值都是标量，防止 integer(0) 或 NULL
safe_scalar <- function(x, default = NA) {
  if (is.null(x) || length(x) == 0) return(default)
  return(x[1])
}

summary_df <- data.frame(
  方法 = c(
    "Kaiser-Guttman (特征值>1)",
    "碎石图 (目视拐点)",
    "平行分析 (PCA)",
    "平行分析 (EFA)",
    "VSS (复杂度最低)",
    "MAP (最小平均偏相关)",
    "Hull/RMSEA (<0.05)",
    "样本量敏感性 (100%子样本)"
  ),
  建议因子数 = c(
    safe_scalar(n_kaiser, NA),
    "目视判断（见 03_scree_plot.png）",
    safe_scalar(pa_pc$ncomp, NA),
    safe_scalar(pa_fa$nfact, NA),
    safe_scalar(n_vss, NA),
    safe_scalar(n_map, NA),
    safe_scalar(ifelse(is.na(n_rmsea), "无", n_rmsea), "无"),
    safe_scalar(sensitivity$pa_suggested[4], NA)
  ),
  stringsAsFactors = FALSE
)

print(summary_df, row.names = FALSE, right = FALSE)

cat("\n==============================================\n")
cat("样本量敏感性分析结果:\n")
print(sensitivity, row.names = FALSE)

cat("\n==============================================\n")
cat("【最终决策提示】\n")
cat("1. 若多种方法一致建议某因子数（如 16），则决策明确。\n")
cat("2. 若平行分析(EFA)与 VSS/MAP 冲突，优先信任平行分析。\n")
cat("3. 请结合 Cattell 16因子理论假设与 Big Five 替代模型综合判断。\n")

# ============================================================
# 9. 保存关键对象供步骤 4 使用
# ============================================================

saveRDS(cor_matrix,      "step3_cor_matrix.rds")
saveRDS(vss_result,      "step3_vss_result.rds")
saveRDS(hull_results,    "step3_hull_results.rds")
saveRDS(pa_fa,           "step3_pa_fa.rds")
saveRDS(summary_df,      "step3_summary.rds")

cat("\n步骤 3 完成！所有中间结果已保存为 step3_*.rds，可直接载入步骤 4 使用。\n")
cat("生成的图表:\n")
cat("  - 03_scree_plot.png\n")
cat("  - 03_parallel_analysis.png\n")
cat("  - 03_vss_plot.png\n")
cat("  - 03_hull_rmsea.png\n")