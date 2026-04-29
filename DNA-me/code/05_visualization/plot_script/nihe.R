# 安装ggplot2（若未安装）

setwd('/home/hyf/rrbs')
# 加载包
library(ggplot2)
# 读取CSV文件（假设文件名为"Age-predictions-Mouse_Blood_Clock.csv"）
data <- read.xlsx("./data/niheData.xlsx")


# ================================================================
#  Age vs sumWeighted_mCpG 拟合曲线绘制与统计分析（R语言）
#  数据文件：niheData.xlsx
# ================================================================

# -------------------- 1. 加载R包 --------------------
if (!require("ggplot2"))    install.packages("ggplot2",    repos = "https://cloud.r-project.org")
if (!require("readxl"))     install.packages("readxl",     repos = "https://cloud.r-project.org")
if (!require("ggpubr"))     install.packages("ggpubr",     repos = "https://cloud.r-project.org")
if (!require("cowplot"))    install.packages("cowplot",    repos = "https://cloud.r-project.org")

library(ggplot2)
library(readxl)
library(ggpubr)
library(cowplot)

# -------------------- 2. 读取数据 --------------------
# 从Excel读取（请确保工作目录正确，或替换为完整路径）
data <- read.xlsx("./data/niheData.xlsx")

# 备选：直接创建数据框（若无法读取Excel，取消下方注释）
# data <- data.frame(
#   sample = c("YM-WT1","YM-WT2","YF-WT1","YF-WT2",
#              "OM-WT1","OM-WT2","OF-WT1","OF-WT2",
#              "mM-WT1","mM-WT2","mF-WT1","mF-WT2"),
#   Age    = c(2.75, 2.25, 2.75, 2.75,
#              26, 26, 26, 26,
#              11, 12, 11, 12),
#   mCpG   = c(0.697747581, 0.483117253, 2.874419212, 2.35531296,
#              10.63583258, 9.254537267, 14.487772743, 15.865764763,
#              4.447047662, 4.720595366, 7.434629407, 6.804543412)
# )

# 查看数据
cat("========== 数据概览 ==========\n")
print(data)
str(data)

# -------------------- 3. 线性回归分析 --------------------
model <- lm(mCpG ~ Age, data = data)

cat("\n==================== 线性回归完整摘要 ====================\n")
print(summary(model))

# 提取关键统计量
slope     <- coef(model)[2]                        # 斜率
intercept <- coef(model)[1]                        # 截距
r_squared <- summary(model)$r.squared              # 决定系数 R²
p_value   <- summary(model)$coefficients[2, 4]     # Age系数的p值
std_err   <- summary(model)$coefficients[2, 2]     # 斜率标准误
t_value   <- summary(model)$coefficients[2, 3]     # t值
f_stat    <- summary(model)$fstatistic              # F统计量
f_pvalue  <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)

# 显著性标记
sig_mark <- ifelse(p_value < 0.001, "***",
                   ifelse(p_value < 0.01,  "**",
                          ifelse(p_value < 0.05,  "*",
                                 "ns")))

# 格式化p值显示
p_display <- ifelse(p_value < 0.001,
                    paste0("p < 0.001 ", sig_mark),
                    paste0("p = ", format(p_value, digits = 4), " ", sig_mark))

cat("\n========== 关键统计指标 ==========\n")
cat("样本量        : n =", nrow(data), "\n")
cat("斜率   :", round(slope, 4), "\n")
cat("斜率标准误 (SE):", round(std_err, 4), "\n")
cat("截距   :", round(intercept, 4), "\n")
cat("t值           :", round(t_value, 4), "\n")
cat("R-squared     :", round(r_squared, 4), "\n")
cat("p值           :", format(p_value, scientific = TRUE, digits = 3), "\n")
cat("F统计量       :", round(f_stat[1], 2), "\n")
cat("F检验p值      :", format(f_pvalue, scientific = TRUE, digits = 3), "\n")
cat("显著性        :", sig_mark, "\n")

# -------------------- 4. 模型诊断 --------------------
cat("\n========== 模型诊断 ==========\n")
# Shapiro-Wilk正态性检验
sw_test <- shapiro.test(residuals(model))
cat("Shapiro-Wilk检验 p值:", round(sw_test$p.value, 4), "\n")
if (sw_test$p.value > 0.05) {
  cat("-> 残差服从正态分布，线性模型假设成立\n")
} else {
  cat("-> 残差可能偏离正态分布，注意谨慎解读\n")
}

# Durbin-Watson自相关检验（需lmtest包）
if (!require("lmtest")) install.packages("lmtest", repos = "https://cloud.r-project.org")
library(lmtest)
dw_test <- dwtest(model)
cat("Durbin-Watson检验 p值:", round(dw_test$p.value, 4), "\n")
if (dw_test$p.value > 0.05) {
  cat("-> 残差无显著自相关\n")
} else {
  cat("-> 残差存在自相关，可能需要考虑混合效应模型\n")
}

# -------------------- 5. 绘制拟合曲线图 --------------------
# 构建标注文本（斜率 + R² + p值）
label_text <- paste0(
  "Slope = ", round(slope, 4), " \u00B1 ", round(std_err, 4), "\n",
  "R\u00B2 = ", round(r_squared, 4), "\n",
  p_display
)

# 构建回归方程文本
eq_text <- paste0(
  "y = ", round(slope, 4), "x + ", round(intercept, 4)
)

# ===== 主图：拟合曲线 =====
p_main <- ggplot(data, aes(x = Age, y = mCpG)) +
  
  # 95%置信区间（浅红色带状区域）
  geom_smooth(
    method  = "lm",
    color   = "#D32F2F",
    fill    = "#FFCDD2",
    alpha   = 0.35,
    linewidth = 1.3,
    se      = TRUE,
    fullrange = TRUE
  ) +
  
  # 散点（深蓝色，白色描边）
  geom_point(
    color   = "#1565C0",
    size    = 4,
    shape   = 16,
    stroke  = 1.2
  ) +
  
  # 每个点添加样本标签
  geom_text(
    aes(label = sample),
    vjust    = -1.8,
    hjust    = 0.5,
    size     = 3.2,
    color    = "#37474F",
    fontface = "bold"
  ) +
  
  # 标注：回归方程
  annotate(
    "text",
    x     = min(data$Age) + 2,
    y     = max(data$mCpG) - 0.3,
    label = eq_text,
    size  = 4,
    color = "#B71C1C",
    fontface = "italic",
    hjust = 0
  ) +
  
  # 标注：统计信息框
  annotate(
    "text",
    x     = min(data$Age) + 2,
    y     = max(data$mCpG) - 2.8,
    label = label_text,
    size  = 4.2,
    color = "#212121",
    fontface = "bold",
    hjust = 0,
    lineheight = 1.4
  ) +
  
  # 标题与坐标轴
  labs(
    title    = expression(paste("Age vs ", sumWeighted, "-mCpG for WT Mice")),
    subtitle = "Linear Regression with 95% Confidence Interval",
    x        = "Age (months)",
    y        = expression(paste(sumWeighted, "-mCpG")),
    caption  = paste0("n = ", nrow(data), " WT samples | ",
                      "Red line: Linear fit | Shaded area: 95% CI")
  ) +
  
  # 主题美化（黑白底，适合论文）
  theme_bw(base_size = 14) +
  theme(
    plot.title       = element_text(face = "bold", size = 17, hjust = 0.5),
    plot.subtitle    = element_text(size = 11, color = "gray50", hjust = 0.5),
    plot.caption     = element_text(size = 9, color = "gray40", hjust = 1),
    panel.grid.major = element_line(linetype = "dashed", color = "gray85"),
    panel.grid.minor = element_blank(),
    axis.title.x     = element_text(face = "bold", size = 13),
    axis.title.y     = element_text(face = "bold", size = 13),
    axis.text        = element_text(size = 11, color = "gray30")
  ) +
  
  # 坐标轴范围
  coord_cartesian(xlim = c(0, 28), ylim = c(-1, 18))


# ===== 残差诊断图（2x2布局）=====
p_resid <- ggplot(data = data, aes(x = Age, y = mCpG)) +
  geom_point(color = "#1565C0", size = 3) +
  geom_smooth(method = "lm", color = "#D32F2F", se = FALSE, linewidth = 0.8) +
  labs(x = "Age (months)", y = expression(paste(sumWeighted, "-mCpG"))) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_line(linetype = "dashed", color = "gray85"))

p_diag <- ggplot_residhist(model) +   # 残差直方图
  theme_bw(base_size = 11) +
  labs(title = "Residuals Histogram")

p_qq <- ggplot_qqplot(model) +        # Q-Q图
  stat_qq_line(color = "#D32F2F", linewidth = 0.8) +
  theme_bw(base_size = 11) +
  labs(title = "Normal Q-Q Plot")

p_resid_fit <- ggplot(model, aes(.fitted, .resid)) +
  geom_point(color = "#1565C0", size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(x = "Fitted Values", y = "Residuals", title = "Residuals vs Fitted") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_line(linetype = "dashed", color = "gray85"))

# 组合诊断图
p_combined_diag <- plot_grid(
  p_resid, p_resid_fit,
  p_qq,    p_diag,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D")
)

# -------------------- 6. 显示图形 --------------------
# 显示主图
print(p_main)

# 显示诊断图
print(p_combined_diag)
p_main
ggsave("./plot/Age_mCpG_WT_diagnostics.pdf", plot = p_main,
       width = 10, height = 8)
# -------------------- 7. 保存图形 --------------------
# 主图保存
ggsave("Age_mCpG_WT_fit_curve.png", plot = p_main,
       width = 8, height = 6, dpi = 300, bg = "white")
ggsave("Age_mCpG_WT_fit_curve.pdf", plot = p_main,
       width = 8, height = 6)

# 诊断图保存
ggsave("Age_mCpG_WT_diagnostics.png", plot = p_combined_diag,
       width = 10, height = 8, dpi = 300, bg = "white")
ggsave("Age_mCpG_WT_diagnostics.pdf", plot = p_combined_diag,
       width = 10, height = 8)

cat("\n======================================================\n")
cat("图形已保存:\n")
cat("  -> Age_mCpG_WT_fit_curve.png / .pdf  (主图)\n")
cat("  -> Age_mCpG_WT_diagnostics.png / .pdf (诊断图)\n")
cat("======================================================\n")
