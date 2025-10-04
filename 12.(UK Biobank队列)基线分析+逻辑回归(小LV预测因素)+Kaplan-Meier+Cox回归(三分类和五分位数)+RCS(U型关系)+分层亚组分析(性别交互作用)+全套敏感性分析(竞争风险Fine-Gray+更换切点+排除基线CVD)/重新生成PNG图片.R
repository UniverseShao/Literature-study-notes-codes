################################################################################
# 重新生成PNG格式图片（解决PDF显示问题）
################################################################################

# 设置工作目录
setwd("f:/桌面/机器学习和预测模型和高分文献/12.(UK Biobank队列)基线分析+逻辑回归(小LV预测因素)+Kaplan-Meier+Cox回归(三分类和五分位数)+RCS(U型关系)+分层亚组分析(性别交互作用)+全套敏感性分析(竞争风险Fine-Gray+更换切点+排除基线CVD)")

cat("\n当前工作目录：", getwd(), "\n\n")

cat("========== 加载之前的分析结果 ==========\n")

# 加载之前保存的工作空间
load("Analysis_Workspace.RData")

# 加载必需的包
library(tidyverse)
library(survival)
library(survminer)
library(rms)
library(ggpubr)

cat("数据加载完成！\n")
cat("样本量：", nrow(data), "\n\n")

# 创建PNG文件夹（使用完整路径）
png_dir <- "Figures_PNG"
if (!dir.exists(png_dir)) {
  dir.create(png_dir)
  cat("已创建文件夹：", png_dir, "\n")
}

cat("PNG文件保存路径：", file.path(getwd(), png_dir), "\n\n")

# =============================================================================
# 重新生成KM曲线（PNG格式）
# =============================================================================

cat("========== 生成Kaplan-Meier曲线（PNG格式）==========\n\n")

# -----------------------------------------------------------------------------
# Figure 1A: MACE的KM曲线
# -----------------------------------------------------------------------------

cat("生成 Figure 1A: MACE KM曲线...\n")

fit_MACE <- survfit(Surv(MACE_time, MACE_event) ~ LV_size_category, data = data)

km_plot_MACE <- ggsurvplot(
  fit_MACE,
  data = data,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  ggtheme = theme_bw(),
  palette = c("#4DBBD5", "#E64B35", "#00A087"),
  legend.title = "LV Size",
  legend.labs = c("Normal", "Small", "Large"),
  xlab = "Time (years)",
  ylab = "Event-free Survival",
  title = "Figure 1A. Kaplan-Meier Curves for MACE by LV Size",
  font.main = c(14, "bold"),
  font.x = c(12, "plain"),
  font.y = c(12, "plain"),
  font.legend = c(11, "plain")
)

# 保存为PNG（高分辨率）
png_file <- file.path(png_dir, "Figure1A_KM_MACE.png")
ggsave(png_file, 
       print(km_plot_MACE$plot), 
       width = 10, height = 8, dpi = 300, bg = "white")

if (file.exists(png_file)) {
  cat("✓ 已保存：", png_file, "\n")
  cat("  文件大小：", file.size(png_file), "字节\n")
} else {
  cat("✗ 保存失败：", png_file, "\n")
}

# -----------------------------------------------------------------------------
# Figure 2: 其他结局的KM曲线
# -----------------------------------------------------------------------------

other_outcomes <- c("CHD", "HF", "stroke", "death")
outcome_labels <- c("Coronary Heart Disease", "Heart Failure", 
                   "Ischemic Stroke", "All-Cause Mortality")

km_plots_list <- list()

for (i in 1:length(other_outcomes)) {
  outcome <- other_outcomes[i]
  label <- outcome_labels[i]
  
  cat(paste0("\n生成 Figure 2", LETTERS[i], ": ", label, " KM曲线...\n"))
  
  event_var <- paste0(outcome, "_event")
  time_var <- paste0(outcome, "_time")
  
  fit <- survfit(as.formula(paste0("Surv(", time_var, ", ", event_var, 
                                   ") ~ LV_size_category")), 
                data = data)
  
  km_plot <- ggsurvplot(
    fit,
    data = data,
    pval = TRUE,
    conf.int = FALSE,
    risk.table = FALSE,
    ggtheme = theme_bw(),
    palette = c("#4DBBD5", "#E64B35", "#00A087"),
    legend.title = "LV Size",
    legend.labs = c("Normal", "Small", "Large"),
    xlab = "Time (years)",
    ylab = "Event-free Survival",
    title = paste0("Panel ", LETTERS[i], ". ", label),
    font.main = c(12, "bold"),
    font.x = c(11, "plain"),
    font.y = c(11, "plain"),
    font.legend = c(10, "plain")
  )
  
  km_plots_list[[i]] <- km_plot$plot
  
  # 保存单独的PNG图形
  png_file <- file.path(png_dir, paste0("Figure2", LETTERS[i], "_KM_", outcome, ".png"))
  ggsave(png_file, km_plot$plot, width = 8, height = 6, dpi = 300, bg = "white")
  
  if (file.exists(png_file)) {
    cat("✓ 已保存：", png_file, "\n")
  } else {
    cat("✗ 保存失败：", png_file, "\n")
  }
}

# 合并四个图形为Figure 2
cat("\n生成 Figure 2 组合图...\n")
figure2_combined <- ggarrange(plotlist = km_plots_list, 
                              ncol = 2, nrow = 2,
                              common.legend = TRUE,
                              legend = "bottom")

png_file <- file.path(png_dir, "Figure2_KM_All_Outcomes.png")
ggsave(png_file, figure2_combined, width = 14, height = 12, dpi = 300, bg = "white")

if (file.exists(png_file)) {
  cat("✓ 已保存：", png_file, "\n")
} else {
  cat("✗ 保存失败：", png_file, "\n")
}

# =============================================================================
# 重新生成RCS曲线（PNG格式）
# =============================================================================

cat("\n========== 生成RCS曲线（PNG格式）==========\n\n")

# 重新设置datadist
dd <- datadist(data)
options(datadist = "dd")

# -----------------------------------------------------------------------------
# Figure 1B: MACE的RCS曲线
# -----------------------------------------------------------------------------

cat("生成 Figure 1B: MACE RCS曲线...\n")

# 构建RCS模型
cox_rcs_MACE <- cph(
  Surv(MACE_time, MACE_event) ~ rcs(LVEDVi, 4) + age + sex + ethnicity +
    education + smoking_status + alcohol_status + BMI + SBP + DBP +
    total_cholesterol + hypertension + diabetes + HbA1c + CRP +
    hemoglobin + albumin + atrial_fibrillation + statin_use +
    antihypertensive_use + LVEF,
  data = data,
  x = TRUE,
  y = TRUE
)

# 预测HR和置信区间
LVEDVi_range <- seq(quantile(data$LVEDVi, 0.01, na.rm = TRUE), 
                    quantile(data$LVEDVi, 0.99, na.rm = TRUE), 
                    length.out = 100)

pred_MACE <- Predict(cox_rcs_MACE, LVEDVi = LVEDVi_range, ref.zero = TRUE)

# 转换为HR
pred_MACE_df <- data.frame(
  LVEDVi = pred_MACE$LVEDVi,
  HR = exp(pred_MACE$yhat),
  CI_lower = exp(pred_MACE$lower),
  CI_upper = exp(pred_MACE$upper)
)

# 绘制RCS图形
rcs_plot_MACE <- ggplot(pred_MACE_df, aes(x = LVEDVi, y = HR)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), 
              fill = "lightblue", alpha = 0.3) +
  geom_line(color = "blue", linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    x = "LVEDVi (mL/m²)",
    y = "Hazard Ratio for MACE",
    title = "Figure 1B. Adjusted Association Between LVEDVi and MACE (RCS)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

# 保存为PNG
png_file <- file.path(png_dir, "Figure1B_RCS_MACE.png")
ggsave(png_file, rcs_plot_MACE, width = 8, height = 6, dpi = 300, bg = "white")

if (file.exists(png_file)) {
  cat("✓ 已保存：", png_file, "\n")
} else {
  cat("✗ 保存失败：", png_file, "\n")
}

# =============================================================================
# 生成组合图片（Figure 1: KM + RCS）
# =============================================================================

cat("\n========== 生成组合图片 ==========\n\n")

# 重新加载Figure 1A的图形（不含风险表）
fit_MACE <- survfit(Surv(MACE_time, MACE_event) ~ LV_size_category, data = data)

km_plot_simple <- ggsurvplot(
  fit_MACE,
  data = data,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = FALSE,
  ggtheme = theme_bw(),
  palette = c("#4DBBD5", "#E64B35", "#00A087"),
  legend.title = "LV Size",
  legend.labs = c("Normal", "Small", "Large"),
  xlab = "Time (years)",
  ylab = "Event-free Survival",
  title = "A. Kaplan-Meier Curves",
  font.main = c(12, "bold")
)

# 修改RCS图形标题
rcs_plot_simple <- ggplot(pred_MACE_df, aes(x = LVEDVi, y = HR)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), 
              fill = "lightblue", alpha = 0.3) +
  geom_line(color = "blue", linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    x = "LVEDVi (mL/m²)",
    y = "Hazard Ratio",
    title = "B. Dose-Response Curve (RCS)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

# 组合Figure 1
cat("生成 Figure 1 组合图...\n")
figure1_combined <- ggarrange(
  km_plot_simple$plot, 
  rcs_plot_simple, 
  ncol = 2, 
  nrow = 1,
  widths = c(1, 1),
  labels = c("", ""),
  common.legend = FALSE
)

png_file <- file.path(png_dir, "Figure1_Combined_MACE.png")
ggsave(png_file, figure1_combined, width = 14, height = 6, dpi = 300, bg = "white")

if (file.exists(png_file)) {
  cat("✓ 已保存：", png_file, "\n")
} else {
  cat("✗ 保存失败：", png_file, "\n")
}

# =============================================================================
# 列出所有生成的文件
# =============================================================================

cat("\n")
cat("=======================================================\n")
cat("       PNG格式图片生成完成！                           \n")
cat("=======================================================\n")
cat("\n检查生成的文件：\n\n")

png_files <- list.files(png_dir, pattern = "\\.png$", full.names = TRUE)

if (length(png_files) > 0) {
  for (f in png_files) {
    size_kb <- round(file.size(f) / 1024, 1)
    cat("  ✓", basename(f), "-", size_kb, "KB\n")
  }
  cat("\n总共生成", length(png_files), "个PNG文件\n")
} else {
  cat("  ✗ 没有找到PNG文件！\n")
}

cat("\n文件保存位置：", file.path(getwd(), png_dir), "\n")
cat("=======================================================\n\n")

################################################################################
# 代码结束
################################################################################
