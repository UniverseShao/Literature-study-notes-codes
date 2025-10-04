################################################################################
# 生成PNG图片（使用传统png设备，更兼容）
################################################################################

cat("\n========== 加载数据和包 ==========\n")

# 加载工作空间
load("Analysis_Workspace.RData")

# 加载必需的包
library(tidyverse)
library(survival)
library(survminer)
library(rms)
library(ggpubr)

cat("数据加载完成！样本量：", nrow(data), "\n\n")

# 创建图形文件夹
if (!dir.exists("Figures")) dir.create("Figures")

# =============================================================================
# 方法1：使用png()设备生成图片（更可靠）
# =============================================================================

cat("========== 生成KM曲线图片 ==========\n")

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
  title = "Figure 1A. Kaplan-Meier Curves for MACE by LV Size"
)

# 使用png设备保存
png("Figures/Figure1A_KM_MACE.png", width = 3000, height = 2400, res = 300)
print(km_plot_MACE)
dev.off()

if (file.exists("Figures/Figure1A_KM_MACE.png")) {
  cat("✓ 成功保存：Figures/Figure1A_KM_MACE.png\n")
  cat("  文件大小：", round(file.size("Figures/Figure1A_KM_MACE.png")/1024, 1), "KB\n")
} else {
  cat("✗ 保存失败\n")
}

# -----------------------------------------------------------------------------
# Figure 2: 其他结局的KM曲线（单独图）
# -----------------------------------------------------------------------------

other_outcomes <- c("CHD", "HF", "stroke", "death")
outcome_labels <- c("Coronary Heart Disease", "Heart Failure", 
                   "Ischemic Stroke", "All-Cause Mortality")

km_plots_list <- list()

for (i in 1:length(other_outcomes)) {
  outcome <- other_outcomes[i]
  label <- outcome_labels[i]
  
  cat(paste0("\n生成 Figure 2", LETTERS[i], ": ", label, "...\n"))
  
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
    title = paste0("Panel ", LETTERS[i], ". ", label)
  )
  
  km_plots_list[[i]] <- km_plot$plot
  
  # 使用png设备保存单独图形
  filename <- paste0("Figures/Figure2", LETTERS[i], "_KM_", outcome, ".png")
  png(filename, width = 2400, height = 1800, res = 300)
  print(km_plot)
  dev.off()
  
  if (file.exists(filename)) {
    cat("✓ 成功保存：", filename, "\n")
  } else {
    cat("✗ 保存失败：", filename, "\n")
  }
}

# 合并Figure 2
cat("\n生成 Figure 2 组合图...\n")
figure2_combined <- ggarrange(plotlist = km_plots_list, 
                              ncol = 2, nrow = 2,
                              common.legend = TRUE,
                              legend = "bottom")

png("Figures/Figure2_KM_All_Outcomes.png", width = 4200, height = 3600, res = 300)
print(figure2_combined)
dev.off()

if (file.exists("Figures/Figure2_KM_All_Outcomes.png")) {
  cat("✓ 成功保存：Figures/Figure2_KM_All_Outcomes.png\n")
} else {
  cat("✗ 保存失败\n")
}

# =============================================================================
# 生成RCS曲线
# =============================================================================

cat("\n========== 生成RCS曲线 ==========\n")

# 重新设置datadist
dd <- datadist(data)
options(datadist = "dd")

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

# 预测
LVEDVi_range <- seq(quantile(data$LVEDVi, 0.01, na.rm = TRUE), 
                    quantile(data$LVEDVi, 0.99, na.rm = TRUE), 
                    length.out = 100)

pred_MACE <- Predict(cox_rcs_MACE, LVEDVi = LVEDVi_range, ref.zero = TRUE)

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

# 使用png设备保存
png("Figures/Figure1B_RCS_MACE.png", width = 2400, height = 1800, res = 300)
print(rcs_plot_MACE)
dev.off()

if (file.exists("Figures/Figure1B_RCS_MACE.png")) {
  cat("✓ 成功保存：Figures/Figure1B_RCS_MACE.png\n")
} else {
  cat("✗ 保存失败\n")
}

# =============================================================================
# 生成组合图（Figure 1: KM + RCS）
# =============================================================================

cat("\n========== 生成Figure 1组合图 ==========\n")

# 简化版KM曲线（无风险表）
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
  title = "A. Kaplan-Meier Curves"
)

# 简化版RCS图形
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

# 组合
figure1_combined <- ggarrange(
  km_plot_simple$plot, 
  rcs_plot_simple, 
  ncol = 2, 
  nrow = 1,
  widths = c(1, 1)
)

png("Figures/Figure1_Combined_MACE.png", width = 4200, height = 1800, res = 300)
print(figure1_combined)
dev.off()

if (file.exists("Figures/Figure1_Combined_MACE.png")) {
  cat("✓ 成功保存：Figures/Figure1_Combined_MACE.png\n")
} else {
  cat("✗ 保存失败\n")
}

# =============================================================================
# 总结
# =============================================================================

cat("\n")
cat("=======================================================\n")
cat("       图片生成完成！                                   \n")
cat("=======================================================\n")
cat("\n生成的文件：\n")

png_files <- list.files("Figures", pattern = "\\.png$", full.names = FALSE)

if (length(png_files) > 0) {
  for (f in png_files) {
    full_path <- file.path("Figures", f)
    size_kb <- round(file.size(full_path) / 1024, 1)
    cat("  ✓", f, "-", size_kb, "KB\n")
  }
  cat("\n总共生成", length(png_files), "个PNG图片文件\n")
} else {
  cat("  ✗ 没有找到PNG文件\n")
}

cat("\n图片位置：Figures文件夹\n")
cat("图片格式：PNG (300 DPI高分辨率)\n")
cat("=======================================================\n\n")

################################################################################
# 代码结束
################################################################################

