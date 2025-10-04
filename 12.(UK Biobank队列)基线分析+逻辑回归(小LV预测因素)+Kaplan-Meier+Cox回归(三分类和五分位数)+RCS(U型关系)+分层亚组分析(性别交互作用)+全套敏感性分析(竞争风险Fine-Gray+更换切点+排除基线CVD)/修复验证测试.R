################################################################################
# 修复验证测试脚本 - 确保无NA值且图形生成成功
# 
# 日期：2025-10-04
################################################################################

cat("开始修复验证测试...\n")

# 加载必需的包
library(tidyverse)
library(survival)
library(survminer)

# 设置随机种子
set.seed(12345)

# 快速生成小样本测试数据
n_test <- 5000

# 生成测试数据
test_data <- data.frame(
  age = rnorm(n_test, mean = 55, sd = 8),
  sex = sample(c("Male", "Female"), n_test, replace = TRUE, prob = c(0.48, 0.52)),
  BMI = rnorm(n_test, mean = 27, sd = 4.5),
  hypertension = rbinom(n_test, 1, 0.30),
  diabetes = rbinom(n_test, 1, 0.08),
  prior_CVD = rbinom(n_test, 1, 0.05),
  smoking_status = sample(c("Never", "Former", "Current"), n_test, replace = TRUE,
                         prob = c(0.55, 0.30, 0.15))
)

# 生成LVEDVi（性别差异）
test_data$LVEDVi <- ifelse(
  test_data$sex == "Male",
  rnorm(n_test, mean = 78, sd = 15),
  rnorm(n_test, mean = 68, sd = 13)
)
test_data$LVEDVi <- pmax(40, pmin(test_data$LVEDVi, 150))

# 使用修复后的分组标准
test_data <- test_data %>%
  mutate(
    LV_size_category = case_when(
      sex == "Male" & LVEDVi < 68 ~ "Small",
      sex == "Male" & LVEDVi >= 68 & LVEDVi <= 88 ~ "Normal", 
      sex == "Male" & LVEDVi > 88 ~ "Large",
      sex == "Female" & LVEDVi < 60 ~ "Small",
      sex == "Female" & LVEDVi >= 60 & LVEDVi <= 76 ~ "Normal",
      sex == "Female" & LVEDVi > 76 ~ "Large"
    ),
    LV_size_category = factor(LV_size_category, levels = c("Normal", "Small", "Large"))
  )

# 生成随访时间
follow_up_years <- runif(n_test, min = 0.5, max = 8)

# 计算基础风险
base_risk <- with(test_data, {
  risk <- 0.01 + 
    0.001 * (age - 55) + 
    0.002 * (BMI > 30) +
    0.003 * hypertension +
    0.004 * diabetes +
    0.005 * prior_CVD +
    0.002 * (smoking_status == "Current")
  risk
})

# 计算LV风险效应（大幅增强）
lv_risk_effect <- with(test_data, {
  ifelse(sex == "Male",
         ifelse(LVEDVi < 68, 0.08, ifelse(LVEDVi > 88, 0.06, 0)),
         ifelse(LVEDVi < 60, 0.08, ifelse(LVEDVi > 76, 0.06, 0))
  )
})

# 生成MACE事件（修复版本 - 无NA）
mace_risk <- pmin(base_risk + lv_risk_effect * 8.0, 0.50)
# 确保概率值有效且生成的事件无NA
mace_prob <- pmax(0, pmin(mace_risk * follow_up_years / 2, 0.95))
test_data$MACE_event <- rbinom(n_test, 1, mace_prob)
test_data$MACE_time <- ifelse(test_data$MACE_event == 1,
                              runif(n_test, 0.1, follow_up_years),
                              follow_up_years)

# 检查是否有NA值
cat("\n=== 数据质量检查 ===\n")
cat("MACE_event中NA的数量：", sum(is.na(test_data$MACE_event)), "\n")
cat("MACE_time中NA的数量：", sum(is.na(test_data$MACE_time)), "\n")
cat("LV_size_category中NA的数量：", sum(is.na(test_data$LV_size_category)), "\n")

# 打印分组统计
cat("\n=== 分组统计 ===\n")
cat("LV分组分布：\n")
print(table(test_data$LV_size_category, test_data$sex))

cat("\n各组MACE事件统计：\n")
event_stats <- test_data %>%
  group_by(LV_size_category) %>%
  summarise(
    n = n(),
    events = sum(MACE_event, na.rm = TRUE),
    event_rate = round(mean(MACE_event, na.rm = TRUE) * 100, 2),
    median_followup = round(median(MACE_time, na.rm = TRUE), 2),
    .groups = 'drop'
  )
print(event_stats)

# 创建并保存KM曲线
cat("\n=== 生成修复后的KM曲线 ===\n")

if (!dir.exists("Fixed_Test_Figures")) dir.create("Fixed_Test_Figures")

# 确保数据没有缺失值
complete_data <- test_data[complete.cases(test_data[c("MACE_time", "MACE_event", "LV_size_category")]), ]
cat("完整数据样本数：", nrow(complete_data), "\n")

fit_test <- survfit(Surv(MACE_time, MACE_event) ~ LV_size_category, data = complete_data)

# 执行对数秩检验
logrank_test <- survdiff(Surv(MACE_time, MACE_event) ~ LV_size_category, data = complete_data)
cat("对数秩检验 p值：", round(1 - pchisq(logrank_test$chisq, df = length(logrank_test$n) - 1), 10), "\n")

# 绘制KM曲线
km_plot_fixed <- ggsurvplot(
  fit_test,
  data = complete_data,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  ggtheme = theme_bw(),
  palette = c("#2E9FDF", "#E74C3C", "#F39C12"),
  legend.title = "LV Size Category",
  legend.labs = c("Normal LV", "Small LV", "Large LV"),
  xlab = "Follow-up Time (years)",
  ylab = "MACE-free Survival Probability", 
  title = "Fixed Test: Kaplan-Meier Curves with No NA Values",
  font.title = 16,
  font.legend = 14,
  linewidth = 2.0,  # 使用新的linewidth参数
  xlim = c(0, 8),
  ylim = c(0, 1),
  break.time.by = 1,
  break.y.by = 0.1,
  surv.scale = "percent"
)

# 保存测试图形
ggsave("Fixed_Test_Figures/Fixed_KM_MACE.png", 
       print(km_plot_fixed$plot), 
       width = 12, height = 10, dpi = 300, bg = "white")

cat("修复后的测试图形已保存：Fixed_Test_Figures/Fixed_KM_MACE.png\n")

# 打印详细生存曲线概要
cat("\n=== 详细生存曲线概要 ===\n")
summary_times <- c(1, 2, 3, 4, 5)
surv_summary <- summary(fit_test, times = summary_times)
print(surv_summary)

# 计算每个时间点的生存率差异
cat("\n=== 各组生存率比较 ===\n")
for (time_point in summary_times) {
  temp_summary <- summary(fit_test, times = time_point)
  if (length(temp_summary$surv) == 3) {
    normal_surv <- temp_summary$surv[1]
    small_surv <- temp_summary$surv[2]  
    large_surv <- temp_summary$surv[3]
    
    cat(sprintf("%d年时点: Normal=%.1f%%, Small=%.1f%%, Large=%.1f%%\n", 
                time_point, normal_surv*100, small_surv*100, large_surv*100))
    cat(sprintf("  差异: Small vs Normal = %.1f%%, Large vs Normal = %.1f%%\n",
                (small_surv - normal_surv)*100, (large_surv - normal_surv)*100))
  }
}

cat("\n=== 测试总结 ===\n")
if (sum(is.na(test_data$MACE_event)) == 0 && sum(is.na(test_data$MACE_time)) == 0) {
  cat("✅ 数据质量检查通过：无NA值\n")
} else {
  cat("❌ 数据质量检查失败：仍有NA值\n")
}

if (file.exists("Fixed_Test_Figures/Fixed_KM_MACE.png")) {
  cat("✅ 图形生成成功\n")
} else {
  cat("❌ 图形生成失败\n")
}

cat("\n修复验证测试完成！\n")