################################################################################
# 极端增强版测试脚本 - 确保KM曲线明显分离
# 
# 日期：2025-10-04
################################################################################

cat("开始极端增强版测试...\n")

# 加载必需的包
library(tidyverse)
library(survival)
library(survminer)

# 设置随机种子
set.seed(12345)

# 快速生成小样本测试数据
n_test <- 3000

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

# 使用极端增强的分组标准
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
follow_up_years <- runif(n_test, min = 0.5, max = 8)  # 缩短随访时间

# 计算基础风险
base_risk <- with(test_data, {
  risk <- 0.02 +  # 增加基础风险
    0.002 * (age - 55) + 
    0.004 * (BMI > 30) +
    0.006 * hypertension +
    0.008 * diabetes +
    0.010 * prior_CVD +
    0.004 * (smoking_status == "Current")
  risk
})

# 计算LV风险效应（极端增强）
lv_risk_effect <- with(test_data, {
  ifelse(sex == "Male",
         ifelse(LVEDVi < 68, 0.15, ifelse(LVEDVi > 88, 0.10, 0)),  # 极大增强
         ifelse(LVEDVi < 60, 0.15, ifelse(LVEDVi > 76, 0.10, 0))
  )
})

# 生成MACE事件（极端增强事件率和组间差异）
mace_risk <- pmin(base_risk + lv_risk_effect * 2.0, 0.80)  # 允许更高的事件率
test_data$MACE_event <- rbinom(n_test, 1, mace_risk * follow_up_years / 1.5)
test_data$MACE_time <- ifelse(test_data$MACE_event == 1,
                              runif(n_test, 0.1, follow_up_years),
                              follow_up_years)

# 打印分组统计
cat("\n=== 极端增强测试数据分组统计 ===\n")
cat("LV分组分布：\n")
print(table(test_data$LV_size_category, test_data$sex))

cat("\n各组MACE事件统计：\n")
event_stats <- test_data %>%
  group_by(LV_size_category) %>%
  summarise(
    n = n(),
    events = sum(MACE_event),
    event_rate = round(mean(MACE_event) * 100, 2),
    median_followup = round(median(MACE_time), 2)
  )
print(event_stats)

# 创建并保存KM曲线
cat("\n=== 生成极端增强KM曲线 ===\n")

if (!dir.exists("Extreme_Test_Figures")) dir.create("Extreme_Test_Figures")

fit_test <- survfit(Surv(MACE_time, MACE_event) ~ LV_size_category, data = test_data)

# 执行对数秩检验
logrank_test <- survdiff(Surv(MACE_time, MACE_event) ~ LV_size_category, data = test_data)
cat("对数秩检验 p值：", round(1 - pchisq(logrank_test$chisq, df = length(logrank_test$n) - 1), 10), "\n")

km_plot_extreme <- ggsurvplot(
  fit_test,
  data = test_data,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  ggtheme = theme_bw(),
  palette = c("#1f77b4", "#d62728", "#ff7f0e"),  # 更对比鲜明的颜色
  legend.title = "LV Size Category",
  legend.labs = c("Normal LV", "Small LV", "Large LV"),
  xlab = "Follow-up Time (years)",
  ylab = "MACE-free Survival Probability", 
  title = "Extreme Test: Kaplan-Meier Curves with Enhanced Separation",
  font.title = 16,
  font.legend = 14,
  size = 2.0,  # 更粗的线条
  xlim = c(0, 8),
  ylim = c(0, 1),
  break.time.by = 1,
  break.y.by = 0.1,
  surv.scale = "percent"
)

# 保存测试图形
ggsave("Extreme_Test_Figures/Extreme_KM_MACE.png", 
       print(km_plot_extreme$plot), 
       width = 12, height = 10, dpi = 300, bg = "white")

cat("极端测试图形已保存：Extreme_Test_Figures/Extreme_KM_MACE.png\n")

# 打印详细生存曲线概要
cat("\n=== 详细生存曲线概要 ===\n")
summary_times <- c(1, 2, 3, 4, 5, 6)
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

cat("\n极端增强测试完成！如果此次曲线仍然分离不明显，可能需要进一步调整参数。\n")