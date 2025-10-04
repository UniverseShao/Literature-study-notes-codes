################################################################################
# 测试KM曲线分离效果 - 快速验证脚本
# 
# 日期：2025-10-04
################################################################################

cat("开始测试KM曲线分离效果...\n")

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

# 使用新的更强效果分组标准
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
follow_up_years <- runif(n_test, min = 0.5, max = 12)

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

# 计算LV风险效应（大幅增强效果）
lv_risk_effect <- with(test_data, {
  ifelse(sex == "Male",
         ifelse(LVEDVi < 68, 0.045, ifelse(LVEDVi > 88, 0.035, 0)),
         ifelse(LVEDVi < 60, 0.045, ifelse(LVEDVi > 76, 0.035, 0))
  )
})

# 生成MACE事件（大幅增强事件率和组间差异）
mace_risk <- pmin(base_risk + lv_risk_effect * 5.0, 0.35)
test_data$MACE_event <- rbinom(n_test, 1, mace_risk * follow_up_years / 3)
test_data$MACE_time <- ifelse(test_data$MACE_event == 1,
                              runif(n_test, 0.5, follow_up_years),
                              follow_up_years)

# 打印分组统计
cat("\n=== 测试数据分组统计 ===\n")
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
cat("\n=== 生成KM曲线 ===\n")

if (!dir.exists("Test_Figures")) dir.create("Test_Figures")

fit_test <- survfit(Surv(MACE_time, MACE_event) ~ LV_size_category, data = test_data)

# 执行对数秩检验
logrank_test <- survdiff(Surv(MACE_time, MACE_event) ~ LV_size_category, data = test_data)
cat("对数秩检验 p值：", round(1 - pchisq(logrank_test$chisq, df = length(logrank_test$n) - 1), 4), "\n")

km_plot_test <- ggsurvplot(
  fit_test,
  data = test_data,
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
  title = "Test: Kaplan-Meier Curves for MACE by LV Size Category",
  font.title = 14,
  font.legend = 12,
  size = 1.2,
  xlim = c(0, 12),
  break.time.by = 2,
  surv.scale = "percent"
)

# 保存测试图形
ggsave("Test_Figures/Test_KM_MACE_Improved.png", 
       print(km_plot_test$plot), 
       width = 10, height = 8, dpi = 300, bg = "white")

cat("测试图形已保存：Test_Figures/Test_KM_MACE_Improved.png\n")

# 打印生存曲线概要
cat("\n=== 生存曲线概要 ===\n")
print(summary(fit_test, times = c(2, 5, 8, 10)))

cat("\n测试完成！请查看生成的图形以验证曲线分离效果。\n")