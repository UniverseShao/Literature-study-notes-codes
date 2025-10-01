# 0. 环境准备-------------------------------------------------------------------

# 安装和加载必要的R包
# 如果您已经安装过，可以跳过 install.packages() 的部分
# 设置CRAN镜像以避免安装错误
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# 尝试加载包，如果失败则安装
if (!require("survival")) {
  install.packages("survival")
  library(survival)
}
if (!require("lme4")) {
  install.packages("lme4")
  library(lme4)
}
if (!require("glmmTMB")) {
  install.packages("glmmTMB")
  library(glmmTMB)
}
if (!require("dplyr")) {
  install.packages("dplyr")
  library(dplyr)
}
if (!require("tibble")) {
  install.packages("tibble")
  library(tibble)
}

library(survival)
library(lme4)
library(glmmTMB)
library(dplyr)
library(tibble)

#创建一个结构上与研究类似的模拟数据集-------------------------------------------
set.seed(42) # 设置随机种子，确保结果可以复现

n_patients <- 2800 # 患者总数
n_centers <- 50    # 临床中心总数

sim_data <- tibble(
  patient_id = 1:n_patients,
  # 随机分配患者到不同的中心
  center_id = factor(sample(1:n_centers, n_patients, replace = TRUE)),
  # 随机分配治疗组
  treatment = factor(sample(c("ERC", "Usual_Care"), n_patients, replace = TRUE)),
  # 随机分配BMI分组
  bmi_group = factor(sample(c("Non_Obese", "Obese"), n_patients, replace = TRUE, prob = c(0.6, 0.4))),
  # 随机分配糖尿病分组
  diabetes = factor(sample(c("No", "Yes"), n_patients, replace = TRUE, prob = c(0.7, 0.3))),
  
  # 为Cox模型模拟生存数据（时间和事件状态）
  time_to_event = rweibull(n_patients, shape = 1.5, scale = 100),
  event_status = rbinom(n_patients, 1, 0.2), # 假设约20%的事件发生率
  
  # 为负二项模型模拟住院天数 (创建过度离散数据)
  hospital_nights = rnbinom(n_patients, mu = 8, size = 0.5), # size参数越小，数据越离散
  
  # 为线性混合模型模拟LVEF（连续变量）
  lvef = rnorm(n_patients, mean = 55, sd = 8),
  
  # 为逻辑斯蒂混合模型模拟2年时是否为窦性心律（二分类变量）
  sinus_rhythm_2y = rbinom(n_patients, 1, 0.6) # 假设60%的患者为窦性心律
)

# 预览生成的数据
cat("模拟数据前6行:\n")
print(head(sim_data))

# 1. 带脆弱性项的Cox回归模型----------------------------------------------------
# Surv(time, event) 定义生存结局
# treatment * bmi_group + treatment * diabetes 包含了主效应和交互效应
# frailty(center_id) 是校正中心效应的关键

cat("\n--- 正在运行带脆弱性项的Cox回归模型 ---\n")
cox_frailty_model <- coxph(
  Surv(time_to_event, event_status) ~ treatment * bmi_group + treatment * diabetes + frailty(center_id),
  data = sim_data
)

# 查看模型摘要
# 关注交互项 (如 treatmentUsual_Care:bmi_groupObese) 的P值
cat("\nCox回归模型结果摘要:\n")
print(summary(cox_frailty_model))

# 2. 负二项混合模型-------------------------------------------------------------
# hospital_nights 是计数型结果变量
# family = nbinom2 指定了负二项分布来处理过度离散
# + (1 | center_id) 是校正中心效应的关键，它将模型转变为混合模型

cat("\n--- 正在运行负二项混合模型 ---\n")
nb_mixed_model <- glmmTMB(
  hospital_nights ~ treatment * bmi_group + treatment * diabetes + (1 | center_id),
  data = sim_data,
  family = nbinom2
)

# 查看模型摘要
# 关注Fixed effects部分中交互项的P值
cat("\n负二项混合模型结果摘要:\n")
print(summary(nb_mixed_model))


# 3. 混合线性回归模型-----------------------------------------------------------
# lvef 是连续型结果变量
# + (1 | center_id) 是校正中心效应的关键

cat("\n--- 正在运行混合线性回归模型 ---\n")
linear_mixed_model <- lmer(
  lvef ~ treatment * bmi_group + treatment * diabetes + (1 | center_id),
  data = sim_data
)

# 查看模型摘要
# 关注Fixed effects部分中交互项的P值
cat("\n混合线性回归模型结果摘要:\n")
print(summary(linear_mixed_model))

# 4. 混合逻辑斯蒂回归模型-------------------------------------------------------

# sinus_rhythm_2y 是二分类结果变量 (0或1)
# family = binomial 指定进行逻辑斯蒂回归
# + (1 | center_id) 是校正中心效应的关键

cat("\n--- 正在运行混合逻辑斯蒂回归模型 ---\n")
logistic_mixed_model <- glmer(
  sinus_rhythm_2y ~ treatment * bmi_group + treatment * diabetes + (1 | center_id),
  data = sim_data,
  family = binomial
)

# 查看模型摘要
# 关注Fixed effects部分中交互项的P值
cat("\n混合逻辑斯蒂回归模型结果摘要:\n")
print(summary( logistic_mixed_model))


# 5. 可视化处理部分-------------------------------------------------------------

cat("\n--- 开始生成可视化图表 ---\n")

# 安装和加载可视化相关包
if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
}
if (!require("survminer")) {
  install.packages("survminer")
  library(survminer)
}
if (!require("gridExtra")) {
  install.packages("gridExtra")
  library(gridExtra)
}
if (!require("broom")) {
  install.packages("broom")
  library(broom)
}
if (!require("broom.mixed")) {
  install.packages("broom.mixed")
  library(broom.mixed)
}

library(ggplot2)
library(survminer)
library(gridExtra)
library(broom)
library(broom.mixed)

# 设置统一的主题
theme_set(theme_bw(base_size = 12))

# 5.1 Cox模型可视化 -----------------------------------------------------------

# 5.1.1 森林图 - 显示风险比和置信区间
cat("\n生成Cox模型森林图...\n")

# 提取Cox模型系数
cox_coef <- broom::tidy(cox_frailty_model, conf.int = TRUE, exponentiate = TRUE)
cox_coef <- cox_coef[!grepl("frailty", cox_coef$term), ]  # 移除frailty项

# 创建森林图
p1 <- ggplot(cox_coef, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_point(size = 3, color = "#2E86AB") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "#2E86AB") +
  labs(
    title = "Forest Plot: Hazard Ratios from Cox Regression Model",
    x = "Hazard Ratio (95% CI)",
    y = ""
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = seq(0, 2, 0.25))

# 5.1.2 Kaplan-Meier生存曲线
cat("生成Kaplan-Meier生存曲线...\n")

# 按治疗组分层
fit_treatment <- survfit(Surv(time_to_event, event_status) ~ treatment, data = sim_data)
p2 <- ggsurvplot(
  fit_treatment,
  data = sim_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  palette = c("#E63946", "#457B9D"),
  title = "Kaplan-Meier Survival Curves by Treatment Group",
  xlab = "Time to Event",
  ylab = "Survival Probability",
  legend.title = "Treatment",
  legend.labs = c("ERC", "Usual Care"),
  ggtheme = theme_bw()
)

# 按治疗和BMI组合分层
sim_data$treat_bmi <- interaction(sim_data$treatment, sim_data$bmi_group)
fit_interaction <- survfit(Surv(time_to_event, event_status) ~ treat_bmi, data = sim_data)
p3 <- ggsurvplot(
  fit_interaction,
  data = sim_data,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.3,
  palette = c("#E63946", "#F4A261", "#2A9D8F", "#264653"),
  title = "Kaplan-Meier Curves: Treatment × BMI Interaction",
  xlab = "Time to Event",
  ylab = "Survival Probability",
  legend.title = "Group",
  legend.labs = c("ERC + Non-Obese", "Usual Care + Non-Obese", 
                  "ERC + Obese", "Usual Care + Obese"),
  ggtheme = theme_bw()
)

# 5.2 负二项模型可视化 -------------------------------------------------------

cat("生成住院天数可视化...\n")

# 5.2.1 箱线图展示住院天数分布
p4 <- ggplot(sim_data, aes(x = treatment, y = hospital_nights, fill = bmi_group)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  scale_fill_manual(values = c("#E9C46A", "#E76F51")) +
  labs(
    title = "Distribution of Hospital Nights by Treatment and BMI",
    x = "Treatment Group",
    y = "Hospital Nights (days)",
    fill = "BMI Group"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

# 5.2.2 交互效应图
p5 <- sim_data %>%
  group_by(treatment, bmi_group) %>%
  summarise(
    mean_nights = mean(hospital_nights),
    se = sd(hospital_nights) / sqrt(n()),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = bmi_group, y = mean_nights, group = treatment, color = treatment)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_nights - se, ymax = mean_nights + se), width = 0.1) +
  scale_color_manual(values = c("#E63946", "#457B9D")) +
  labs(
    title = "Interaction Plot: Treatment × BMI on Hospital Nights",
    x = "BMI Group",
    y = "Mean Hospital Nights (± SE)",
    color = "Treatment"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

# 5.3 线性混合模型可视化 -----------------------------------------------------

cat("生成LVEF可视化...\n")

# 5.3.1 小提琴图展示LVEF分布
p6 <- ggplot(sim_data, aes(x = treatment, y = lvef, fill = diabetes)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#06D6A0", "#EF476F")) +
  labs(
    title = "Distribution of LVEF by Treatment and Diabetes Status",
    x = "Treatment Group",
    y = "Left Ventricular Ejection Fraction (%)",
    fill = "Diabetes"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

# 5.3.2 交互效应图
p7 <- sim_data %>%
  group_by(treatment, diabetes) %>%
  summarise(
    mean_lvef = mean(lvef),
    se = sd(lvef) / sqrt(n()),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = diabetes, y = mean_lvef, group = treatment, color = treatment)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_lvef - se, ymax = mean_lvef + se), width = 0.1) +
  scale_color_manual(values = c("#E63946", "#457B9D")) +
  labs(
    title = "Interaction Plot: Treatment × Diabetes on LVEF",
    x = "Diabetes Status",
    y = "Mean LVEF (± SE)",
    color = "Treatment"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

# 5.4 逻辑斯蒂模型可视化 -----------------------------------------------------

cat("生成窦性心律可视化...\n")

# 5.4.1 窦性心律比例条形图
p8 <- sim_data %>%
  group_by(treatment, diabetes) %>%
  summarise(
    sinus_rate = mean(sinus_rhythm_2y),
    n = n(),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = treatment, y = sinus_rate, fill = diabetes)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", sinus_rate * 100)),
            position = position_dodge(0.9), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("#06D6A0", "#EF476F")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.8)) +
  labs(
    title = "Sinus Rhythm Rate at 2 Years by Treatment and Diabetes",
    x = "Treatment Group",
    y = "Proportion with Sinus Rhythm",
    fill = "Diabetes"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

# 5.4.2 优势比森林图
cat("生成逻辑斯蒂模型森林图...\n")

logistic_coef <- broom.mixed::tidy(logistic_mixed_model, conf.int = TRUE, exponentiate = TRUE, effects = "fixed")

p9 <- ggplot(logistic_coef, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_point(size = 3, color = "#8338EC") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "#8338EC") +
  labs(
    title = "Forest Plot: Odds Ratios from Logistic Regression Model",
    x = "Odds Ratio (95% CI)",
    y = ""
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# 5.5 保存所有图表 -----------------------------------------------------------

cat("\n保存所有可视化图表...\n")

# 保存Cox森林图
ggsave("Figure1_Cox_Forest_Plot.png", plot = p1, width = 10, height = 6, dpi = 300)
cat("已保存: Figure1_Cox_Forest_Plot.png\n")

# 保存KM曲线 - 治疗组（使用特殊方法保存survminer对象）
png("Figure2_KM_Treatment.png", width = 3000, height = 2400, res = 300)
print(p2)
dev.off()
cat("已保存: Figure2_KM_Treatment.png\n")

# 保存KM曲线 - 交互效应（使用特殊方法保存survminer对象）
png("Figure3_KM_Interaction.png", width = 3600, height = 2400, res = 300)
print(p3)
dev.off()
cat("已保存: Figure3_KM_Interaction.png\n")

# 保存住院天数箱线图
ggsave("Figure4_Hospital_Nights_Boxplot.png", plot = p4, width = 10, height = 6, dpi = 300)
cat("已保存: Figure4_Hospital_Nights_Boxplot.png\n")

# 保存住院天数交互图
ggsave("Figure5_Hospital_Interaction.png", plot = p5, width = 8, height = 6, dpi = 300)
cat("已保存: Figure5_Hospital_Interaction.png\n")

# 保存LVEF小提琴图
ggsave("Figure6_LVEF_Violin.png", plot = p6, width = 10, height = 6, dpi = 300)
cat("已保存: Figure6_LVEF_Violin.png\n")

# 保存LVEF交互图
ggsave("Figure7_LVEF_Interaction.png", plot = p7, width = 8, height = 6, dpi = 300)
cat("已保存: Figure7_LVEF_Interaction.png\n")

# 保存窦性心律条形图
ggsave("Figure8_Sinus_Rhythm_Barplot.png", plot = p8, width = 10, height = 6, dpi = 300)
cat("已保存: Figure8_Sinus_Rhythm_Barplot.png\n")

# 保存逻辑斯蒂森林图
ggsave("Figure9_Logistic_Forest_Plot.png", plot = p9, width = 10, height = 6, dpi = 300)
cat("已保存: Figure9_Logistic_Forest_Plot.png\n")

# 创建综合图表
cat("\n创建综合对比图...\n")

# 综合展示：所有交互效应
p_combined <- grid.arrange(p5, p7, ncol = 2, 
                           top = "Treatment Interaction Effects Comparison")
ggsave("Figure10_Combined_Interactions.png", plot = p_combined, width = 14, height = 6, dpi = 300)
cat("已保存: Figure10_Combined_Interactions.png\n")

cat("\n========================================\n")
cat("✓ 所有可视化图表生成完成！\n")
cat("========================================\n")
cat("\n生成的图表列表：\n")
cat("  1. Figure1_Cox_Forest_Plot.png - Cox模型风险比森林图\n")
cat("  2. Figure2_KM_Treatment.png - KM生存曲线（按治疗分组）\n")
cat("  3. Figure3_KM_Interaction.png - KM生存曲线（治疗×BMI交互）\n")
cat("  4. Figure4_Hospital_Nights_Boxplot.png - 住院天数箱线图\n")
cat("  5. Figure5_Hospital_Interaction.png - 住院天数交互效应图\n")
cat("  6. Figure6_LVEF_Violin.png - LVEF分布小提琴图\n")
cat("  7. Figure7_LVEF_Interaction.png - LVEF交互效应图\n")
cat("  8. Figure8_Sinus_Rhythm_Barplot.png - 窦性心律比例条形图\n")
cat("  9. Figure9_Logistic_Forest_Plot.png - 逻辑斯蒂模型优势比森林图\n")
cat(" 10. Figure10_Combined_Interactions.png - 综合交互效应对比图\n")
cat("\n所有图表已保存至当前工作目录！\n")





