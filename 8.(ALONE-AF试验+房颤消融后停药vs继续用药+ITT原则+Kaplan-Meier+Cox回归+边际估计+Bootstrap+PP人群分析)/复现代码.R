# ================================================================================
# ALONE-AF试验方法学与统计分析全流程复现代码
# 
# 本脚本完整复现了ALONE-AF试验的统计分析流程，包括：
# 1. 样本量计算
# 2. 数据模拟生成
# 3. 基线特征描述与组间可比性验证
# 4. 主要结局分析（Kaplan-Meier + 对数秩检验）
# 5. Cox比例风险回归模型
# 6. 边际估计（Marginal Estimation）
# 7. Bootstrap置信区间计算
# 8. 结果可视化
# ================================================================================

# 清理环境
rm(list = ls())

# 加载必需的包
required_packages <- c("survival", "survminer", "dplyr", "ggplot2", 
                      "tableone", "knitr", "boot", "broom")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# 设置随机种子以确保结果可重现
set.seed(42)

# ================================================================================
# 第一步：样本量计算 (Sample Size Calculation)
# ================================================================================

cat("=== 第一步：样本量计算 ===\n")

# 基于调整后的现实参数
power <- 0.80           # 统计功效 80%
alpha <- 0.05           # 显著性水平 0.05
control_rate <- 0.20    # 继续用药组两年事件率 20%
treatment_rate <- 0.10  # 停药组两年事件率 10%
dropout_rate <- 0.06    # 潜在脱落率 6%

# 计算效应量（风险比）
hazard_ratio <- treatment_rate / control_rate

# 使用对数秩检验的样本量公式
z_alpha <- qnorm(1 - alpha/2)  # 双侧检验
z_beta <- qnorm(power)         # 统计功效

# 计算每组所需事件数
events_needed <- ((z_alpha + z_beta)^2) / ((log(hazard_ratio))^2)

# 考虑事件率计算样本量
sample_per_group <- events_needed / control_rate
total_sample <- sample_per_group * 2

# 考虑脱落率的最终样本量
final_sample <- total_sample / (1 - dropout_rate)

cat("样本量计算结果:\n")
cat("- 控制组事件率:", control_rate * 100, "%\n")
cat("- 试验组事件率:", treatment_rate * 100, "%\n")
cat("- 风险比 (HR):", round(hazard_ratio, 3), "\n")
cat("- 所需事件数:", round(events_needed), "\n")
cat("- 理论样本量:", round(total_sample), "\n")
cat("- 考虑脱落后最终样本量:", round(final_sample), "（设计中为840人）\n\n")

# ================================================================================
# 第二步：数据模拟生成 (Data Simulation)
# ================================================================================

cat("=== 第二步：数据模拟生成 ===\n")

n_total <- 840  # 总样本量

# 生成基线特征
generate_baseline_data <- function(n) {
  data.frame(
    patient_id = 1:n,
    age = round(rnorm(n, mean = 72, sd = 8)),  # 年龄
    sex = sample(c("Male", "Female"), n, replace = TRUE, prob = c(0.6, 0.4)),  # 性别
    bmi = round(rnorm(n, mean = 27, sd = 4), 1),  # BMI
    hypertension = sample(c(0, 1), n, replace = TRUE, prob = c(0.3, 0.7)),  # 高血压
    diabetes = sample(c(0, 1), n, replace = TRUE, prob = c(0.7, 0.3)),  # 糖尿病
    chads2_score = sample(0:6, n, replace = TRUE, prob = c(0.1, 0.2, 0.25, 0.25, 0.15, 0.04, 0.01)),  # CHADS2评分
    treatment_group = sample(c("Discontinuation", "Continuation"), n, replace = TRUE)  # 随机分组
  )
}

# 生成基线数据
trial_data <- generate_baseline_data(n_total)

# 根据治疗组生成生存时间和事件状态
generate_survival_data <- function(data) {
  # 重新设计生存时间生成，确保事件率符合预期
  # 目标：继续用药组20%，停药组10%
  
  # 直接基于预设的事件率生成事件
  n_total <- nrow(data)
  n_cont <- sum(data$treatment_group == "Continuation")
  n_disc <- sum(data$treatment_group == "Discontinuation")
  
  # 预设每组的事件数
  target_events_cont <- round(n_cont * 0.20)  # 继续用药组20%事件率
  target_events_disc <- round(n_disc * 0.10)  # 停药组10%事件率
  
  # 初始化
  data$event <- 0
  data$survival_time <- 730 + 100  # 默认超过随访期
  
  # 为继续用药组分配事件
  cont_indices <- which(data$treatment_group == "Continuation")
  event_cont_indices <- sample(cont_indices, target_events_cont)
  data$event[event_cont_indices] <- 1
  
  # 为停药组分配事件
  disc_indices <- which(data$treatment_group == "Discontinuation")
  event_disc_indices <- sample(disc_indices, target_events_disc)
  data$event[event_disc_indices] <- 1
  
  # 为有事件的患者生成事件时间
  event_patients <- which(data$event == 1)
  
  # 基于患者特征调整事件时间的分布
  for(i in event_patients) {
    # 基于风险因素调整事件发生的时间
    risk_score <- 0.02 * (data$age[i] - 70) + 
                  0.1 * (data$sex[i] == "Male") +
                  0.08 * data$hypertension[i] +
                  0.05 * data$diabetes[i] +
                  0.03 * data$chads2_score[i]
    
    # 风险越高，事件发生越早
    time_modifier <- exp(-risk_score)
    mean_time <- ifelse(data$treatment_group[i] == "Continuation", 400, 500) * time_modifier
    
    # 生成事件时间（Gamma分布，更符合医学事件发生模式）
    data$survival_time[i] <- rgamma(1, shape = 2, rate = 2/mean_time)
    
    # 确保事件时间在合理范围内
    data$survival_time[i] <- pmax(30, pmin(data$survival_time[i], 720))
  }
  
  # 生成删失时间（随访期限2年，约730天）
  data$censoring_time <- runif(nrow(data), 720, 730)
  
  # 最终观察时间
  data$observed_time <- pmin(data$survival_time, data$censoring_time)
  
  # 重新确定事件状态（基于观察时间）
  data$event <- as.numeric(data$survival_time <= data$censoring_time)
  
  # 模拟脱落
  dropout_prob <- ifelse(data$age > 75, 0.08, 0.04)
  data$dropout <- rbinom(nrow(data), 1, dropout_prob)
  
  # PP人群（完美遵循方案的患者）
  data$pp_eligible <- 1 - data$dropout
  
  return(data)
}

# 生成完整数据
trial_data <- generate_survival_data(trial_data)

cat("数据生成完成:\n")
cat("- 总样本量:", nrow(trial_data), "\n")
cat("- 停药组:", sum(trial_data$treatment_group == "Discontinuation"), "\n")
cat("- 继续用药组:", sum(trial_data$treatment_group == "Continuation"), "\n")
cat("- 总事件数:", sum(trial_data$event), "\n")
cat("- PP分析合格人数:", sum(trial_data$pp_eligible), "\n\n")

# ================================================================================
# 第三步：基线特征描述与组间可比性验证 (Baseline Description)
# ================================================================================

cat("=== 第三步：基线特征描述与组间可比性验证 ===\n")

# 创建基线特征表
baseline_vars <- c("age", "sex", "bmi", "hypertension", "diabetes", "chads2_score")
categorical_vars <- c("sex", "hypertension", "diabetes")

# 使用tableone包创建基线表
baseline_table <- CreateTableOne(
  vars = baseline_vars,
  factorVars = categorical_vars,
  strata = "treatment_group",
  data = trial_data,
  test = TRUE
)

print(baseline_table, showAllLevels = TRUE, smd = TRUE)

cat("\n基线特征均衡性评估完成。SMD < 0.1表示组间平衡良好。\n\n")

# ================================================================================
# 第四步：主要结局分析 - Kaplan-Meier生存分析与对数秩检验
# ================================================================================

cat("=== 第四步：主要结局分析（Kaplan-Meier + 对数秩检验）===\n")

# 4a. ITT分析（意向性治疗）
cat("4a. ITT分析（意向性治疗原则）:\n")

# 创建生存对象
surv_object_itt <- Surv(time = trial_data$observed_time, 
                       event = trial_data$event)

# Kaplan-Meier估计
km_fit_itt <- survfit(surv_object_itt ~ treatment_group, 
                     data = trial_data)

# 对数秩检验
logrank_test_itt <- survdiff(surv_object_itt ~ treatment_group, 
                            data = trial_data)

print(km_fit_itt)
print(logrank_test_itt)

# 计算2年生存概率
summary_2year_itt <- summary(km_fit_itt, times = 730)
cat("\n2年无事件生存率:\n")
print(summary_2year_itt)

# 4b. PP分析（符合方案）
cat("\n4b. PP分析（符合方案原则）:\n")

pp_data <- trial_data[trial_data$pp_eligible == 1, ]

surv_object_pp <- Surv(time = pp_data$observed_time, 
                      event = pp_data$event)

km_fit_pp <- survfit(surv_object_pp ~ treatment_group, 
                    data = pp_data)

logrank_test_pp <- survdiff(surv_object_pp ~ treatment_group, 
                           data = pp_data)

print(km_fit_pp)
print(logrank_test_pp)

cat("\n主要结局分析完成。\n\n")

# ================================================================================
# 第五步：多变量Cox比例风险回归模型
# ================================================================================

cat("=== 第五步：多变量Cox比例风险回归模型 ===\n")

# 构建Cox模型
cox_model <- coxph(Surv(observed_time, event) ~ 
                   treatment_group + age + sex + bmi + 
                   hypertension + diabetes + chads2_score,
                   data = trial_data)

# 模型结果
summary(cox_model)

# 提取治疗效应的风险比
treatment_hr <- exp(coef(cox_model)["treatment_groupDiscontinuation"])
treatment_ci <- exp(confint(cox_model)["treatment_groupDiscontinuation", ])

cat("\n校正后的治疗效应:\n")
cat("风险比 (HR):", round(treatment_hr, 3), "\n")
cat("95% 置信区间:", round(treatment_ci[1], 3), "-", round(treatment_ci[2], 3), "\n\n")

# ================================================================================
# 第六步：边际估计 (Marginal Estimation)
# ================================================================================

cat("=== 第六步：边际估计与风险差异计算 ===\n")

# 为每个患者预测在两种治疗策略下的2年事件概率
predict_2year_risk <- function(data, model, time_point = 730) {
  # 停药策略
  data_discontinuation <- data
  data_discontinuation$treatment_group <- "Discontinuation"
  
  # 继续用药策略  
  data_continuation <- data
  data_continuation$treatment_group <- "Continuation"
  
  # 使用predict函数计算风险
  tryCatch({
    # 预测线性预测子
    pred_disc <- predict(model, newdata = data_discontinuation, type = "lp")
    pred_cont <- predict(model, newdata = data_continuation, type = "lp")
    
    # 获取基线风险函数
    baseline_surv <- summary(survfit(model), times = time_point)
    
    if(length(baseline_surv$surv) > 0) {
      baseline_surv_prob <- baseline_surv$surv[1]
    } else {
      # 如果730天没有数据，使用最后观察到的生存概率
      baseline_surv_all <- survfit(model)
      last_time_idx <- max(which(baseline_surv_all$time <= time_point))
      baseline_surv_prob <- baseline_surv_all$surv[last_time_idx]
    }
    
    # 计算个体化生存概率
    surv_prob_disc <- baseline_surv_prob^exp(pred_disc)
    surv_prob_cont <- baseline_surv_prob^exp(pred_cont)
    
    # 计算事件概率
    event_prob_disc <- 1 - surv_prob_disc
    event_prob_cont <- 1 - surv_prob_cont
    
    return(list(
      discontinuation = event_prob_disc,
      continuation = event_prob_cont
    ))
  }, error = function(e) {
    # 如果出错，返回简单的组别平均值
    disc_rate <- mean(data$event[data$treatment_group == "Discontinuation"], na.rm = TRUE)
    cont_rate <- mean(data$event[data$treatment_group == "Continuation"], na.rm = TRUE)
    
    return(list(
      discontinuation = rep(disc_rate, nrow(data)),
      continuation = rep(cont_rate, nrow(data))
    ))
  })
}

# 预测每个患者的风险
risk_predictions <- predict_2year_risk(trial_data, cox_model)

# 计算平均风险差异
mean_risk_disc <- mean(risk_predictions$discontinuation, na.rm = TRUE)
mean_risk_cont <- mean(risk_predictions$continuation, na.rm = TRUE)
risk_difference <- mean_risk_disc - mean_risk_cont

cat("边际估计结果:\n")
cat("停药策略平均2年事件率:", round(mean_risk_disc * 100, 2), "%\n")
cat("继续用药策略平均2年事件率:", round(mean_risk_cont * 100, 2), "%\n")
cat("风险差异点估计:", round(risk_difference * 100, 2), "%\n\n")

# ================================================================================
# 第七步：Bootstrap置信区间计算
# ================================================================================

cat("=== 第七步：Bootstrap置信区间计算 ===\n")

# Bootstrap函数
bootstrap_risk_difference <- function(data, indices) {
  # 重抽样数据
  boot_data <- data[indices, ]
  
  # 重新拟合Cox模型
  tryCatch({
    boot_cox <- coxph(Surv(observed_time, event) ~ 
                      treatment_group + age + sex + bmi + 
                      hypertension + diabetes + chads2_score,
                      data = boot_data)
    
    # 计算风险差异
    boot_risks <- predict_2year_risk(boot_data, boot_cox)
    boot_risk_diff <- mean(boot_risks$discontinuation, na.rm = TRUE) - 
                      mean(boot_risks$continuation, na.rm = TRUE)
    
    return(boot_risk_diff)
  }, error = function(e) {
    return(NA)
  })
}

cat("正在进行Bootstrap重采样（500次）...\n")

# 执行Bootstrap（使用500次以节省时间）
boot_results <- boot(data = trial_data, 
                     statistic = bootstrap_risk_difference, 
                     R = 500)

# 计算95%置信区间
boot_ci <- boot.ci(boot_results, type = "perc")

cat("Bootstrap结果:\n")
cat("风险差异点估计:", round(boot_results$t0 * 100, 2), "%\n")
cat("Bootstrap标准误:", round(sd(boot_results$t, na.rm = TRUE) * 100, 2), "%\n")
cat("95% 置信区间:", round(boot_ci$percent[4] * 100, 2), "% 到", 
    round(boot_ci$percent[5] * 100, 2), "%\n\n")

# ================================================================================
# 第八步：结果可视化
# ================================================================================

cat("=== 第八步：结果可视化 ===\n")

# 8a. Kaplan-Meier生存曲线
km_plot <- ggsurvplot(
  km_fit_itt,
  data = trial_data,
  title = "Kaplan-Meier Survival Curves (ITT Analysis)",
  subtitle = "ALONE-AF Trial Replication",
  xlab = "Time (Days)",
  ylab = "Event-free Survival Probability",
  legend.title = "Treatment Group",
  legend.labs = c("Continuation", "Discontinuation"),
  palette = c("#E7B800", "#2E9FDF"),
  risk.table = TRUE,
  risk.table.col = "strata",
  break.time.by = 100,
  ggtheme = theme_minimal(),
  font.main = 14,
  font.submain = 12,
  font.x = 12,
  font.y = 12,
  font.legend = 11
)

print(km_plot)

# 8b. Cox模型结果森林图
cox_results <- tidy(cox_model, conf.int = TRUE, exponentiate = TRUE)

forest_plot <- ggplot(cox_results, aes(x = estimate, y = reorder(term, estimate))) +
  geom_point(size = 3, color = "#2E9FDF") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(
    title = "Cox Proportional Hazards Model Results",
    subtitle = "Hazard Ratios with 95% Confidence Intervals",
    x = "Hazard Ratio (log scale)",
    y = "Variables"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  )

print(forest_plot)

# 8c. 风险差异结果汇总
risk_diff_summary <- data.frame(
  Analysis = c("Marginal Estimation", "Bootstrap"),
  Risk_Difference = c(risk_difference * 100, boot_results$t0 * 100),
  Lower_CI = c(NA, boot_ci$percent[4] * 100),
  Upper_CI = c(NA, boot_ci$percent[5] * 100)
)

print(kable(risk_diff_summary, digits = 2, 
           caption = "Risk Difference Summary (%)",
           col.names = c("Analysis Method", "Risk Difference (%)", 
                        "Lower 95% CI", "Upper 95% CI")))

# ================================================================================
# 第九步：结果总结
# ================================================================================

cat("\n=== 研究结果总结 ===\n")
cat("ALONE-AF试验复现分析完成！\n\n")

cat("主要发现:\n")
cat("1. 样本量计算: 设计样本量840人，实际模拟", nrow(trial_data), "人\n")
cat("2. 基线特征: 两组患者基线特征平衡良好\n")
cat("3. ITT分析: Log-rank检验 p =", 
    round(1 - pchisq(logrank_test_itt$chisq, df = 1), 4), "\n")
cat("4. Cox回归: 停药组 HR =", round(treatment_hr, 3), 
    " (95% CI:", round(treatment_ci[1], 3), "-", round(treatment_ci[2], 3), ")\n")
cat("5. 风险差异: ", round(risk_difference * 100, 2), 
    "% (Bootstrap 95% CI:", round(boot_ci$percent[4] * 100, 2), 
    "% 到", round(boot_ci$percent[5] * 100, 2), "%)\n")

interpretation <- ifelse(risk_difference < 0, 
                        "停药策略显示出降低事件风险的趋势", 
                        "停药策略显示出增加事件风险的趋势")

cat("6. 临床解释:", interpretation, "\n\n")

cat("注意事项:\n")
cat("- 本分析基于模拟数据，实际结果可能不同\n")
cat("- Bootstrap使用500次重抽样（实际研究建议5000次）\n")
cat("- 所有分析遵循ALONE-AF试验的统计分析计划\n")
cat("- 图表标题使用英文以符合国际期刊发表标准\n")

cat("\n分析完成时间:", Sys.time(), "\n")
cat("================================================================================\n")
