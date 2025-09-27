#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. 准备工作：加载所需的R包
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 如果您尚未安装这些包，请先运行 install.packages("...")
# install.packages(c("survival", "ggplot2", "dplyr", "patchwork"))

library(survival)  # 用于Cox回归模型
library(ggplot2)   # 用于数据可视化
library(dplyr)     # 用于数据处理
library(patchwork) # 用于将多个图组合在一起

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. 数据准备
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 从HR随age变化曲线.R中读取模拟数据
# 注意：这里我们假设数据已经存在，直接使用
set.seed(42)

# 创建基础数据集（与原脚本一致）
sim_data <- data.frame(
  id = 1:3000,
  age = runif(3000, 35, 95),
  treatment = sample(0:1, 3000, replace = TRUE)
)

# --- 定义生成事件时间的函数 ---
generate_survival_time <- function(data, beta_treat, beta_age, beta_interact, shape = 1.5) {
  linear_predictor <- beta_treat * data$treatment + 
                      beta_age * (data$age - 65) + 
                      beta_interact * data$treatment * (data$age - 65)
  u <- runif(nrow(data))
  lambda <- 0.005  # 恢复lambda值
  time <- (-log(u) / (lambda * exp(linear_predictor)))^(1 / shape)
  return(time)
}

# --- 模拟三个不同的终点事件 ---
# A. 疗效终点
sim_data$time_efficacy <- generate_survival_time(
  sim_data,
  beta_treat = 0.5,  # 增加治疗效应系数
  beta_age = 0.02,
  beta_interact = -0.03
)
sim_data$status_efficacy <- ifelse(sim_data$time_efficacy <= 10, 1, 0)

# B. 安全性终点
sim_data$time_safety <- generate_survival_time(
  sim_data,
  beta_treat = -0.8,
  beta_age = 0.01,
  beta_interact = 0.05
)
sim_data$status_safety <- ifelse(sim_data$time_safety <= 10, 1, 0)

# C. 净获益终点
sim_data$time_net <- generate_survival_time(
  sim_data,
  beta_treat = -0.3,
  beta_age = 0.03,
  beta_interact = -0.02
)
sim_data$status_net <- ifelse(sim_data$time_net <= 10, 1, 0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. 绘制累计风险曲线
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# --- 创建一个通用的绘图函数 ---
plot_cumulative_incidence <- function(data, time_var, status_var, title, age_group) {
  # 根据年龄组过滤数据
  filtered_data <- data %>%
    filter(
      if(age_group == "<70") age < 70 else 
      if(age_group == "70-74") age >= 70 & age <= 74 else 
      if(age_group == "75-79") age >= 75 & age <= 79 else 
      age >= 80
    )
  
  # 检查过滤后的数据是否为空
  if(nrow(filtered_data) == 0) {
    stop(paste("No data for age group:", age_group))
  }
  
  # 确保数据中有两种治疗方法
  # 检查treatment变量的分布
  treatment_counts <- table(filtered_data$treatment)
  print(paste("Treatment distribution:", paste(treatment_counts, collapse=", ")))
  
  # 确保从时间点0开始
  times_seq <- c(0, seq(1, min(36, max(filtered_data[[time_var]], na.rm=TRUE)), by = 1))
  
  # 为两种治疗方法分别创建数据集
  mono_data <- filtered_data %>% filter(treatment == 0)
  combo_data <- filtered_data %>% filter(treatment == 1)
  
  # 分别计算两种治疗方法的生存曲线
  surv_fit_mono <- survfit(Surv(get(time_var), get(status_var)) ~ 1, data = mono_data)
  surv_fit_combo <- survfit(Surv(get(time_var), get(status_var)) ~ 1, data = combo_data)
  
  # 提取结果
  results_mono <- summary(surv_fit_mono, times = times_seq)
  results_combo <- summary(surv_fit_combo, times = times_seq)
  
  # 获取生存率数据
  surv_mono <- results_mono$surv
  surv_combo <- results_combo$surv
  
  # 确保两个结果有相同的时间点
  common_times <- intersect(results_mono$time, results_combo$time)
  
  # 构建数据框用于绘图
  plot_data <- rbind(
    data.frame(
      time = results_mono$time,
      cumulative_incidence = (1 - surv_mono) * 100,  # 转换为百分比
      treatment = "Rivaroxaban monotherapy",
      group = "monotherapy"
    ),
    data.frame(
      time = results_combo$time,
      cumulative_incidence = (1 - surv_combo) * 100,  # 转换为百分比
      treatment = "Rivaroxaban + antiplatelet agent therapy",
      group = "combination"
    )
  )

  # 绘图
  p <- ggplot(plot_data, aes(x = time, y = cumulative_incidence, color = treatment)) +
    geom_step(size = 1) +
    labs(
      title = paste(title, age_group),
      x = "Time, mo",
      y = "Cumulative incidence, %"
    ) +
    scale_y_continuous(limits = c(0, 60), expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(0, 36, by = 6)) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      legend.position = "top",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )

  return(p)
}

# --- 生成四张疗效终点图 ---
plot_a_eff <- plot_cumulative_incidence(sim_data, "time_efficacy", "status_efficacy", "Primary efficacy end point", "<70")
plot_b_eff <- plot_cumulative_incidence(sim_data, "time_efficacy", "status_efficacy", "Primary efficacy end point", "70-74")
plot_c_eff <- plot_cumulative_incidence(sim_data, "time_efficacy", "status_efficacy", "Primary efficacy end point", "75-79")
plot_d_eff <- plot_cumulative_incidence(sim_data, "time_efficacy", "status_efficacy", "Primary efficacy end point", "≥80")

# --- 生成四张安全性终点图 ---
plot_a_saf <- plot_cumulative_incidence(sim_data, "time_safety", "status_safety", "Primary safety end point", "<70")
plot_b_saf <- plot_cumulative_incidence(sim_data, "time_safety", "status_safety", "Primary safety end point", "70-74")
plot_c_saf <- plot_cumulative_incidence(sim_data, "time_safety", "status_safety", "Primary safety end point", "75-79")
plot_d_saf <- plot_cumulative_incidence(sim_data, "time_safety", "status_safety", "Primary safety end point", "≥80")

# --- 组合图形 ---
final_plot_eff <- plot_a_eff / plot_b_eff / plot_c_eff / plot_d_eff
final_plot_saf <- plot_a_saf / plot_b_saf / plot_c_saf / plot_d_saf

# 显示最终的组合图
final_plot_eff
final_plot_saf

# --- 保存为PDF文件 ---
pdf("累计分线曲线.pdf", width = 10, height = 12)
print(final_plot_eff)
print(final_plot_saf)
dev.off()