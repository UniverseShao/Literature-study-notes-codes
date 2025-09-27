#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 森林图.R
# 用于绘制与文献中相同的森林图
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 加载所需的R包
library(survival)  # 用于Cox回归模型
library(ggplot2)   # 用于数据可视化
library(dplyr)     # 用于数据处理
library(gridExtra) # 用于组合多个图形

# --- 重新生成与HR随age变化曲线.R中完全相同的模拟数据 ---
set.seed(42)
n <- 3000

sim_data <- data.frame(
  id = 1:n,
  age = runif(n, 35, 95),
  treatment = sample(0:1, n, replace = TRUE)
)

# 定义生成事件时间的函数
generate_survival_time <- function(data, beta_treat, beta_age, beta_interact, shape = 1.5) {
  linear_predictor <- beta_treat * data$treatment + 
                      beta_age * (data$age - 65) + 
                      beta_interact * data$treatment * (data$age - 65)
  u <- runif(nrow(data))
  lambda <- 0.005
  time <- (-log(u) / (lambda * exp(linear_predictor)))^(1 / shape)
  return(time)
}

# 模拟三个不同的终点事件
sim_data$time_efficacy <- generate_survival_time(
  sim_data,
  beta_treat = 0.2,
  beta_age = 0.02,
  beta_interact = -0.03
)
sim_data$status_efficacy <- ifelse(sim_data$time_efficacy <= 10, 1, 0)

sim_data$time_safety <- generate_survival_time(
  sim_data,
  beta_treat = -0.8,
  beta_age = 0.01,
  beta_interact = 0.05
)
sim_data$status_safety <- ifelse(sim_data$time_safety <= 10, 1, 0)

sim_data$time_net <- generate_survival_time(
  sim_data,
  beta_treat = -0.3,
  beta_age = 0.03,
  beta_interact = -0.02
)
sim_data$status_net <- ifelse(sim_data$time_net <= 10, 1, 0)

# --- 创建一个函数来计算每个年龄组的HR和CI ---
calculate_hr_by_age_group <- function(data, time_var, status_var, age_groups) {
  # 创建一个空的数据框来存储结果
  results <- data.frame()
  
  # 遍历每个年龄组
  for (group in age_groups) {
    # 定义年龄范围
    if (group == "<70") {
      age_filter <- data$age < 70
    } else if (group == "70-74") {
      age_filter <- data$age >= 70 & data$age <= 74
    } else if (group == "75-79") {
      age_filter <- data$age >= 75 & data$age <= 79
    } else if (group == ">=80") {
      age_filter <- data$age >= 80
    }
    
    # 提取当前年龄组的数据
    group_data <- data[age_filter, ]
    
    # 拟合Cox模型
    model <- coxph(Surv(get(time_var), get(status_var)) ~ treatment, data = group_data)
    
    # 提取HR和CI
    hr <- exp(coef(model)["treatment"])
    ci_low <- exp(confint(model)["treatment", "2.5 %"])
    ci_high <- exp(confint(model)["treatment", "97.5 %"])
    p_value <- coef(summary(model))["treatment", "Pr(>|z|)"]
    
    # 添加到结果数据框
    results <- rbind(results, data.frame(
      age_group = group,
      hr = hr,
      ci_low = ci_low,
      ci_high = ci_high,
      p_value = p_value
    ))
  }
  
  return(results)
}

# --- 计算三个终点事件的HR ---
age_groups <- c("<70", "70-74", "75-79", ">=80")

# A. 疗效终点
efficacy_results <- calculate_hr_by_age_group(sim_data, "time_efficacy", "status_efficacy", age_groups)

# B. 安全性终点
safety_results <- calculate_hr_by_age_group(sim_data, "time_safety", "status_safety", age_groups)

# C. 净获益终点
net_results <- calculate_hr_by_age_group(sim_data, "time_net", "status_net", age_groups)

# --- 创建森林图 ---

# 通用的森林图绘图函数
plot_forest_plot <- function(data, title, y_label) {
  # 创建ggplot对象
  p <- ggplot(data, aes(x = hr, y = age_group)) +
    geom_point(size = 3) +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
    scale_x_log10(breaks = c(0.1, 0.5, 1, 2), labels = c("0.1", "0.5", "1", "2")) +
    labs(title = title, x = y_label, y = "") +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_text(size = 9),
      axis.text.y = element_text(size = 9)
    )
  
  return(p)
}

# --- 绘制三个森林图 ---

# A. 疗效终点森林图
forest_plot_a <- plot_forest_plot(efficacy_results, "Primary efficacy end point", "HR (95% CI)")

# B. 安全性终点森林图
forest_plot_b <- plot_forest_plot(safety_results, "Primary safety end point", "HR (95% CI)")

# C. 净获益终点森林图
forest_plot_c <- plot_forest_plot(net_results, "Net MACE + major bleeding", "HR (95% CI)")

# --- 使用gridExtra包组合三个图 ---
# --- 将森林图保存为PDF文件 ---
# 设置PDF输出参数
pdf("森林图.pdf", width = 8, height = 10)
final_forest_plot <- grid.arrange(forest_plot_a, forest_plot_b, forest_plot_c, ncol = 1)
dev.off()

# 显示最终的组合图
print(final_forest_plot)