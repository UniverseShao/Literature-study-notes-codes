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
# 2. 数据模拟
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 我们的目标是创建一个数据集，其中治疗效果与年龄存在交互作用
# 从而复现图A, B, C中的趋势

set.seed(42) # 设置随机种子以保证结果可重复
n <- 3000    # 模拟的患者数量

# 创建基础数据集
sim_data <- data.frame(
  id = 1:n,
  # 年龄均匀分布在35至95岁之间
  age = runif(n, 35, 95),
  # 随机分配治疗方案 (0 = 联合治疗, 1 = 单药治疗)
  treatment = sample(0:1, n, replace = TRUE)
)

# --- 定义生成事件时间的函数 ---
# 这个函数的核心是Cox模型的线性预测部分 (log-hazard)
# log(h(t)) = log(h0(t)) + beta_treat * treatment + beta_age * age + beta_interact * treatment * age
# 生存时间可以通过 Weibull 分布的逆CDF从这个log-hazard生成
generate_survival_time <- function(data, beta_treat, beta_age, beta_interact, shape = 1.5) {
  # 计算线性预测部分
  linear_predictor <- beta_treat * data$treatment + 
                      beta_age * (data$age - 65) + # 中心化年龄以提高模型稳定性
                      beta_interact * data$treatment * (data$age - 65)
  
  # 从Weibull分布生成生存时间
  # T = (-log(U) / (lambda * exp(LP)))^(1/shape)
  u <- runif(nrow(data))
  lambda <- 0.005 # 基线风险参数
  time <- (-log(u) / (lambda * exp(linear_predictor)))^(1 / shape)
  return(time)
}

# --- 模拟三个不同的终点事件 ---

# A. 疗效终点 (Primary Efficacy Endpoint)
# 趋势：单药治疗在高龄时更优 (HR随年龄增长而下降)
# 这意味着交互项系数 (beta_interact) 必须为负
sim_data$time_efficacy <- generate_survival_time(
  sim_data,
  beta_treat = 0.2,    # 在65岁时，单药治疗风险稍高
  beta_age = 0.02,     # 年龄本身是风险因素
  beta_interact = -0.03  # 关键：交互项为负，使单药效果随年龄增长而变好
)
sim_data$status_efficacy <- ifelse(sim_data$time_efficacy <= 10, 1, 0) # 假设10年随访期

# B. 安全性终点 (Primary Safety Endpoint)
# 趋势：单药治疗在低龄时更优 (HR随年龄增长而上升)
# 这意味着交互项系数 (beta_interact) 必须为正
sim_data$time_safety <- generate_survival_time(
  sim_data,
  beta_treat = -0.8,   # 在65岁时，单药治疗风险显著更低
  beta_age = 0.01,     # 年龄本身是风险因素
  beta_interact = 0.05   # 关键：交互项为正，使单药安全性优势随年龄增长而减弱
)
sim_data$status_safety <- ifelse(sim_data$time_safety <= 10, 1, 0)

# C. 净获益终点 (Net MACE + Major Bleeding)
# 趋势：单药治疗在所有年龄段风险都较低，但在高龄时优势更显著
# 这意味着交互项系数为负，且主效应也为负
sim_data$time_net <- generate_survival_time(
  sim_data,
  beta_treat = -0.3,   # 在65岁时，单药治疗有优势
  beta_age = 0.03,     # 年龄是强风险因素
  beta_interact = -0.02  # 关键：交互项为负，使单药优势随年龄增长而更强
)
sim_data$status_net <- ifelse(sim_data$time_net <= 10, 1, 0)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. 建模与绘图
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- 创建一个通用的绘图函数 ---
# 这个函数将自动化建模、预测和绘图的过程
plot_hr_interaction <- function(data, time_var, status_var, title) {
  
  # 1. 拟合Cox模型，包含治疗与年龄的交互项
  formula <- as.formula(sprintf("Surv(%s, %s) ~ treatment * age", time_var, status_var))
  model <- coxph(formula, data = data)
  
  # 2. 创建用于预测的新数据集
  # 我们将在35到95岁之间，每隔1岁计算一次HR
  pred_grid <- data.frame(
    age = seq(35, 95, by = 1),
    treatment = 1  # 我们关心的是单药治疗(1)相对于对照(0)的效果
  )
  
  # 3. 进行预测
  # `predict`可以计算log(HR)及其标准误
  # `type="lp"`给出了线性预测值(log-hazard)
  # `reference="zero"`可以计算与所有协变量为0时的差异，但更手动的方法更清晰
  
  # 手动计算log(HR)和SE，以获得最大的清晰度和控制
  # log(HR(age)) = (β_treat*1 + β_age*age + β_interact*1*age) - (β_treat*0 + β_age*age + β_interact*0*age)
  #             = β_treat + β_interact * age
  
  coefs <- coef(model)
  vcov_mat <- vcov(model)
  
  # 提取所需的系数和方差/协方差
  b_treat <- coefs["treatment"]
  b_interact <- coefs["treatment:age"]
  var_treat <- vcov_mat["treatment", "treatment"]
  var_interact <- vcov_mat["treatment:age", "treatment:age"]
  cov_treat_interact <- vcov_mat["treatment", "treatment:age"]
  
  # 计算每个年龄点的log(HR)和标准误
  pred_results <- pred_grid %>%
    mutate(
      log_hr = b_treat + b_interact * age,
      # Var(aX + bY) = a^2*Var(X) + b^2*Var(Y) + 2ab*Cov(X,Y)
      # 这里 Var(log_hr) = Var(b_treat + b_interact*age)
      se_log_hr = sqrt(var_treat + age^2 * var_interact + 2 * age * cov_treat_interact),
      
      # 计算HR和95%置信区间
      hr = exp(log_hr),
      ci_low = exp(log_hr - 1.96 * se_log_hr),
      ci_high = exp(log_hr + 1.96 * se_log_hr)
    )

  # 4. 使用ggplot2绘图
  p <- ggplot(pred_results, aes(x = age, y = hr)) +
    # 添加95% CI的带状区域
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "darkcyan", alpha = 0.3) +
    # 添加HR的点估计曲线
    geom_line(color = "darkcyan", linewidth = 1) +
    # 添加HR=1的参考线
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    # **关键：使用对数刻度Y轴，以匹配原始图形的视觉效果**
    scale_y_log10(
      breaks = c(0.01, 0.02, 0.04, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 10),
      labels = c("0.01","0.02","0.04","0.08","0.1","0.2","0.4","0.6","0.8","1", "2", "4", "6", "8", "10")
    ) +
    scale_x_continuous(breaks = seq(35, 95, by = 10)) +
    labs(
      title = title,
      x = "Age, y",
      y = "HR of rivaroxaban to rivaroxaban plus\nantiplatelet therapy (95% CI)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      axis.title.y = element_text(size = 9)
    )
  
  return(p)
}

# --- 生成三张图 ---
plot_a <- plot_hr_interaction(sim_data, "time_efficacy", "status_efficacy", "A  Primary efficacy end point")
plot_b <- plot_hr_interaction(sim_data, "time_safety", "status_safety", "B  Primary safety end point")
plot_c <- plot_hr_interaction(sim_data, "time_net", "status_net", "C  Net MACE + major bleeding")

# --- 使用patchwork包将三张图组合在一起 ---
final_plot <- plot_a / plot_b / plot_c

# 显示最终的组合图
final_plot

# --- 保存为PDF文件 ---
# 设置PDF设备并保存图形
pdf("HR随age变化曲线.pdf", width = 8, height = 10)
print(final_plot)
dev.off()

# --- 添加模型结果输出 ---
# 显示三个模型的详细结果

# A. 疗效终点模型结果
print("A. 主要疗效终点模型结果")
model_efficacy <- coxph(Surv(time_efficacy, status_efficacy) ~ treatment * age, data = sim_data)
print(summary(model_efficacy))

# B. 安全性终点模型结果
print("B. 主要安全性终点模型结果")
model_safety <- coxph(Surv(time_safety, status_safety) ~ treatment * age, data = sim_data)
print(summary(model_safety))

# C. 净获益终点模型结果
print("C. 净获益终点模型结果")
model_net <- coxph(Surv(time_net, status_net) ~ treatment * age, data = sim_data)
print(summary(model_net))