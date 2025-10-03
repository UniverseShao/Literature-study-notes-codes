# -----------------------------------------------------------------------------
# 步骤 0: 加载必要的R包
# -----------------------------------------------------------------------------
library(tidyverse)
library(survival)
library(rms)       # 用于样条Cox回归 (cph, rcs, Predict)
library(patchwork) # 用于拼接图形

# -----------------------------------------------------------------------------
# 步骤 1: 创建一个模拟的、逼真的数据集
# -----------------------------------------------------------------------------
set.seed(42)
n_patients <- 3000

sim_data <- tibble(
  id = 1:n_patients,
  
  # 模拟连续预测变量 LVEDD
  LVEDD = rnorm(n_patients, mean = 48, sd = 8),
  
  # 创建6个二分类亚组变量
  sex = factor(sample(c("Male", "Female"), n_patients, replace = TRUE, prob = c(0.6, 0.4))),
  age_group = factor(ifelse(rnorm(n_patients, 65, 10) >= 65, "Age >=65 Years", "Age <65 Years")),
  bmi_group = factor(ifelse(rnorm(n_patients, 25, 3) >= 24, "BMI >=24kg/m^2", "BMI <24kg/m^2")),
  ablation = factor(sample(c("Ablation", "No-Ablation"), n_patients, replace = TRUE, prob = c(0.4, 0.6))),
  diabetes = factor(sample(c("DM", "No-DM"), n_patients, replace = TRUE, prob = c(0.3, 0.7))),
  lv_hypertrophy = factor(sample(c("LV Hypertrophy", "No-LV Hypertrophy"), n_patients, replace = TRUE, prob = c(0.2, 0.8))),
  
  # 模拟一个真实的"U型"风险 + 与性别有轻微交互作用
  true_log_hazard = 0.005 * (LVEDD - 48)^2 + 
                    0.2 * (sex == "Female") + 
                    0.002 * (LVEDD - 48)^2 * (sex == "Female"),
  
  # 模拟生存时间
  time_to_event = rexp(n_patients, rate = 0.05 * exp(true_log_hazard)),
  censoring_time = runif(n_patients, 1, 10),
  observed_time = pmin(time_to_event, censoring_time),
  event_status = as.numeric(time_to_event <= censoring_time)
)

# 使用rms包需要先设定数据分布环境
ddist <- datadist(sim_data)
options(datadist = 'ddist')

# -----------------------------------------------------------------------------
# 步骤 2: 创建LVEDD分布数据（用于柱状图背景）
# -----------------------------------------------------------------------------
# 计算LVEDD的直方图数据
lvedd_min <- floor(min(sim_data$LVEDD))
lvedd_max <- ceiling(max(sim_data$LVEDD))
lvedd_breaks <- seq(lvedd_min, lvedd_max, by = 2)  # 每2mm一个bin
lvedd_hist <- hist(sim_data$LVEDD, breaks = lvedd_breaks, plot = FALSE)

# 创建柱状图数据框
hist_data <- tibble(
  midpoint = lvedd_hist$mids,
  percentage = (lvedd_hist$counts / sum(lvedd_hist$counts)) * 100
)

# -----------------------------------------------------------------------------
# 步骤 3: 创建一个可复用的、强大的绘图函数
# -----------------------------------------------------------------------------

#' @title 创建亚组交互作用的样条图（带LVEDD分布柱状图背景）
#' @param data 数据集
#' @param subgroup_var 亚组变量的字符串名, e.g., "sex"
#' @param subgroup_labels 图例标签, e.g., c("Female", "Male")
#' @param subgroup_colors 曲线颜色, e.g., c("blue", "red")
#' @param plot_title 图的标题, e.g., "A"
#' @param hist_data LVEDD分布柱状图数据
#' @return 一个ggplot对象

create_subgroup_plot <- function(data, subgroup_var, subgroup_labels, subgroup_colors, plot_title, hist_data) {
  
  # 将字符串变量名转换为公式中可用的符号
  subgroup_sym <- as.symbol(subgroup_var)
  
  # --- 模型拟合与P值计算 ---
  
  # 模型1: 包含交互项 (用于绘图和计算交互P值)
  # 使用cph(coxph的增强版)和rcs(限制性立方样条), 设定4个节点
  formula_interaction <- as.formula(paste("Surv(observed_time, event_status) ~ rcs(LVEDD, 4) *", subgroup_var))
  model_interaction <- cph(formula_interaction, data = data, x=TRUE, y=TRUE)
  
  # 模型2: 不含交互项 (用于计算Overall P和Nonlinear P)
  formula_main <- as.formula(paste("Surv(observed_time, event_status) ~ rcs(LVEDD, 4) +", subgroup_var))
  model_main <- cph(formula_main, data = data, x=TRUE, y=TRUE)
  
  # 使用 anova() 比较两个模型来获得交互作用P值
  anova_interaction <- anova(model_interaction)
  # 查找包含交互项的行（通常是最后一行或包含" * "的行）
  interaction_rows <- grep("\\*", rownames(anova_interaction))
  if(length(interaction_rows) > 0) {
    p_interaction <- anova_interaction[tail(interaction_rows, 1), "P"]
  } else {
    p_interaction <- NA
  }
  
  # 从主效应模型中提取 LVEDD 的 P 值
  p_values_main <- anova(model_main)
  # 查找LVEDD相关的行
  lvedd_rows <- grep("LVEDD", rownames(p_values_main), ignore.case = TRUE)
  nonlinear_rows <- grep("Nonlinear", rownames(p_values_main), ignore.case = TRUE)
  
  p_overall <- if(length(lvedd_rows) > 0) p_values_main[lvedd_rows[1], "P"] else NA
  p_nonlinear <- if(length(nonlinear_rows) > 0) p_values_main[nonlinear_rows[1], "P"] else NA
  
  # 保存统计结果供后续输出
  stats_results <- list(
    subgroup = plot_title,
    subgroup_var = subgroup_var,
    p_overall = p_overall,
    p_nonlinear = p_nonlinear,
    p_interaction = p_interaction,
    model_interaction = model_interaction,
    model_main = model_main
  )
  
  # --- 数据预测 ---
  
  # 使用rms::Predict函数从交互模型中生成预测值 (HR和95% CI)
  # 这是最关键的一步！
  # 根据亚组变量动态生成预测数据
  if(subgroup_var == "sex") {
    pred_data <- Predict(model_interaction, LVEDD, sex, fun = exp, ref.zero = TRUE) %>% as_tibble()
  } else if(subgroup_var == "age_group") {
    pred_data <- Predict(model_interaction, LVEDD, age_group, fun = exp, ref.zero = TRUE) %>% as_tibble()
  } else if(subgroup_var == "bmi_group") {
    pred_data <- Predict(model_interaction, LVEDD, bmi_group, fun = exp, ref.zero = TRUE) %>% as_tibble()
  } else if(subgroup_var == "ablation") {
    pred_data <- Predict(model_interaction, LVEDD, ablation, fun = exp, ref.zero = TRUE) %>% as_tibble()
  } else if(subgroup_var == "diabetes") {
    pred_data <- Predict(model_interaction, LVEDD, diabetes, fun = exp, ref.zero = TRUE) %>% as_tibble()
  } else if(subgroup_var == "lv_hypertrophy") {
    pred_data <- Predict(model_interaction, LVEDD, lv_hypertrophy, fun = exp, ref.zero = TRUE) %>% as_tibble()
  }
  
  # --- 绘制图形（带柱状图背景和双Y轴） ---
  
  # 确定Y轴范围
  max_hr <- max(pred_data$upper, na.rm = TRUE)
  y_limit <- min(max(ceiling(max_hr), 3), 8)  # 最大不超过8
  
  # 确定柱状图的缩放比例（使柱状图最高不超过Y轴的1/4）
  max_percentage <- max(hist_data$percentage)
  scale_factor <- y_limit * 0.25 / max_percentage
  
  # 创建注释文本
  annotation_text <- paste0(
    "P nonlinear = ", format.pval(p_nonlinear, digits = 3, eps = 0.001), "\n",
    "P overall = ", format.pval(p_overall, digits = 3, eps = 0.001), "\n\n",
    "P for interaction = ", format.pval(p_interaction, digits = 3, eps = 0.001)
  )
  
  # 使用ggplot2绘图
  p <- ggplot() +
    # 1. 绘制柱状图（背景层）
    geom_col(data = hist_data, 
             aes(x = midpoint, y = percentage * scale_factor),
             fill = "#6BAED6", alpha = 0.6, width = 1.8) +
    
    # 2. 绘制HR=1的参考线
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
    
    # 3. 绘制95%置信区间（阴影）- 保留亚组差异
    geom_ribbon(data = pred_data,
                aes(x = LVEDD, ymin = lower, ymax = upper, 
                    fill = .data[[subgroup_var]], group = .data[[subgroup_var]]),
                alpha = 0.15) +
    
    # 4. 绘制HR估计值曲线（实线）- 保留亚组差异
    geom_line(data = pred_data,
              aes(x = LVEDD, y = yhat, color = .data[[subgroup_var]], 
                  group = .data[[subgroup_var]]),
              linewidth = 1.2) +
    
    # 5. 设置颜色
    scale_color_manual(values = subgroup_colors, labels = subgroup_labels, name = NULL) +
    scale_fill_manual(values = subgroup_colors, labels = subgroup_labels, name = NULL) +
    
    # 6. 添加统计注释
    annotate("text", x = Inf, y = Inf, 
             label = annotation_text, 
             hjust = 1.05, vjust = 1.3, 
             size = 3.5, lineheight = 0.9) +
    
    # 7. 设置坐标轴
    scale_x_continuous(
      name = "LVEDD (mm)",
      breaks = seq(30, 70, by = 10),
      limits = c(lvedd_min - 2, lvedd_max + 2)
    ) +
    scale_y_continuous(
      name = "HR",
      limits = c(0, y_limit),
      breaks = seq(0, y_limit, by = 1),
      sec.axis = sec_axis(
        ~ . / scale_factor, 
        name = "Percentage of Population (%)",
        breaks = seq(0, 20, by = 5)
      )
    ) +
    
    # 8. 添加标题和主题
    labs(title = plot_title) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      legend.position = c(0.15, 0.15),  # 图例在左下角
      legend.background = element_blank(),
      legend.key.size = unit(0.8, "lines"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 10, r = 15, b = 10, l = 10, unit = "pt")
    )
  
  # 返回图形和统计结果
  return(list(plot = p, stats = stats_results))
}

# -----------------------------------------------------------------------------
# 步骤 4: 调用函数生成6个图并拼接
# -----------------------------------------------------------------------------

# A: 按性别分组
result_A <- create_subgroup_plot(sim_data, "sex", c("Female", "Male"), c("#377eb8", "#e41a1c"), "A", hist_data)

# B: 按年龄分组
result_B <- create_subgroup_plot(sim_data, "age_group", c("Age <65 Years", "Age >=65 Years"), c("#377eb8", "#e41a1c"), "B", hist_data)

# C: 按BMI分组
result_C <- create_subgroup_plot(sim_data, "bmi_group", c("BMI <24kg/m^2", "BMI >=24kg/m^2"), c("#377eb8", "#e41a1c"), "C", hist_data)

# D: 按消融分组
result_D <- create_subgroup_plot(sim_data, "ablation", c("No-Ablation", "Ablation"), c("#377eb8", "#e41a1c"), "D", hist_data)

# E: 按糖尿病分组
result_E <- create_subgroup_plot(sim_data, "diabetes", c("No-DM", "DM"), c("#377eb8", "#e41a1c"), "E", hist_data)

# F: 按左心室肥厚分组
result_F <- create_subgroup_plot(sim_data, "lv_hypertrophy", c("No-LV Hypertrophy", "LV Hypertrophy"), c("#377eb8", "#e41a1c"), "F", hist_data)

# 提取图形对象
pA <- result_A$plot
pB <- result_B$plot
pC <- result_C$plot
pD <- result_D$plot
pE <- result_E$plot
pF <- result_F$plot

# -----------------------------------------------------------------------------
# 步骤 5: 汇总并输出统计结果
# -----------------------------------------------------------------------------

# 收集所有亚组的统计结果
all_results <- bind_rows(
  tibble(
    Subgroup = result_A$stats$subgroup,
    Variable = result_A$stats$subgroup_var,
    P_Overall = result_A$stats$p_overall,
    P_Nonlinear = result_A$stats$p_nonlinear,
    P_Interaction = result_A$stats$p_interaction
  ),
  tibble(
    Subgroup = result_B$stats$subgroup,
    Variable = result_B$stats$subgroup_var,
    P_Overall = result_B$stats$p_overall,
    P_Nonlinear = result_B$stats$p_nonlinear,
    P_Interaction = result_B$stats$p_interaction
  ),
  tibble(
    Subgroup = result_C$stats$subgroup,
    Variable = result_C$stats$subgroup_var,
    P_Overall = result_C$stats$p_overall,
    P_Nonlinear = result_C$stats$p_nonlinear,
    P_Interaction = result_C$stats$p_interaction
  ),
  tibble(
    Subgroup = result_D$stats$subgroup,
    Variable = result_D$stats$subgroup_var,
    P_Overall = result_D$stats$p_overall,
    P_Nonlinear = result_D$stats$p_nonlinear,
    P_Interaction = result_D$stats$p_interaction
  ),
  tibble(
    Subgroup = result_E$stats$subgroup,
    Variable = result_E$stats$subgroup_var,
    P_Overall = result_E$stats$p_overall,
    P_Nonlinear = result_E$stats$p_nonlinear,
    P_Interaction = result_E$stats$p_interaction
  ),
  tibble(
    Subgroup = result_F$stats$subgroup,
    Variable = result_F$stats$subgroup_var,
    P_Overall = result_F$stats$p_overall,
    P_Nonlinear = result_F$stats$p_nonlinear,
    P_Interaction = result_F$stats$p_interaction
  )
)

# 打印统计结果摘要
cat("\n========================================\n")
cat("限制性立方样条分亚组分析统计结果汇总\n")
cat("========================================\n\n")
print(all_results, n = Inf)

# 保存统计结果为CSV文件
write.csv(all_results, "RCS_Subgroup_Statistics.csv", row.names = FALSE)
cat("\n统计结果已保存为: RCS_Subgroup_Statistics.csv\n")

# 打印详细的模型结果
cat("\n========================================\n")
cat("各亚组详细模型结果\n")
cat("========================================\n\n")

cat("\n--- Subgroup A: Sex ---\n")
print(result_A$stats$model_interaction)

cat("\n--- Subgroup B: Age Group ---\n")
print(result_B$stats$model_interaction)

cat("\n--- Subgroup C: BMI Group ---\n")
print(result_C$stats$model_interaction)

cat("\n--- Subgroup D: Ablation ---\n")
print(result_D$stats$model_interaction)

cat("\n--- Subgroup E: Diabetes ---\n")
print(result_E$stats$model_interaction)

cat("\n--- Subgroup F: LV Hypertrophy ---\n")
print(result_F$stats$model_interaction)

# -----------------------------------------------------------------------------
# 步骤 6: 拼接图形并保存
# -----------------------------------------------------------------------------

# 使用 patchwork 拼接所有图形，并增加垂直和水平间距
final_plot <- (pA | pB) / (pC | pD) / (pE | pF) + 
  plot_layout(heights = c(1, 1, 1)) & 
  theme(plot.margin = margin(t = 20, r = 15, b = 20, l = 15, unit = "pt"))

# 保存图形为PNG文件 - 使用png设备，增加宽度保持合理纵横比
png("RCS_Subgroup_Analysis.png", width = 24*300, height = 28*300, res = 300, bg = "white")
print(final_plot)
dev.off()

# 也保存为PDF格式以便于编辑
ggsave("RCS_Subgroup_Analysis.pdf", final_plot, width = 24, height = 28)

# 打印完成信息
cat("\n========================================\n")
cat("所有文件已成功保存\n")
cat("========================================\n")
cat("图形文件:\n")
cat("- RCS_Subgroup_Analysis.png\n")
cat("- RCS_Subgroup_Analysis.pdf\n")
cat("\n统计结果文件:\n")
cat("- RCS_Subgroup_Statistics.csv\n")
cat("========================================\n\n")

# 显示最终图形
final_plot