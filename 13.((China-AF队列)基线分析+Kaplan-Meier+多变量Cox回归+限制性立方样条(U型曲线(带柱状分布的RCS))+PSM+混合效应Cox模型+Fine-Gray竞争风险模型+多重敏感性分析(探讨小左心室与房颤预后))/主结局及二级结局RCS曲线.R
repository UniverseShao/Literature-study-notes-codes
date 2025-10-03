# =============================================================================
# 主结局及二级结局RCS曲线绘制脚本
# 功能：绘制LVEDD与多个心血管结局的限制性立方样条曲线，并叠加LVEDD分布柱状图
# =============================================================================

# -----------------------------------------------------------------------------
# 步骤 1: 加载必要的R包
# -----------------------------------------------------------------------------
library(tidyverse)
library(survival)
library(rms)       # 用于限制性立方样条Cox回归
library(patchwork) # 用于图形拼接

# -----------------------------------------------------------------------------
# 步骤 2: 创建模拟数据集
# -----------------------------------------------------------------------------
set.seed(2024)
n_patients <- 3000

# 创建模拟的China-AF队列数据
sim_data <- tibble(
  id = 1:n_patients,
  
  # 模拟LVEDD（左心室舒张末期内径）
  # 均值48mm，标准差8mm，符合真实临床分布
  LVEDD = rnorm(n_patients, mean = 48, sd = 8),
  
  # 模拟年龄、性别等协变量
  age = rnorm(n_patients, mean = 65, sd = 10),
  sex = factor(sample(c("Male", "Female"), n_patients, replace = TRUE, prob = c(0.6, 0.4))),
  
  # 模拟U型风险关系的真实对数风险
  # LVEDD偏离48mm时风险增加
  true_log_hazard = 0.005 * (LVEDD - 48)^2,
  
  # ----------- 结局1: 主要复合结局（心血管死亡或大出血）-----------
  # 事件发生率较高，U型关系明显
  time_to_event1 = rexp(n_patients, rate = 0.08 * exp(true_log_hazard)),
  censoring_time1 = runif(n_patients, 1, 10),
  time1 = pmin(time_to_event1, censoring_time1),
  event1 = as.numeric(time_to_event1 <= censoring_time1),
  
  # ----------- 结局2: 心血管死亡 -----------
  # 事件发生率中等，U型关系显著
  time_to_event2 = rexp(n_patients, rate = 0.06 * exp(true_log_hazard * 1.2)),
  censoring_time2 = runif(n_patients, 1, 10),
  time2 = pmin(time_to_event2, censoring_time2),
  event2 = as.numeric(time_to_event2 <= censoring_time2),
  
  # ----------- 结局3: 缺血性卒中/系统性栓塞 -----------
  # 事件发生率较低，关系较弱（P值不显著）
  time_to_event3 = rexp(n_patients, rate = 0.03 * exp(true_log_hazard * 0.3)),
  censoring_time3 = runif(n_patients, 1, 10),
  time3 = pmin(time_to_event3, censoring_time3),
  event3 = as.numeric(time_to_event3 <= censoring_time3),
  
  # ----------- 结局4: 大出血 -----------
  # 事件发生率较低，关系边缘显著
  time_to_event4 = rexp(n_patients, rate = 0.04 * exp(true_log_hazard * 0.5)),
  censoring_time4 = runif(n_patients, 1, 10),
  time4 = pmin(time_to_event4, censoring_time4),
  event4 = as.numeric(time_to_event4 <= censoring_time4)
)

# 设置rms包的数据分布环境
ddist <- datadist(sim_data)
options(datadist = 'ddist')

# -----------------------------------------------------------------------------
# 步骤 3: 创建LVEDD分布数据（用于柱状图）
# -----------------------------------------------------------------------------
# 计算LVEDD的直方图数据
# 动态设置breaks以覆盖所有数据
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
# 步骤 4: 创建绘图函数（RCS曲线 + LVEDD分布柱状图）
# -----------------------------------------------------------------------------

#' @title 创建RCS曲线与LVEDD分布组合图
#' @param data 数据集
#' @param time_var 生存时间变量名（字符串）
#' @param event_var 事件变量名（字符串）
#' @param hist_data LVEDD分布数据框
#' @param plot_title 图的标题（如"A"）
#' @param outcome_label 结局标签
#' @return ggplot对象

create_rcs_histogram_plot <- function(data, time_var, event_var, hist_data, 
                                      plot_title, outcome_label) {
  
  # 构建Cox模型公式
  formula_str <- paste0("Surv(", time_var, ", ", event_var, ") ~ rcs(LVEDD, 4)")
  model_formula <- as.formula(formula_str)
  
  # 拟合包含RCS的Cox模型
  model_rcs <- cph(model_formula, data = data, x = TRUE, y = TRUE)
  
  # 拟合线性模型（用于比较）
  formula_linear <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ LVEDD"))
  model_linear <- cph(formula_linear, data = data, x = TRUE, y = TRUE)
  
  # 提取P值
  anova_result <- anova(model_rcs)
  
  # P overall: LVEDD总体效应的P值
  lvedd_rows <- grep("LVEDD", rownames(anova_result), ignore.case = TRUE)
  p_overall <- if(length(lvedd_rows) > 0) anova_result[lvedd_rows[1], "P"] else NA
  
  # P nonlinear: 非线性效应的P值
  nonlinear_rows <- grep("Nonlinear", rownames(anova_result), ignore.case = TRUE)
  p_nonlinear <- if(length(nonlinear_rows) > 0) anova_result[nonlinear_rows[1], "P"] else NA
  
  # 使用Predict函数生成预测数据
  # ref.zero = TRUE 表示在参考点（中位数或均值）HR=1
  pred_data <- Predict(model_rcs, LVEDD, fun = exp, ref.zero = TRUE) %>% 
    as_tibble() %>%
    rename(HR = yhat, lower_ci = lower, upper_ci = upper)
  
  # 确定Y轴范围
  max_hr <- max(pred_data$upper_ci, na.rm = TRUE)
  y_limit <- min(max(ceiling(max_hr), 3), 8)  # 最大不超过8
  
  # 确定柱状图的缩放比例（使柱状图最高不超过Y轴的1/4）
  max_percentage <- max(hist_data$percentage)
  scale_factor <- y_limit * 0.25 / max_percentage
  
  # 创建注释文本
  annotation_text <- paste0(
    "P overall = ", format.pval(p_overall, digits = 3, eps = 0.001), "\n",
    "P nonlinear = ", format.pval(p_nonlinear, digits = 3, eps = 0.001)
  )
  
  # 创建组合图
  p <- ggplot() +
    # 1. 绘制柱状图（背景层）
    geom_col(data = hist_data, 
             aes(x = midpoint, y = percentage * scale_factor),
             fill = "#6BAED6", alpha = 0.6, width = 1.8) +
    
    # 2. 绘制HR=1的参考线
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.5) +
    
    # 3. 绘制95%置信区间（虚线）
    geom_line(data = pred_data, 
              aes(x = LVEDD, y = lower_ci), 
              linetype = "dashed", color = "#D6604D", size = 0.7) +
    geom_line(data = pred_data, 
              aes(x = LVEDD, y = upper_ci), 
              linetype = "dashed", color = "#D6604D", size = 0.7) +
    
    # 4. 绘制HR估计值曲线（实线）
    geom_line(data = pred_data, 
              aes(x = LVEDD, y = HR), 
              color = "#D6604D", size = 1.2) +
    
    # 5. 添加统计注释
    annotate("text", x = Inf, y = Inf, 
             label = annotation_text, 
             hjust = 1.05, vjust = 1.3, 
             size = 3.5, lineheight = 0.9) +
    
    # 6. 设置坐标轴
    scale_x_continuous(
      name = "LVEDD (mm)",
      breaks = seq(30, 70, by = 10),
      limits = c(30, 75)
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
    
    # 7. 添加标题和主题
    labs(title = plot_title) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10, color = "black"),
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 10, r = 15, b = 10, l = 10, unit = "pt")
    )
  
  # 返回图形和统计结果
  return(list(
    plot = p,
    p_overall = p_overall,
    p_nonlinear = p_nonlinear,
    model = model_rcs
  ))
}

# -----------------------------------------------------------------------------
# 步骤 5: 为每个结局生成RCS曲线图
# -----------------------------------------------------------------------------

cat("\n========================================\n")
cat("开始绘制主结局及二级结局RCS曲线\n")
cat("========================================\n\n")

# Panel A: 主要复合结局（心血管死亡或大出血）
cat("正在分析结局A: 主要复合结局（心血管死亡或大出血）...\n")
result_A <- create_rcs_histogram_plot(
  data = sim_data,
  time_var = "time1",
  event_var = "event1",
  hist_data = hist_data,
  plot_title = "A",
  outcome_label = "Cardiovascular death or major bleeding"
)

# Panel B: 心血管死亡
cat("正在分析结局B: 心血管死亡...\n")
result_B <- create_rcs_histogram_plot(
  data = sim_data,
  time_var = "time2",
  event_var = "event2",
  hist_data = hist_data,
  plot_title = "B",
  outcome_label = "Cardiovascular death"
)

# Panel C: 缺血性卒中/系统性栓塞
cat("正在分析结局C: 缺血性卒中/系统性栓塞...\n")
result_C <- create_rcs_histogram_plot(
  data = sim_data,
  time_var = "time3",
  event_var = "event3",
  hist_data = hist_data,
  plot_title = "C",
  outcome_label = "Ischemic stroke/systemic embolism"
)

# Panel D: 大出血
cat("正在分析结局D: 大出血...\n")
result_D <- create_rcs_histogram_plot(
  data = sim_data,
  time_var = "time4",
  event_var = "event4",
  hist_data = hist_data,
  plot_title = "D",
  outcome_label = "Major bleeding"
)

# -----------------------------------------------------------------------------
# 步骤 6: 汇总统计结果
# -----------------------------------------------------------------------------

cat("\n========================================\n")
cat("统计结果汇总\n")
cat("========================================\n\n")

# 创建统计结果表
results_summary <- tibble(
  Outcome = c(
    "A: CV death or major bleeding",
    "B: Cardiovascular death",
    "C: Ischemic stroke/SE",
    "D: Major bleeding"
  ),
  Events = c(
    sum(sim_data$event1),
    sum(sim_data$event2),
    sum(sim_data$event3),
    sum(sim_data$event4)
  ),
  Event_Rate = sprintf("%.1f%%", 100 * Events / n_patients),
  P_Overall = c(
    result_A$p_overall,
    result_B$p_overall,
    result_C$p_overall,
    result_D$p_overall
  ),
  P_Nonlinear = c(
    result_A$p_nonlinear,
    result_B$p_nonlinear,
    result_C$p_nonlinear,
    result_D$p_nonlinear
  )
)

print(results_summary, n = Inf)

# 保存统计结果
write.csv(results_summary, "主结局及二级结局RCS统计结果.csv", row.names = FALSE)
cat("\n统计结果已保存为: 主结局及二级结局RCS统计结果.csv\n")

# -----------------------------------------------------------------------------
# 步骤 7: 拼接图形并保存
# -----------------------------------------------------------------------------

cat("\n========================================\n")
cat("正在拼接图形并保存...\n")
cat("========================================\n\n")

# 使用patchwork拼接4个图形（2行2列）
final_plot <- (result_A$plot | result_B$plot) / (result_C$plot | result_D$plot) +
  plot_annotation(
    title = "FIGURE 2  Adjusted Association of LVEDD and Cardiovascular Outcomes",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0)
    )
  )

# 保存为PNG（高分辨率）
png("主结局及二级结局RCS曲线.png", 
    width = 3200, height = 2800, res = 400, bg = "white")
print(final_plot)
dev.off()

# 保存为PDF（矢量图，可编辑）
ggsave("主结局及二级结局RCS曲线.pdf", 
       final_plot, 
       width = 8, height = 7, 
       device = "pdf")

# 保存为高分辨率TIFF（适合投稿）
ggsave("主结局及二级结局RCS曲线.tiff", 
       final_plot, 
       width = 8, height = 7, 
       dpi = 600, 
       device = "tiff", 
       compression = "lzw")

cat("\n========================================\n")
cat("所有文件已成功保存\n")
cat("========================================\n")
cat("图形文件:\n")
cat("- 主结局及二级结局RCS曲线.png (高分辨率，400 DPI)\n")
cat("- 主结局及二级结局RCS曲线.pdf (矢量图，可编辑)\n")
cat("- 主结局及二级结局RCS曲线.tiff (投稿用，600 DPI)\n")
cat("\n统计结果文件:\n")
cat("- 主结局及二级结局RCS统计结果.csv\n")
cat("========================================\n\n")

# 打印详细的模型结果
cat("\n========================================\n")
cat("详细模型结果\n")
cat("========================================\n\n")

cat("\n--- 结局A: 主要复合结局 ---\n")
print(result_A$model)

cat("\n--- 结局B: 心血管死亡 ---\n")
print(result_B$model)

cat("\n--- 结局C: 缺血性卒中/系统性栓塞 ---\n")
print(result_C$model)

cat("\n--- 结局D: 大出血 ---\n")
print(result_D$model)

cat("\n========================================\n")
cat("分析完成！\n")
cat("========================================\n\n")

# 显示最终图形
print(final_plot)

