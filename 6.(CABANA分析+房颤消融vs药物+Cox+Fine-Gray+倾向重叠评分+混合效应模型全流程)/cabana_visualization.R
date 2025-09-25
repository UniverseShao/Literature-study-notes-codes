# ===============================================================================
# CABANA试验可视化脚本 - 生成高质量的分析图表
# 配合主分析脚本af_cabana_analysis.R使用
# ===============================================================================

# 加载必要的包
library(tidyverse)
library(survival)
library(survminer)
library(cmprsk)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(scales)
library(forestplot)
library(patchwork)

# 设置主题
theme_set(theme_minimal() + 
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  ))

# 定义配色方案
colors_treatment <- c("Drug" = "#E74C3C", "Ablation" = "#3498DB")

# ===============================================================================
# 函数定义：创建各种图表
# ===============================================================================

# 1. 基线特征对比图
create_baseline_comparison <- function(data) {
  
  # 连续变量的分布对比
  continuous_vars <- c("age", "af_duration", "cha2ds2_vasc")
  
  plot_list <- list()
  
  for (var in continuous_vars) {
    var_label <- switch(var,
                       "age" = "年龄 (岁)",
                       "af_duration" = "房颤病程 (年)", 
                       "cha2ds2_vasc" = "CHA₂DS₂-VASc评分")
    
    p <- ggplot(data, aes(x = !!sym(var), fill = treatment)) +
      geom_density(alpha = 0.7) +
      scale_fill_manual(values = colors_treatment, 
                       name = "治疗组", 
                       labels = c("药物治疗", "导管消融")) +
      labs(title = paste0(var_label, "分布对比"),
           x = var_label,
           y = "密度") +
      theme(legend.position = "bottom")
    
    plot_list[[var]] <- p
  }
  
  # 分类变量的条形图
  categorical_vars <- c("sex", "af_type", "heart_failure")
  
  for (var in categorical_vars) {
    var_label <- switch(var,
                       "sex" = "性别",
                       "af_type" = "房颤类型",
                       "heart_failure" = "心力衰竭史")
    
    prop_data <- data %>%
      group_by(treatment, !!sym(var)) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(treatment) %>%
      mutate(prop = n / sum(n) * 100)
    
    p <- ggplot(prop_data, aes(x = !!sym(var), y = prop, fill = treatment)) +
      geom_col(position = "dodge", alpha = 0.8) +
      scale_fill_manual(values = colors_treatment,
                       name = "治疗组",
                       labels = c("药物治疗", "导管消融")) +
      labs(title = paste0(var_label, "分布对比"),
           x = var_label,
           y = "百分比 (%)") +
      theme(legend.position = "bottom")
    
    plot_list[[paste0(var, "_cat")]] <- p
  }
  
  return(plot_list)
}

# 2. Kaplan-Meier生存曲线
create_km_curves <- function(data) {
  
  plot_list <- list()
  
  # 主要复合终点
  fit_primary <- survfit(Surv(primary_time/365.25, primary_event) ~ treatment, data = data)
  
  p1 <- ggsurvplot(
    fit_primary,
    data = data,
    title = "主要复合终点的生存曲线",
    xlab = "时间 (年)",
    ylab = "无事件生存率",
    legend.title = "治疗组",
    legend.labs = c("药物治疗", "导管消融"),
    palette = c("#E74C3C", "#3498DB"),
    risk.table = TRUE,
    risk.table.height = 0.3,
    pval = TRUE,
    conf.int = TRUE,
    ggtheme = theme_minimal()
  )
  
  plot_list[["primary"]] <- p1
  
  # 房颤复发
  fit_recurrence <- survfit(Surv(af_recur_time/365.25, af_recur_event) ~ treatment, data = data)
  
  p2 <- ggsurvplot(
    fit_recurrence,
    data = data,
    title = "房颤复发的生存曲线",
    xlab = "时间 (年)",
    ylab = "无复发生存率",
    legend.title = "治疗组",
    legend.labs = c("药物治疗", "导管消融"),
    palette = c("#E74C3C", "#3498DB"),
    risk.table = TRUE,
    risk.table.height = 0.3,
    pval = TRUE,
    conf.int = TRUE,
    ggtheme = theme_minimal()
  )
  
  plot_list[["recurrence"]] <- p2
  
  # 全因死亡
  fit_death <- survfit(Surv(death_time/365.25, death_event) ~ treatment, data = data)
  
  p3 <- ggsurvplot(
    fit_death,
    data = data,
    title = "全因死亡的生存曲线",
    xlab = "时间 (年)",
    ylab = "生存率",
    legend.title = "治疗组",
    legend.labs = c("药物治疗", "导管消融"),
    palette = c("#E74C3C", "#3498DB"),
    risk.table = TRUE,
    risk.table.height = 0.3,
    pval = TRUE,
    conf.int = TRUE,
    ggtheme = theme_minimal()
  )
  
  plot_list[["death"]] <- p3
  
  return(plot_list)
}

# 3. 竞争风险累积发生函数图
create_competing_risk_plot <- function(data) {
  
  # 准备竞争风险数据
  comp_risk_data <- data %>%
    mutate(
      time_years = af_recur_comp_time / 365.25,
      status = af_recur_comp_event
    )
  
  # 计算累积发生函数
  cif_result <- cuminc(ftime = comp_risk_data$time_years,
                       fstatus = comp_risk_data$status,
                       group = comp_risk_data$treatment)
  
  # 转换为ggplot可用的数据格式
  cif_df <- data.frame()
  
  for (group in names(cif_result)) {
    if (group %in% c("Drug 1", "Ablation 1", "Drug 2", "Ablation 2")) {
      treatment_group <- strsplit(group, " ")[[1]][1]
      event_type <- ifelse(strsplit(group, " ")[[1]][2] == "1", "房颤复发", "死亡")
      
      temp_df <- data.frame(
        time = cif_result[[group]]$time,
        est = cif_result[[group]]$est,
        treatment = treatment_group,
        event = event_type
      )
      
      cif_df <- rbind(cif_df, temp_df)
    }
  }
  
  # 创建图表
  p <- ggplot(cif_df, aes(x = time, y = est, color = treatment, linetype = event)) +
    geom_step(size = 1) +
    scale_color_manual(values = colors_treatment,
                      name = "治疗组",
                      labels = c("导管消融", "药物治疗")) +
    scale_linetype_manual(values = c("房颤复发" = "solid", "死亡" = "dashed"),
                         name = "事件类型") +
    labs(title = "竞争风险累积发生函数",
         subtitle = "房颤复发 vs 死亡",
         x = "时间 (年)",
         y = "累积发生率") +
    scale_y_continuous(labels = percent_format()) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# 4. Forest Plot - 风险比可视化
create_forest_plot <- function(results_data) {
  
  # 准备数据
  forest_data <- results_data %>%
    filter(!grepl("未调整", Outcome)) %>%
    mutate(
      outcome_clean = case_when(
        grepl("主要复合终点.*多变量", Outcome) ~ "主要复合终点 (多变量Cox)",
        grepl("主要复合终点.*PS", Outcome) ~ "主要复合终点 (PS加权Cox)",
        grepl("房颤复发.*多变量", Outcome) ~ "房颤复发 (多变量Cox)",
        grepl("房颤复发.*PS加权Cox", Outcome) ~ "房颤复发 (PS加权Cox)",
        grepl("房颤复发.*Fine-Gray", Outcome) ~ "房颤复发 (Fine-Gray)",
        grepl("房颤复发.*PS加权Fine", Outcome) ~ "房颤复发 (PS加权Fine-Gray)",
        grepl("全因死亡", Outcome) ~ "全因死亡 (多变量Cox)",
        TRUE ~ Outcome
      ),
      ci_text = paste0(HR, " (", CI_Lower, "-", CI_Upper, ")"),
      significant = ifelse(P_Value < 0.05, "显著", "不显著")
    ) %>%
    arrange(desc(row_number()))
  
  # 创建forest plot
  p <- ggplot(forest_data, aes(y = reorder(outcome_clean, row_number()))) +
    geom_point(aes(x = HR, color = significant), size = 3) +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper, color = significant), 
                   height = 0.2, size = 1) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
    scale_color_manual(values = c("显著" = "#E74C3C", "不显著" = "#95A5A6"),
                      name = "统计学意义") +
    scale_x_log10(breaks = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 2),
                  labels = c("0.25", "0.50", "0.75", "1.00", "1.25", "1.50", "2.00")) +
    labs(title = "导管消融 vs 药物治疗的风险比",
         subtitle = "风险比 < 1 表示导管消融降低风险",
         x = "风险比 (HR/sHR)",
         y = "") +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    ) +
    annotate("text", x = 0.3, y = nrow(forest_data) + 0.5, 
             label = "导管消融更优", size = 3, color = "#2ECC71") +
    annotate("text", x = 1.8, y = nrow(forest_data) + 0.5, 
             label = "药物治疗更优", size = 3, color = "#E74C3C")
  
  return(p)
}

# 5. MAFSI症状评分轨迹图
create_mafsi_trajectory <- function(mafsi_data) {
  
  # 计算每个时间点的均值和标准误
  mafsi_summary <- mafsi_data %>%
    group_by(time_months, treatment) %>%
    summarise(
      mean_score = mean(mafsi_score, na.rm = TRUE),
      se_score = sd(mafsi_score, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      ci_lower = mean_score - 1.96 * se_score,
      ci_upper = mean_score + 1.96 * se_score
    )
  
  # 创建轨迹图
  p <- ggplot(mafsi_summary, aes(x = time_months, y = mean_score, color = treatment)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = treatment), 
                alpha = 0.2, color = NA) +
    scale_color_manual(values = colors_treatment,
                      name = "治疗组",
                      labels = c("导管消融", "药物治疗")) +
    scale_fill_manual(values = colors_treatment,
                     name = "治疗组",
                     labels = c("导管消融", "药物治疗")) +
    labs(title = "MAFSI症状评分随时间变化轨迹",
         subtitle = "评分越低表示症状越轻",
         x = "随访时间 (月)",
         y = "MAFSI评分 (均值 ± 95%CI)") +
    scale_x_continuous(breaks = c(0, 6, 12, 24, 36, 48, 60)) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# 6. 倾向性评分分布图
create_ps_distribution <- function(data, ps_weights) {
  
  # 获取倾向性评分
  ps_data <- data.frame(
    treatment = data$treatment,
    ps_score = ps_weights$ps,
    weights = ps_weights$weights
  )
  
  # 倾向性评分分布对比
  p1 <- ggplot(ps_data, aes(x = ps_score, fill = treatment)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = colors_treatment,
                     name = "治疗组",
                     labels = c("药物治疗", "导管消融")) +
    labs(title = "倾向性评分分布",
         x = "倾向性评分",
         y = "密度") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # 权重分布
  p2 <- ggplot(ps_data, aes(x = weights, fill = treatment)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = colors_treatment,
                     name = "治疗组",
                     labels = c("药物治疗", "导管消融")) +
    labs(title = "重叠权重分布",
         x = "权重",
         y = "密度") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # 结合图
  combined_plot <- p1 / p2
  
  return(list(ps_dist = p1, weight_dist = p2, combined = combined_plot))
}

# ===============================================================================
# 主要可视化函数
# ===============================================================================

create_all_visualizations <- function(data, mafsi_data, results_summary, ps_weights) {
  
  cat("开始创建可视化图表...\n")
  
  # 1. 基线特征对比
  cat("创建基线特征对比图...\n")
  baseline_plots <- create_baseline_comparison(data)
  
  # 2. Kaplan-Meier曲线
  cat("创建Kaplan-Meier生存曲线...\n")
  km_plots <- create_km_curves(data)
  
  # 3. 竞争风险图
  cat("创建竞争风险累积发生函数图...\n")
  competing_risk_plot <- create_competing_risk_plot(data)
  
  # 4. Forest plot
  cat("创建Forest plot...\n")
  forest_plot <- create_forest_plot(results_summary)
  
  # 5. MAFSI轨迹图
  cat("创建MAFSI症状评分轨迹图...\n")
  mafsi_plot <- create_mafsi_trajectory(mafsi_data)
  
  # 6. 倾向性评分分布
  cat("创建倾向性评分分布图...\n")
  ps_plots <- create_ps_distribution(data, ps_weights)
  
  # 保存所有图表
  plot_list <- list(
    baseline = baseline_plots,
    kaplan_meier = km_plots,
    competing_risk = competing_risk_plot,
    forest = forest_plot,
    mafsi = mafsi_plot,
    propensity_score = ps_plots
  )
  
  return(plot_list)
}

# ===============================================================================
# 保存图表函数
# ===============================================================================

save_all_plots <- function(plot_list, output_dir = "plots") {
  
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat(sprintf("保存图表到目录: %s\n", output_dir))
  
  # 保存基线特征对比图
  for (i in seq_along(plot_list$baseline)) {
    plot_name <- names(plot_list$baseline)[i]
    filename <- file.path(output_dir, paste0("baseline_", plot_name, ".png"))
    ggsave(filename, plot_list$baseline[[i]], width = 10, height = 6, dpi = 300)
  }
  
  # 保存Kaplan-Meier曲线
  for (i in seq_along(plot_list$kaplan_meier)) {
    outcome_name <- names(plot_list$kaplan_meier)[i]
    filename <- file.path(output_dir, paste0("km_", outcome_name, ".png"))
    ggsave(filename, plot_list$kaplan_meier[[i]]$plot, width = 12, height = 8, dpi = 300)
  }
  
  # 保存竞争风险图
  ggsave(file.path(output_dir, "competing_risk.png"), 
         plot_list$competing_risk, width = 10, height = 6, dpi = 300)
  
  # 保存Forest plot
  ggsave(file.path(output_dir, "forest_plot.png"), 
         plot_list$forest, width = 12, height = 8, dpi = 300)
  
  # 保存MAFSI轨迹图
  ggsave(file.path(output_dir, "mafsi_trajectory.png"), 
         plot_list$mafsi, width = 10, height = 6, dpi = 300)
  
  # 保存倾向性评分图
  ggsave(file.path(output_dir, "ps_distribution.png"), 
         plot_list$propensity_score$ps_dist, width = 10, height = 6, dpi = 300)
  ggsave(file.path(output_dir, "weight_distribution.png"), 
         plot_list$propensity_score$weight_dist, width = 10, height = 6, dpi = 300)
  ggsave(file.path(output_dir, "ps_combined.png"), 
         plot_list$propensity_score$combined, width = 10, height = 10, dpi = 300)
  
  cat("所有图表保存完成！\n")
}

# ===============================================================================
# 使用示例（如果直接运行此脚本）
# ===============================================================================

if (FALSE) {  # 设置为TRUE来运行示例
  
  # 首先运行主分析脚本
  source("af_cabana_analysis.R")
  
  # 创建所有可视化
  all_plots <- create_all_visualizations(
    data = baseline_data,
    mafsi_data = mafsi_data,
    results_summary = results_summary,
    ps_weights = ps_weights
  )
  
  # 保存所有图表
  save_all_plots(all_plots)
  
  # 显示部分图表（可选）
  print(all_plots$forest)
  print(all_plots$mafsi)
}

cat("\n可视化脚本加载完成！\n")
cat("使用方法：\n")
cat("1. 先运行主分析脚本: source('af_cabana_analysis.R')\n")
cat("2. 创建可视化: all_plots <- create_all_visualizations(baseline_data, mafsi_data, results_summary, ps_weights)\n")
cat("3. 保存图表: save_all_plots(all_plots)\n")
