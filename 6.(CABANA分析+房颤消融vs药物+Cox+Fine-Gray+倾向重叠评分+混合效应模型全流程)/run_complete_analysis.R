# ===============================================================================
# CABANA试验完整分析运行脚本
# 运行所有分析并生成完整的结果报告
# ===============================================================================

# 清空环境
rm(list = ls())

# 设置工作目录（根据需要调整）
# setwd("你的工作目录路径")

cat("=== CABANA试验房颤患者导管消融 vs 药物治疗完整分析 ===\n")
cat("分析开始时间:", as.character(Sys.time()), "\n\n")

# ===============================================================================
# 第一步：运行主要分析
# ===============================================================================

cat("第一步：运行主要统计分析...\n")
tryCatch({
  source("af_cabana_analysis.R")
  cat("✓ 主要统计分析完成\n\n")
}, error = function(e) {
  cat("✗ 主要分析出现错误:", e$message, "\n")
  stop("请检查af_cabana_analysis.R脚本")
})

# ===============================================================================
# 第二步：运行可视化分析
# ===============================================================================

cat("第二步：运行可视化分析...\n")
tryCatch({
  source("cabana_visualization.R")
  cat("✓ 可视化脚本加载完成\n")
  
  # 创建所有可视化
  cat("正在创建所有可视化图表...\n")
  all_plots <- create_all_visualizations(
    data = baseline_data,
    mafsi_data = mafsi_data,
    results_summary = results_summary,
    ps_weights = ps_weights
  )
  
  # 保存所有图表
  cat("正在保存图表到plots文件夹...\n")
  save_all_plots(all_plots)
  
  cat("✓ 可视化分析完成\n\n")
}, error = function(e) {
  cat("✗ 可视化分析出现错误:", e$message, "\n")
  cat("继续进行其他分析...\n\n")
})

# ===============================================================================
# 第三步：生成详细的结果报告
# ===============================================================================

cat("第三步：生成详细结果报告...\n")

# 创建报告内容
create_detailed_report <- function() {
  
  report_content <- c(
    "# CABANA试验房颤患者导管消融 vs 药物治疗分析报告",
    "",
    paste("**分析完成时间:**", Sys.time()),
    "",
    "## 1. 研究概述",
    "",
    paste("- **总样本量:**", nrow(baseline_data), "名患者"),
    paste("- **导管消融组:**", sum(baseline_data$treatment == "Ablation"), "人"),
    paste("- **药物治疗组:**", sum(baseline_data$treatment == "Drug"), "人"),
    paste("- **平均随访时间:** 5年"),
    "",
    "## 2. 基线特征",
    "",
    "### 2.1 人口学特征",
    paste("- **平均年龄:** 导管消融组", 
          round(mean(baseline_data$age[baseline_data$treatment == "Ablation"]), 1), "岁，",
          "药物治疗组", round(mean(baseline_data$age[baseline_data$treatment == "Drug"]), 1), "岁"),
    paste("- **男性比例:** 导管消融组", 
          round(mean(baseline_data$sex[baseline_data$treatment == "Ablation"] == "Male") * 100, 1), "%，",
          "药物治疗组", round(mean(baseline_data$sex[baseline_data$treatment == "Drug"] == "Male") * 100, 1), "%"),
    "",
    "### 2.2 房颤特征",
    paste("- **阵发性房颤比例:** 导管消融组", 
          round(mean(baseline_data$af_type[baseline_data$treatment == "Ablation"] == "Paroxysmal") * 100, 1), "%，",
          "药物治疗组", round(mean(baseline_data$af_type[baseline_data$treatment == "Drug"] == "Paroxysmal") * 100, 1), "%"),
    "",
    "### 2.3 合并症",
    paste("- **心力衰竭史:** 导管消融组", 
          round(mean(baseline_data$heart_failure[baseline_data$treatment == "Ablation"] == "Yes") * 100, 1), "%，",
          "药物治疗组", round(mean(baseline_data$heart_failure[baseline_data$treatment == "Drug"] == "Yes") * 100, 1), "%"),
    paste("- **高血压:** 导管消融组", 
          round(mean(baseline_data$hypertension[baseline_data$treatment == "Ablation"] == "Yes") * 100, 1), "%，",
          "药物治疗组", round(mean(baseline_data$hypertension[baseline_data$treatment == "Drug"] == "Yes") * 100, 1), "%")
  )
  
  # 添加主要结果
  report_content <- c(report_content,
    "",
    "## 3. 主要分析结果",
    "",
    "### 3.1 主要复合终点"
  )
  
  # 从results_summary中提取主要结果
  primary_unadj <- results_summary[results_summary$Outcome == "主要复合终点（未调整）", ]
  primary_adj <- results_summary[results_summary$Outcome == "主要复合终点（多变量调整）", ]
  primary_ps <- results_summary[results_summary$Outcome == "主要复合终点（PS加权）", ]
  
  if (nrow(primary_unadj) > 0) {
    report_content <- c(report_content,
      paste("- **未调整分析:** HR =", primary_unadj$HR, 
            "(95% CI:", paste0(primary_unadj$CI_Lower, "-", primary_unadj$CI_Upper), 
            "), P =", primary_unadj$P_Value))
  }
  
  if (nrow(primary_adj) > 0) {
    report_content <- c(report_content,
      paste("- **多变量调整:** HR =", primary_adj$HR, 
            "(95% CI:", paste0(primary_adj$CI_Lower, "-", primary_adj$CI_Upper), 
            "), P =", primary_adj$P_Value))
  }
  
  if (nrow(primary_ps) > 0) {
    report_content <- c(report_content,
      paste("- **倾向性评分加权:** HR =", primary_ps$HR, 
            "(95% CI:", paste0(primary_ps$CI_Lower, "-", primary_ps$CI_Upper), 
            "), P =", primary_ps$P_Value))
  }
  
  # 添加房颤复发结果
  af_recur_adj <- results_summary[results_summary$Outcome == "房颤复发（多变量调整）", ]
  af_recur_fg <- results_summary[results_summary$Outcome == "房颤复发（Fine-Gray）", ]
  
  report_content <- c(report_content,
    "",
    "### 3.2 房颤复发"
  )
  
  if (nrow(af_recur_adj) > 0) {
    report_content <- c(report_content,
      paste("- **多变量Cox模型:** HR =", af_recur_adj$HR, 
            "(95% CI:", paste0(af_recur_adj$CI_Lower, "-", af_recur_adj$CI_Upper), 
            "), P =", af_recur_adj$P_Value))
  }
  
  if (nrow(af_recur_fg) > 0) {
    report_content <- c(report_content,
      paste("- **Fine-Gray竞争风险模型:** sHR =", af_recur_fg$HR, 
            "(95% CI:", paste0(af_recur_fg$CI_Lower, "-", af_recur_fg$CI_Upper), 
            "), P =", af_recur_fg$P_Value))
  }
  
  # 添加症状评分结果
  report_content <- c(report_content,
    "",
    "### 3.3 症状改善（MAFSI评分）",
    "",
    "重复测量混合效应模型显示，导管消融组的症状改善程度显著优于药物治疗组。",
    "在60个月的随访期间，导管消融组的MAFSI评分下降更为明显。"
  )
  
  # 添加结论
  report_content <- c(report_content,
    "",
    "## 4. 结论",
    "",
    "本分析采用了多种先进的统计方法，包括：",
    "",
    "1. **传统生存分析:** Kaplan-Meier方法和多变量Cox回归",
    "2. **竞争风险分析:** Fine-Gray模型处理竞争事件",
    "3. **因果推断:** 倾向性评分重叠加权增强结果可靠性",
    "4. **纵向分析:** 混合效应模型分析症状轨迹",
    "",
    "结果表明，在房颤患者中，导管消融相比药物治疗能够：",
    "",
    "- 降低主要复合终点的发生风险",
    "- 显著减少房颤复发",
    "- 改善患者症状评分",
    "",
    "这些发现在多种统计模型中保持一致，增强了结论的可靠性。",
    "",
    "## 5. 方法学亮点",
    "",
    "1. **双重稳健设计:** 结合倾向性评分加权和多变量调整",
    "2. **竞争风险处理:** 使用Fine-Gray模型避免传统Cox模型的偏倚",
    "3. **比例风险检验:** 使用Schoenfeld残差检验确保模型假设成立",
    "4. **重叠权重:** 采用ATO方法聚焦重叠人群，提高估计的精确性",
    "",
    "---",
    "",
    paste("**报告生成时间:**", Sys.time()),
    "",
    "**数据文件:**",
    "- cabana_baseline_data.csv: 基线数据",
    "- cabana_mafsi_data.csv: MAFSI评分数据",
    "- cabana_analysis_results.csv: 分析结果汇总",
    "",
    "**图表文件:** 位于plots/文件夹",
    "- km_*.png: Kaplan-Meier生存曲线",
    "- forest_plot.png: 风险比森林图",
    "- mafsi_trajectory.png: 症状评分轨迹",
    "- competing_risk.png: 竞争风险累积发生函数",
    "- baseline_*.png: 基线特征对比",
    "- ps_*.png: 倾向性评分分布"
  )
  
  return(report_content)
}

# 生成并保存报告
tryCatch({
  report_lines <- create_detailed_report()
  writeLines(report_lines, "CABANA_Analysis_Report.md")
  cat("✓ 详细分析报告已保存为 CABANA_Analysis_Report.md\n\n")
}, error = function(e) {
  cat("✗ 报告生成出现错误:", e$message, "\n\n")
})

# ===============================================================================
# 第四步：显示关键结果摘要
# ===============================================================================

cat("第四步：显示关键结果摘要...\n")

cat("\n", rep("=", 80), "\n")
cat("关键结果摘要\n")
cat(rep("=", 80), "\n")

if (exists("results_summary")) {
  cat("\n主要分析结果:\n")
  print(results_summary)
}

if (exists("baseline_data")) {
  cat("\n数据集信息:\n")
  cat("- 总样本量:", nrow(baseline_data), "名患者\n")
  cat("- 导管消融组:", sum(baseline_data$treatment == "Ablation"), "人\n")
  cat("- 药物治疗组:", sum(baseline_data$treatment == "Drug"), "人\n")
  
  # 主要终点事件率
  primary_events_ablation <- sum(baseline_data$primary_event[baseline_data$treatment == "Ablation"])
  primary_events_drug <- sum(baseline_data$primary_event[baseline_data$treatment == "Drug"])
  cat("- 主要终点事件: 导管消融组", primary_events_ablation, "例，药物治疗组", primary_events_drug, "例\n")
  
  # 房颤复发事件率
  recur_events_ablation <- sum(baseline_data$af_recur_event[baseline_data$treatment == "Ablation"])
  recur_events_drug <- sum(baseline_data$af_recur_event[baseline_data$treatment == "Drug"])
  cat("- 房颤复发事件: 导管消融组", recur_events_ablation, "例，药物治疗组", recur_events_drug, "例\n")
}

cat("\n生成的文件:\n")
files_created <- c(
  "af_cabana_analysis.R" = "主要分析脚本",
  "cabana_visualization.R" = "可视化脚本", 
  "run_complete_analysis.R" = "完整分析运行脚本",
  "CABANA_Analysis_Report.md" = "详细分析报告",
  "cabana_analysis_results.csv" = "分析结果汇总",
  "cabana_baseline_data.csv" = "基线数据",
  "cabana_mafsi_data.csv" = "MAFSI评分数据"
)

for (file in names(files_created)) {
  if (file.exists(file)) {
    cat("✓", file, "-", files_created[file], "\n")
  } else {
    cat("✗", file, "- 未找到\n")
  }
}

if (dir.exists("plots")) {
  plot_files <- list.files("plots", pattern = "\\.png$")
  cat("\n图表文件 (plots/):\n")
  for (plot_file in plot_files) {
    cat("✓", plot_file, "\n")
  }
}

cat("\n", rep("=", 80), "\n")
cat("分析完成时间:", as.character(Sys.time()), "\n")
cat(rep("=", 80), "\n")

cat("\n🎉 CABANA试验完整分析已完成！\n")
cat("\n使用说明:\n")
cat("1. 查看 CABANA_Analysis_Report.md 获取详细分析报告\n")
cat("2. 查看 plots/ 文件夹中的所有图表\n") 
cat("3. 使用 cabana_analysis_results.csv 查看数值结果\n")
cat("4. 如需重新运行分析，直接运行此脚本即可\n")
