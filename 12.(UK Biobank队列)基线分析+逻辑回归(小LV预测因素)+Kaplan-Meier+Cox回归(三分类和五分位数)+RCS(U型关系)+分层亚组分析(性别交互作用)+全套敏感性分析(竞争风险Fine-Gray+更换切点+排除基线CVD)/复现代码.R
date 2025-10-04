################################################################################
# UK Biobank 左心室大小与心血管结局关联研究 - 完整复现代码
# 
# 研究设计：基线分析 + 逻辑回归 + Cox回归 + RCS + 分层亚组 + 敏感性分析
# 
# 作者：根据研究思路自动生成
# 日期：2025-10-02
################################################################################

# =============================================================================
# 第一部分：环境准备与数据加载
# =============================================================================

# 1.1 清空环境
rm(list = ls())
gc()

# 1.2 加载必需的R包
required_packages <- c(
  "tidyverse",      # 数据处理和可视化
  "survival",       # 生存分析（Cox回归、KM曲线）
  "survminer",      # 生存曲线可视化
  "rms",            # 限制性立方样条（RCS）
  "cmprsk",         # 竞争风险分析
  "forestplot",     # 森林图
  "tableone",       # 基线特征表
  "gtsummary",      # 表格美化
  "ggpubr",         # 图形组合
  "grid",           # 图形布局
  "gridExtra"       # 图形排列
)

# 安装并加载包
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 1.3 设置工作目录（请根据实际情况修改）
# setwd("your_working_directory")

# 1.4 加载数据（假设数据文件名为ukb_data.csv或ukb_data.RData）
# data <- read.csv("ukb_data.csv")
# 或者
# load("ukb_data.RData")

# 示例：创建模拟数据结构（实际使用时请替换为真实数据）
# data <- read_data_function()


# =============================================================================
# 第二部分：数据预处理
# =============================================================================

# 2.1 创建左心室大小分类变量
# 根据LVEDVi（左心室舒张末期容积指数）的标准切点分组
# 参考指南标准：男性 < 61 mL/m² 为小LV，> 96 mL/m² 为大LV
#               女性 < 53 mL/m² 为小LV，> 84 mL/m² 为大LV

data <- data %>%
  mutate(
    LV_size_category = case_when(
      sex == "Male" & LVEDVi < 61 ~ "Small",
      sex == "Male" & LVEDVi >= 61 & LVEDVi <= 96 ~ "Normal",
      sex == "Male" & LVEDVi > 96 ~ "Large",
      sex == "Female" & LVEDVi < 53 ~ "Small",
      sex == "Female" & LVEDVi >= 53 & LVEDVi <= 84 ~ "Normal",
      sex == "Female" & LVEDVi > 84 ~ "Large"
    ),
    LV_size_category = factor(LV_size_category, 
                              levels = c("Small", "Normal", "Large"))
  )

# 2.2 创建二分类变量用于逻辑回归（小LV vs 非小LV）
data <- data %>%
  mutate(small_LV = ifelse(LV_size_category == "Small", 1, 0))

# 2.3 创建LVEDVi五分位数分组
data <- data %>%
  mutate(LVEDVi_quintile = ntile(LVEDVi, 5),
         LVEDVi_quintile = factor(LVEDVi_quintile))

# 2.4 确保结局变量和时间变量正确设置
# 假设结局变量命名为：MACE_event, CHD_event, HF_event, stroke_event, death_event
# 时间变量命名为：MACE_time, CHD_time, HF_time, stroke_time, death_time


# =============================================================================
# 阶段一：基线队列描述与"小左心室"预测因素探索
# =============================================================================

cat("\n========== 阶段一：基线特征分析 ==========\n")

# -----------------------------------------------------------------------------
# 3.1 描述性统计（Table 1）
# -----------------------------------------------------------------------------

# 定义需要展示的变量
baseline_vars <- c(
  "age", "sex", "ethnicity", "education",
  "BMI", "SBP", "DBP", "total_cholesterol",
  "smoking_status", "alcohol_status",
  "hypertension", "diabetes", "prior_CVD",
  "HbA1c", "CRP", "hemoglobin", "albumin",
  "atrial_fibrillation", "LVEF",
  "statin_use", "antihypertensive_use"
)

# 创建分组基线特征表
table1 <- CreateTableOne(
  vars = baseline_vars,
  strata = "LV_size_category",
  data = data,
  test = TRUE,
  addOverall = TRUE
)

# 打印表格
print(table1, showAllLevels = TRUE, cramVars = "sex")

# 保存为CSV
table1_output <- print(table1, printToggle = FALSE, noSpaces = TRUE)
write.csv(table1_output, "Table1_Baseline_Characteristics.csv")

# 使用gtsummary生成更美观的表格
table1_gt <- data %>%
  select(LV_size_category, all_of(baseline_vars)) %>%
  tbl_summary(
    by = LV_size_category,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no"
  ) %>%
  add_p(
    test = list(
      all_continuous() ~ "kruskal.test",
      all_categorical() ~ "chisq.test"
    )
  ) %>%
  add_overall() %>%
  modify_header(label = "**Variable**") %>%
  bold_labels()

# 保存表格
table1_gt %>%
  as_gt() %>%
  gt::gtsave("Table1_Baseline_Characteristics.html")


# -----------------------------------------------------------------------------
# 3.2 逻辑回归分析：预测"小左心室"的因素（Table S1）
# -----------------------------------------------------------------------------

cat("\n----- 3.2 逻辑回归分析：小LV预测因素 -----\n")

# 定义预测变量
predictor_vars <- c(
  "age", "sex", "ethnicity", "BMI", "SBP", "DBP",
  "total_cholesterol", "smoking_status", "alcohol_status",
  "hypertension", "diabetes", "HbA1c", "CRP", "LVEF"
)

# 3.2.1 单变量逻辑回归
univariate_logistic_results <- data.frame()

for (var in predictor_vars) {
  formula_str <- paste0("small_LV ~ ", var)
  model <- glm(as.formula(formula_str), 
               data = data, 
               family = binomial(link = "logit"))
  
  coef_summary <- summary(model)$coefficients
  
  # 提取第二行（第一个预测变量的系数）
  if (nrow(coef_summary) > 1) {
    OR <- exp(coef_summary[2, "Estimate"])
    CI_lower <- exp(coef_summary[2, "Estimate"] - 1.96 * coef_summary[2, "Std. Error"])
    CI_upper <- exp(coef_summary[2, "Estimate"] + 1.96 * coef_summary[2, "Std. Error"])
    p_value <- coef_summary[2, "Pr(>|z|)"]
    
    result_row <- data.frame(
      Variable = var,
      OR = OR,
      CI_95_lower = CI_lower,
      CI_95_upper = CI_upper,
      P_value = p_value
    )
    
    univariate_logistic_results <- rbind(univariate_logistic_results, result_row)
  }
}

print(univariate_logistic_results)
write.csv(univariate_logistic_results, "Univariate_Logistic_Regression.csv", row.names = FALSE)

# 3.2.2 多变量逻辑回归
multivariate_formula <- paste("small_LV ~", paste(predictor_vars, collapse = " + "))
multivariate_logistic <- glm(
  as.formula(multivariate_formula),
  data = data,
  family = binomial(link = "logit")
)

# 提取结果
multivariate_results <- summary(multivariate_logistic)$coefficients
multivariate_results_df <- data.frame(
  Variable = rownames(multivariate_results)[-1],
  OR = exp(multivariate_results[-1, "Estimate"]),
  CI_95_lower = exp(multivariate_results[-1, "Estimate"] - 1.96 * multivariate_results[-1, "Std. Error"]),
  CI_95_upper = exp(multivariate_results[-1, "Estimate"] + 1.96 * multivariate_results[-1, "Std. Error"]),
  P_value = multivariate_results[-1, "Pr(>|z|)"]
)

print(multivariate_results_df)
write.csv(multivariate_results_df, "Multivariate_Logistic_Regression.csv", row.names = FALSE)


# =============================================================================
# 阶段二：核心假设检验 - LV大小与心血管结局的关联
# =============================================================================

cat("\n========== 阶段二：Cox回归分析 ==========\n")

# 定义所有结局
outcomes <- c("MACE", "CHD", "HF", "stroke", "death")

# -----------------------------------------------------------------------------
# 4.1 Cox比例风险回归模型（Table 2 - 核心结果表）
# -----------------------------------------------------------------------------

# 定义三个逐步调整的模型
# Model 1: 基础调整（年龄 + 性别 + 种族）
# Model 2: Model 1 + 生活方式 + 主要CVD风险因素
# Model 3: Model 2 + 生物标志物 + 药物使用 + LVEF

cox_results_all <- list()

for (outcome in outcomes) {
  
  cat("\n----- 分析结局:", outcome, "-----\n")
  
  # 构建变量名
  event_var <- paste0(outcome, "_event")
  time_var <- paste0(outcome, "_time")
  
  # Model 1: 基础模型
  model1_formula <- as.formula(paste0(
    "Surv(", time_var, ", ", event_var, ") ~ LV_size_category + age + sex + ethnicity"
  ))
  
  cox_model1 <- coxph(model1_formula, data = data)
  
  # Model 2: 扩展模型
  model2_formula <- as.formula(paste0(
    "Surv(", time_var, ", ", event_var, ") ~ LV_size_category + age + sex + ethnicity + ",
    "education + smoking_status + alcohol_status + BMI + SBP + DBP + ",
    "total_cholesterol + hypertension + diabetes"
  ))
  
  cox_model2 <- coxph(model2_formula, data = data)
  
  # Model 3: 完全调整模型
  model3_formula <- as.formula(paste0(
    "Surv(", time_var, ", ", event_var, ") ~ LV_size_category + age + sex + ethnicity + ",
    "education + smoking_status + alcohol_status + BMI + SBP + DBP + ",
    "total_cholesterol + hypertension + diabetes + HbA1c + CRP + ",
    "hemoglobin + albumin + atrial_fibrillation + statin_use + ",
    "antihypertensive_use + LVEF"
  ))
  
  cox_model3 <- coxph(model3_formula, data = data)
  
  # 提取结果
  extract_cox_results <- function(model, model_name) {
    coef_summary <- summary(model)$conf.int
    
    # 提取小LV和大LV的HR
    small_LV_row <- grep("Small", rownames(coef_summary))
    large_LV_row <- grep("Large", rownames(coef_summary))
    
    results <- data.frame(
      Outcome = outcome,
      Model = model_name,
      Group = c("Small LV", "Large LV"),
      HR = c(coef_summary[small_LV_row, "exp(coef)"], 
             coef_summary[large_LV_row, "exp(coef)"]),
      CI_lower = c(coef_summary[small_LV_row, "lower .95"], 
                   coef_summary[large_LV_row, "lower .95"]),
      CI_upper = c(coef_summary[small_LV_row, "upper .95"], 
                   coef_summary[large_LV_row, "upper .95"]),
      P_value = c(summary(model)$coefficients[small_LV_row, "Pr(>|z|)"],
                  summary(model)$coefficients[large_LV_row, "Pr(>|z|)"])
    )
    return(results)
  }
  
  # 合并三个模型的结果
  results_model1 <- extract_cox_results(cox_model1, "Model 1")
  results_model2 <- extract_cox_results(cox_model2, "Model 2")
  results_model3 <- extract_cox_results(cox_model3, "Model 3")
  
  outcome_results <- rbind(results_model1, results_model2, results_model3)
  cox_results_all[[outcome]] <- outcome_results
  
  print(outcome_results)
}

# 合并所有结局的结果
table2_cox_results <- do.call(rbind, cox_results_all)
write.csv(table2_cox_results, "Table2_Cox_Regression_Results.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# 4.2 Kaplan-Meier生存曲线（Figure 1A, Figure 2）
# -----------------------------------------------------------------------------

cat("\n----- 4.2 Kaplan-Meier生存曲线 -----\n")

# 创建保存图形的文件夹
if (!dir.exists("Figures")) dir.create("Figures")

# 绘制MACE的KM曲线（Figure 1A）
fit_MACE <- survfit(Surv(MACE_time, MACE_event) ~ LV_size_category, data = data)

km_plot_MACE <- ggsurvplot(
  fit_MACE,
  data = data,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  ggtheme = theme_bw(),
  palette = c("#E64B35", "#4DBBD5", "#00A087"),  # Small=红, Normal=蓝, Large=绿
  legend.title = "LV Size",
  legend.labs = c("Small", "Normal", "Large"),
  xlab = "Time (years)",
  ylab = "Cumulative Incidence of MACE",
  title = "Figure 1A. Kaplan-Meier Curves for MACE by LV Size"
)

# 保存图形
ggsave("Figures/Figure1A_KM_MACE.pdf", 
       print(km_plot_MACE), 
       width = 10, height = 8)

# 绘制其他结局的KM曲线（Figure 2: A-D）
other_outcomes <- c("CHD", "HF", "stroke", "death")
outcome_labels <- c("Coronary Heart Disease", "Heart Failure", 
                   "Ischemic Stroke", "All-Cause Mortality")

km_plots_list <- list()

for (i in 1:length(other_outcomes)) {
  outcome <- other_outcomes[i]
  label <- outcome_labels[i]
  
  event_var <- paste0(outcome, "_event")
  time_var <- paste0(outcome, "_time")
  
  fit <- survfit(as.formula(paste0("Surv(", time_var, ", ", event_var, 
                                   ") ~ LV_size_category")), 
                data = data)
  
  km_plot <- ggsurvplot(
    fit,
    data = data,
    pval = TRUE,
    pval.method = TRUE,
    conf.int = FALSE,
    risk.table = FALSE,
    ggtheme = theme_bw(),
    palette = c("#E64B35", "#4DBBD5", "#00A087"),
    legend.title = "LV Size",
    legend.labs = c("Small", "Normal", "Large"),
    xlab = "Time (years)",
    ylab = paste0("Cumulative Incidence of ", label),
    title = paste0("Panel ", LETTERS[i], ". ", label)
  )
  
  km_plots_list[[i]] <- km_plot$plot
  
  # 保存单独图形
  ggsave(paste0("Figures/Figure2", LETTERS[i], "_KM_", outcome, ".pdf"),
         print(km_plot),
         width = 8, height = 6)
}

# 合并四个图形为Figure 2
figure2_combined <- ggarrange(plotlist = km_plots_list, 
                              ncol = 2, nrow = 2,
                              common.legend = TRUE,
                              legend = "bottom")

ggsave("Figures/Figure2_KM_All_Outcomes.pdf", 
       figure2_combined, 
       width = 14, height = 12)


# =============================================================================
# 阶段三：深度模式挖掘 - "U型"关系验证
# =============================================================================

cat("\n========== 阶段三：U型关系验证 ==========\n")

# -----------------------------------------------------------------------------
# 5.1 五分位数分析（Table 3）
# -----------------------------------------------------------------------------

cat("\n----- 5.1 五分位数分析 -----\n")

quintile_results_all <- list()

for (outcome in outcomes) {
  
  event_var <- paste0(outcome, "_event")
  time_var <- paste0(outcome, "_time")
  
  # 以第3组（中间组）为参照
  data$LVEDVi_quintile <- relevel(data$LVEDVi_quintile, ref = "3")
  
  # Model 3（完全调整）
  quintile_formula <- as.formula(paste0(
    "Surv(", time_var, ", ", event_var, ") ~ LVEDVi_quintile + age + sex + ethnicity + ",
    "education + smoking_status + alcohol_status + BMI + SBP + DBP + ",
    "total_cholesterol + hypertension + diabetes + HbA1c + CRP + ",
    "hemoglobin + albumin + atrial_fibrillation + statin_use + ",
    "antihypertensive_use + LVEF"
  ))
  
  cox_quintile <- coxph(quintile_formula, data = data)
  
  # 提取结果
  coef_summary <- summary(cox_quintile)$conf.int
  quintile_rows <- grep("LVEDVi_quintile", rownames(coef_summary))
  
  results <- data.frame(
    Outcome = outcome,
    Quintile = c("Q1", "Q2", "Q3 (Ref)", "Q4", "Q5"),
    HR = c(coef_summary[quintile_rows[1:2], "exp(coef)"], 1.0,
           coef_summary[quintile_rows[3:4], "exp(coef)"]),
    CI_lower = c(coef_summary[quintile_rows[1:2], "lower .95"], NA,
                 coef_summary[quintile_rows[3:4], "lower .95"]),
    CI_upper = c(coef_summary[quintile_rows[1:2], "upper .95"], NA,
                 coef_summary[quintile_rows[3:4], "upper .95"]),
    P_value = c(summary(cox_quintile)$coefficients[quintile_rows[1:2], "Pr(>|z|)"], 
                NA,
                summary(cox_quintile)$coefficients[quintile_rows[3:4], "Pr(>|z|)"])
  )
  
  quintile_results_all[[outcome]] <- results
  print(results)
}

table3_quintile <- do.call(rbind, quintile_results_all)
write.csv(table3_quintile, "Table3_Quintile_Analysis.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# 5.2 限制性立方样条（RCS）分析（Figure 1B, Figure 3）
# -----------------------------------------------------------------------------

cat("\n----- 5.2 限制性立方样条（RCS）分析 -----\n")

# 使用rms包进行RCS分析
# 设置数据分布用于rms
dd <- datadist(data)
options(datadist = "dd")

rcs_plots_list <- list()

# 绘制MACE的RCS曲线（Figure 1B）
cat("\n绘制 MACE 的 RCS 曲线...\n")

# 构建完全调整的Cox模型，使用RCS
cox_rcs_MACE <- cph(
  Surv(MACE_time, MACE_event) ~ rcs(LVEDVi, 4) + age + sex + ethnicity +
    education + smoking_status + alcohol_status + BMI + SBP + DBP +
    total_cholesterol + hypertension + diabetes + HbA1c + CRP +
    hemoglobin + albumin + atrial_fibrillation + statin_use +
    antihypertensive_use + LVEF,
  data = data,
  x = TRUE,
  y = TRUE
)

# 预测HR和置信区间
LVEDVi_range <- seq(quantile(data$LVEDVi, 0.01), 
                    quantile(data$LVEDVi, 0.99), 
                    length.out = 100)

# 设置参考值（通常为中位数）
ref_value <- median(data$LVEDVi)

pred_MACE <- Predict(cox_rcs_MACE, LVEDVi = LVEDVi_range, ref.zero = TRUE)

# 转换为HR（exp(log HR)）
pred_MACE_df <- data.frame(
  LVEDVi = pred_MACE$LVEDVi,
  HR = exp(pred_MACE$yhat),
  CI_lower = exp(pred_MACE$lower),
  CI_upper = exp(pred_MACE$upper)
)

# 绘制图形
rcs_plot_MACE <- ggplot(pred_MACE_df, aes(x = LVEDVi, y = HR)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), 
              fill = "lightblue", alpha = 0.3) +
  geom_line(color = "blue", size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    x = "LVEDVi (mL/m²)",
    y = "Hazard Ratio for MACE",
    title = "Figure 1B. Adjusted Association Between LVEDVi and MACE (RCS)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

print(rcs_plot_MACE)
ggsave("Figures/Figure1B_RCS_MACE.pdf", rcs_plot_MACE, width = 8, height = 6)

# 绘制其他结局的RCS曲线（Figure 3: A-D）
for (i in 1:length(other_outcomes)) {
  outcome <- other_outcomes[i]
  label <- outcome_labels[i]
  
  cat("\n绘制", outcome, "的 RCS 曲线...\n")
  
  event_var <- paste0(outcome, "_event")
  time_var <- paste0(outcome, "_time")
  
  # 构建RCS模型
  formula_rcs <- as.formula(paste0(
    "Surv(", time_var, ", ", event_var, ") ~ rcs(LVEDVi, 4) + age + sex + ethnicity + ",
    "education + smoking_status + alcohol_status + BMI + SBP + DBP + ",
    "total_cholesterol + hypertension + diabetes + HbA1c + CRP + ",
    "hemoglobin + albumin + atrial_fibrillation + statin_use + ",
    "antihypertensive_use + LVEF"
  ))
  
  cox_rcs <- cph(formula_rcs, data = data, x = TRUE, y = TRUE)
  
  # 预测
  pred <- Predict(cox_rcs, LVEDVi = LVEDVi_range, ref.zero = TRUE)
  
  pred_df <- data.frame(
    LVEDVi = pred$LVEDVi,
    HR = exp(pred$yhat),
    CI_lower = exp(pred$lower),
    CI_upper = exp(pred$upper)
  )
  
  # 绘图
  rcs_plot <- ggplot(pred_df, aes(x = LVEDVi, y = HR)) +
    geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), 
                fill = "lightblue", alpha = 0.3) +
    geom_line(color = "blue", size = 1.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(
      x = "LVEDVi (mL/m²)",
      y = paste0("Hazard Ratio for ", label),
      title = paste0("Panel ", LETTERS[i], ". ", label)
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9)
    )
  
  rcs_plots_list[[i]] <- rcs_plot
  
  # 保存单独图形
  ggsave(paste0("Figures/Figure3", LETTERS[i], "_RCS_", outcome, ".pdf"),
         rcs_plot, width = 7, height = 5)
}

# 合并为Figure 3
figure3_combined <- ggarrange(plotlist = rcs_plots_list,
                              ncol = 2, nrow = 2)

ggsave("Figures/Figure3_RCS_All_Outcomes.pdf",
       figure3_combined,
       width = 14, height = 12)


# =============================================================================
# 阶段四：稳健性验证 - 敏感性分析
# =============================================================================

cat("\n========== 阶段四：敏感性分析 ==========\n")

# -----------------------------------------------------------------------------
# 6.1 敏感性分析1：使用UK Biobank特有切点值（Table S3）
# -----------------------------------------------------------------------------

cat("\n----- 6.1 使用UKB内部切点值 -----\n")

# 筛选健康参考队列（无CVD、高血压、糖尿病等）
healthy_cohort <- data %>%
  filter(
    prior_CVD == 0,
    hypertension == 0,
    diabetes == 0,
    BMI >= 18.5 & BMI < 30,
    smoking_status == "Never"
  )

# 计算性别特异性的2.5和97.5百分位数
cutoffs_male <- quantile(healthy_cohort$LVEDVi[healthy_cohort$sex == "Male"], 
                        c(0.025, 0.975), na.rm = TRUE)
cutoffs_female <- quantile(healthy_cohort$LVEDVi[healthy_cohort$sex == "Female"], 
                          c(0.025, 0.975), na.rm = TRUE)

cat("Male cutoffs:", cutoffs_male, "\n")
cat("Female cutoffs:", cutoffs_female, "\n")

# 使用新切点重新分组
data <- data %>%
  mutate(
    LV_size_UKB = case_when(
      sex == "Male" & LVEDVi < cutoffs_male[1] ~ "Small",
      sex == "Male" & LVEDVi >= cutoffs_male[1] & LVEDVi <= cutoffs_male[2] ~ "Normal",
      sex == "Male" & LVEDVi > cutoffs_male[2] ~ "Large",
      sex == "Female" & LVEDVi < cutoffs_female[1] ~ "Small",
      sex == "Female" & LVEDVi >= cutoffs_female[1] & LVEDVi <= cutoffs_female[2] ~ "Normal",
      sex == "Female" & LVEDVi > cutoffs_female[2] ~ "Large"
    ),
    LV_size_UKB = factor(LV_size_UKB, levels = c("Small", "Normal", "Large"))
  )

# 重复主分析
sensitivity1_results <- list()

for (outcome in outcomes) {
  event_var <- paste0(outcome, "_event")
  time_var <- paste0(outcome, "_time")
  
  model_formula <- as.formula(paste0(
    "Surv(", time_var, ", ", event_var, ") ~ LV_size_UKB + age + sex + ethnicity + ",
    "education + smoking_status + alcohol_status + BMI + SBP + DBP + ",
    "total_cholesterol + hypertension + diabetes + HbA1c + CRP + ",
    "hemoglobin + albumin + atrial_fibrillation + statin_use + ",
    "antihypertensive_use + LVEF"
  ))
  
  cox_model <- coxph(model_formula, data = data)
  
  coef_summary <- summary(cox_model)$conf.int
  small_row <- grep("Small", rownames(coef_summary))
  large_row <- grep("Large", rownames(coef_summary))
  
  results <- data.frame(
    Outcome = outcome,
    Group = c("Small LV", "Large LV"),
    HR = c(coef_summary[small_row, "exp(coef)"], 
           coef_summary[large_row, "exp(coef)"]),
    CI_lower = c(coef_summary[small_row, "lower .95"], 
                 coef_summary[large_row, "lower .95"]),
    CI_upper = c(coef_summary[small_row, "upper .95"], 
                 coef_summary[large_row, "upper .95"]),
    P_value = c(summary(cox_model)$coefficients[small_row, "Pr(>|z|)"],
                summary(cox_model)$coefficients[large_row, "Pr(>|z|)"])
  )
  
  sensitivity1_results[[outcome]] <- results
}

tableS3 <- do.call(rbind, sensitivity1_results)
write.csv(tableS3, "TableS3_Sensitivity_UKB_Cutoffs.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# 6.2 敏感性分析2：竞争风险模型（Table S4）
# -----------------------------------------------------------------------------

cat("\n----- 6.2 竞争风险分析（Fine-Gray模型）-----\n")

# 使用Fine-Gray模型分析心血管结局，将非CVD死亡作为竞争风险

# 创建竞争风险变量（0=未发生, 1=目标事件, 2=竞争风险事件）
# 假设有non_CVD_death变量

sensitivity2_results <- list()

for (outcome in c("MACE", "CHD", "HF", "stroke")) {
  
  cat("\n分析结局:", outcome, "（竞争风险模型）\n")
  
  event_var <- paste0(outcome, "_event")
  time_var <- paste0(outcome, "_time")
  
  # 创建竞争风险状态变量
  data$competing_status <- with(data, {
    status <- rep(0, nrow(data))
    status[get(event_var) == 1] <- 1  # 目标事件
    status[non_CVD_death == 1 & get(event_var) == 0] <- 2  # 竞争风险
    status
  })
  
  # 构建协变量矩阵
  covariates <- c("age", "sex", "ethnicity", "education", "smoking_status", 
                 "alcohol_status", "BMI", "SBP", "DBP", "total_cholesterol",
                 "hypertension", "diabetes", "HbA1c", "CRP", "hemoglobin",
                 "albumin", "atrial_fibrillation", "statin_use", 
                 "antihypertensive_use", "LVEF")
  
  # 创建哑变量
  data$LV_small <- ifelse(data$LV_size_category == "Small", 1, 0)
  data$LV_large <- ifelse(data$LV_size_category == "Large", 1, 0)
  
  # Fine-Gray模型（小LV vs 正常LV）
  fg_model_small <- crr(
    ftime = data[[time_var]],
    fstatus = data$competing_status,
    cov1 = data[, c("LV_small", covariates)],
    failcode = 1,
    cencode = 0
  )
  
  # Fine-Gray模型（大LV vs 正常LV）
  fg_model_large <- crr(
    ftime = data[[time_var]],
    fstatus = data$competing_status,
    cov1 = data[, c("LV_large", covariates)],
    failcode = 1,
    cencode = 0
  )
  
  # 提取结果
  results <- data.frame(
    Outcome = outcome,
    Group = c("Small LV", "Large LV"),
    HR = c(exp(fg_model_small$coef[1]), exp(fg_model_large$coef[1])),
    CI_lower = c(exp(fg_model_small$coef[1] - 1.96 * sqrt(fg_model_small$var[1,1])),
                 exp(fg_model_large$coef[1] - 1.96 * sqrt(fg_model_large$var[1,1]))),
    CI_upper = c(exp(fg_model_small$coef[1] + 1.96 * sqrt(fg_model_small$var[1,1])),
                 exp(fg_model_large$coef[1] + 1.96 * sqrt(fg_model_large$var[1,1])))
  )
  
  sensitivity2_results[[outcome]] <- results
}

tableS4 <- do.call(rbind, sensitivity2_results)
write.csv(tableS4, "TableS4_Sensitivity_Competing_Risk.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# 6.3 敏感性分析3：排除基线已患CVD者（Table S5）
# -----------------------------------------------------------------------------

cat("\n----- 6.3 排除基线CVD患者 -----\n")

# 筛选无基线CVD的人群
data_no_CVD <- data %>% filter(prior_CVD == 0)

sensitivity3_results <- list()

for (outcome in outcomes) {
  event_var <- paste0(outcome, "_event")
  time_var <- paste0(outcome, "_time")
  
  model_formula <- as.formula(paste0(
    "Surv(", time_var, ", ", event_var, ") ~ LV_size_category + age + sex + ethnicity + ",
    "education + smoking_status + alcohol_status + BMI + SBP + DBP + ",
    "total_cholesterol + hypertension + diabetes + HbA1c + CRP + ",
    "hemoglobin + albumin + atrial_fibrillation + statin_use + ",
    "antihypertensive_use + LVEF"
  ))
  
  cox_model <- coxph(model_formula, data = data_no_CVD)
  
  coef_summary <- summary(cox_model)$conf.int
  small_row <- grep("Small", rownames(coef_summary))
  large_row <- grep("Large", rownames(coef_summary))
  
  results <- data.frame(
    Outcome = outcome,
    Group = c("Small LV", "Large LV"),
    HR = c(coef_summary[small_row, "exp(coef)"], 
           coef_summary[large_row, "exp(coef)"]),
    CI_lower = c(coef_summary[small_row, "lower .95"], 
                 coef_summary[large_row, "lower .95"]),
    CI_upper = c(coef_summary[small_row, "upper .95"], 
                 coef_summary[large_row, "upper .95"]),
    P_value = c(summary(cox_model)$coefficients[small_row, "Pr(>|z|)"],
                summary(cox_model)$coefficients[large_row, "Pr(>|z|)"])
  )
  
  sensitivity3_results[[outcome]] <- results
}

tableS5 <- do.call(rbind, sensitivity3_results)
write.csv(tableS5, "TableS5_Sensitivity_Exclude_Baseline_CVD.csv", row.names = FALSE)


# =============================================================================
# 阶段五：结论普适性探讨 - 分层亚组分析
# =============================================================================

cat("\n========== 阶段五：分层亚组分析 ==========\n")

# -----------------------------------------------------------------------------
# 7.1 亚组分析与交互作用检验（Figure 4 - 森林图）
# -----------------------------------------------------------------------------

cat("\n----- 7.1 亚组分析（MACE结局）-----\n")

# 定义亚组变量
subgroup_vars <- c("sex", "age_group", "hypertension", "diabetes", 
                  "BMI_category", "smoking_status")

# 创建亚组分类变量
data <- data %>%
  mutate(
    age_group = ifelse(age < 60, "<60 years", "≥60 years"),
    BMI_category = case_when(
      BMI < 25 ~ "Normal weight",
      BMI >= 25 & BMI < 30 ~ "Overweight",
      BMI >= 30 ~ "Obese"
    )
  )

subgroup_results <- data.frame()

for (subvar in subgroup_vars) {
  
  cat("\n亚组变量:", subvar, "\n")
  
  # 获取亚组水平
  sublevels <- unique(data[[subvar]])
  
  for (sublevel in sublevels) {
    
    # 筛选亚组数据
    data_sub <- data %>% filter(get(subvar) == sublevel)
    
    # Model 3
    cox_sub <- coxph(
      Surv(MACE_time, MACE_event) ~ LV_size_category + age + sex + ethnicity +
        education + smoking_status + alcohol_status + BMI + SBP + DBP +
        total_cholesterol + hypertension + diabetes + HbA1c + CRP +
        hemoglobin + albumin + atrial_fibrillation + statin_use +
        antihypertensive_use + LVEF,
      data = data_sub
    )
    
    coef_summary <- summary(cox_sub)$conf.int
    small_row <- grep("Small", rownames(coef_summary))
    large_row <- grep("Large", rownames(coef_summary))
    
    # 提取小LV的结果
    result_small <- data.frame(
      Subgroup_var = subvar,
      Subgroup_level = sublevel,
      LV_group = "Small LV",
      HR = coef_summary[small_row, "exp(coef)"],
      CI_lower = coef_summary[small_row, "lower .95"],
      CI_upper = coef_summary[small_row, "upper .95"]
    )
    
    # 提取大LV的结果
    result_large <- data.frame(
      Subgroup_var = subvar,
      Subgroup_level = sublevel,
      LV_group = "Large LV",
      HR = coef_summary[large_row, "exp(coef)"],
      CI_lower = coef_summary[large_row, "lower .95"],
      CI_upper = coef_summary[large_row, "upper .95"]
    )
    
    subgroup_results <- rbind(subgroup_results, result_small, result_large)
  }
  
  # 计算交互作用P值
  interaction_formula <- as.formula(paste0(
    "Surv(MACE_time, MACE_event) ~ LV_size_category * ", subvar, 
    " + age + sex + ethnicity + education + smoking_status + alcohol_status + ",
    "BMI + SBP + DBP + total_cholesterol + hypertension + diabetes + ",
    "HbA1c + CRP + hemoglobin + albumin + atrial_fibrillation + ",
    "statin_use + antihypertensive_use + LVEF"
  ))
  
  cox_interaction <- coxph(interaction_formula, data = data)
  
  # 提取交互项P值
  interaction_terms <- grep(":", names(coef(cox_interaction)), value = TRUE)
  p_interaction <- summary(cox_interaction)$coefficients[interaction_terms[1], "Pr(>|z|)"]
  
  cat("交互作用 P 值:", p_interaction, "\n")
}

write.csv(subgroup_results, "Subgroup_Analysis_Results.csv", row.names = FALSE)

# 绘制森林图（Figure 4）
# 这里使用forestplot包绘制
# 需要根据实际数据格式调整

library(forestplot)

# 准备森林图数据（示例）
forest_data <- subgroup_results %>%
  filter(LV_group == "Small LV") %>%
  mutate(
    HR_text = sprintf("%.2f (%.2f-%.2f)", HR, CI_lower, CI_upper),
    Subgroup = paste(Subgroup_var, Subgroup_level, sep = ": ")
  )

# 绘制森林图（简化版）
pdf("Figures/Figure4_Forest_Plot_Subgroups.pdf", width = 12, height = 10)

forestplot(
  labeltext = as.matrix(forest_data[, c("Subgroup", "HR_text")]),
  mean = forest_data$HR,
  lower = forest_data$CI_lower,
  upper = forest_data$CI_upper,
  xlog = TRUE,
  xlab = "Hazard Ratio",
  title = "Figure 4. Subgroup Analysis for MACE (Small LV vs Normal LV)"
)

dev.off()


# -----------------------------------------------------------------------------
# 7.2 按性别分层的完整分析（Table S6, Table S7, Figure 5）
# -----------------------------------------------------------------------------

cat("\n----- 7.2 按性别分层的完整分析 -----\n")

# 7.2.1 按性别分层的基线特征（Table S6）
tableS6_male <- CreateTableOne(
  vars = baseline_vars,
  strata = "LV_size_category",
  data = data %>% filter(sex == "Male"),
  test = TRUE
)

tableS6_female <- CreateTableOne(
  vars = baseline_vars,
  strata = "LV_size_category",
  data = data %>% filter(sex == "Female"),
  test = TRUE
)

write.csv(print(tableS6_male, printToggle = FALSE), 
         "TableS6_Baseline_Male.csv")
write.csv(print(tableS6_female, printToggle = FALSE), 
         "TableS6_Baseline_Female.csv")

# 7.2.2 按性别分层的Cox回归（Table S7）
sex_stratified_results <- list()

for (sex_group in c("Male", "Female")) {
  
  data_sex <- data %>% filter(sex == sex_group)
  
  for (outcome in outcomes) {
    
    event_var <- paste0(outcome, "_event")
    time_var <- paste0(outcome, "_time")
    
    # Model 3
    model_formula <- as.formula(paste0(
      "Surv(", time_var, ", ", event_var, ") ~ LV_size_category + age + ethnicity + ",
      "education + smoking_status + alcohol_status + BMI + SBP + DBP + ",
      "total_cholesterol + hypertension + diabetes + HbA1c + CRP + ",
      "hemoglobin + albumin + atrial_fibrillation + statin_use + ",
      "antihypertensive_use + LVEF"
    ))
    
    cox_model <- coxph(model_formula, data = data_sex)
    
    coef_summary <- summary(cox_model)$conf.int
    small_row <- grep("Small", rownames(coef_summary))
    large_row <- grep("Large", rownames(coef_summary))
    
    results <- data.frame(
      Sex = sex_group,
      Outcome = outcome,
      Group = c("Small LV", "Large LV"),
      HR = c(coef_summary[small_row, "exp(coef)"], 
             coef_summary[large_row, "exp(coef)"]),
      CI_lower = c(coef_summary[small_row, "lower .95"], 
                   coef_summary[large_row, "lower .95"]),
      CI_upper = c(coef_summary[small_row, "upper .95"], 
                   coef_summary[large_row, "upper .95"]),
      P_value = c(summary(cox_model)$coefficients[small_row, "Pr(>|z|)"],
                  summary(cox_model)$coefficients[large_row, "Pr(>|z|)"])
    )
    
    sex_stratified_results[[paste(sex_group, outcome, sep = "_")]] <- results
  }
}

tableS7 <- do.call(rbind, sex_stratified_results)
write.csv(tableS7, "TableS7_Sex_Stratified_Cox_Results.csv", row.names = FALSE)

# 7.2.3 按性别分层的RCS曲线（Figure 5B - 最重要的图）
cat("\n绘制按性别分层的RCS曲线...\n")

# 男性
data_male <- data %>% filter(sex == "Male")
dd_male <- datadist(data_male)
options(datadist = "dd_male")

cox_rcs_male <- cph(
  Surv(MACE_time, MACE_event) ~ rcs(LVEDVi, 4) + age + ethnicity +
    education + smoking_status + alcohol_status + BMI + SBP + DBP +
    total_cholesterol + hypertension + diabetes + HbA1c + CRP +
    hemoglobin + albumin + atrial_fibrillation + statin_use +
    antihypertensive_use + LVEF,
  data = data_male,
  x = TRUE,
  y = TRUE
)

pred_male <- Predict(cox_rcs_male, LVEDVi = LVEDVi_range, ref.zero = TRUE)
pred_male_df <- data.frame(
  LVEDVi = pred_male$LVEDVi,
  HR = exp(pred_male$yhat),
  CI_lower = exp(pred_male$lower),
  CI_upper = exp(pred_male$upper),
  Sex = "Male"
)

# 女性
data_female <- data %>% filter(sex == "Female")
dd_female <- datadist(data_female)
options(datadist = "dd_female")

cox_rcs_female <- cph(
  Surv(MACE_time, MACE_event) ~ rcs(LVEDVi, 4) + age + ethnicity +
    education + smoking_status + alcohol_status + BMI + SBP + DBP +
    total_cholesterol + hypertension + diabetes + HbA1c + CRP +
    hemoglobin + albumin + atrial_fibrillation + statin_use +
    antihypertensive_use + LVEF,
  data = data_female,
  x = TRUE,
  y = TRUE
)

pred_female <- Predict(cox_rcs_female, LVEDVi = LVEDVi_range, ref.zero = TRUE)
pred_female_df <- data.frame(
  LVEDVi = pred_female$LVEDVi,
  HR = exp(pred_female$yhat),
  CI_lower = exp(pred_female$lower),
  CI_upper = exp(pred_female$upper),
  Sex = "Female"
)

# 合并数据
pred_sex_combined <- rbind(pred_male_df, pred_female_df)

# 绘制Figure 5B
figure5B <- ggplot(pred_sex_combined, aes(x = LVEDVi, y = HR, color = Sex, fill = Sex)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2, linetype = 0) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Male" = "#0073C2", "Female" = "#EFC000")) +
  scale_fill_manual(values = c("Male" = "#0073C2", "Female" = "#EFC000")) +
  labs(
    x = "LVEDVi (mL/m²)",
    y = "Hazard Ratio for MACE",
    title = "Figure 5B. Sex-Specific Association Between LVEDVi and MACE"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  )

print(figure5B)
ggsave("Figures/Figure5B_RCS_Sex_Stratified.pdf", figure5B, width = 10, height = 7)


# =============================================================================
# 第六部分：补充分析
# =============================================================================

cat("\n========== 补充分析 ==========\n")

# -----------------------------------------------------------------------------
# 8.1 小LV和大LV的患病率与年龄的关系（Figure S2）
# -----------------------------------------------------------------------------

prevalence_by_age <- data %>%
  mutate(age_group_5yr = cut(age, breaks = seq(40, 80, by = 5))) %>%
  group_by(age_group_5yr, LV_size_category) %>%
  summarise(n = n()) %>%
  group_by(age_group_5yr) %>%
  mutate(prevalence = n / sum(n) * 100)

figureS2 <- ggplot(prevalence_by_age, 
                   aes(x = age_group_5yr, y = prevalence, fill = LV_size_category)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Small" = "#E64B35", 
                               "Normal" = "#4DBBD5", 
                               "Large" = "#00A087")) +
  labs(
    x = "Age Group (years)",
    y = "Prevalence (%)",
    fill = "LV Size",
    title = "Figure S2. Prevalence of LV Size Categories by Age Group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("Figures/FigureS2_Prevalence_by_Age.pdf", figureS2, width = 10, height = 6)


# -----------------------------------------------------------------------------
# 8.2 LV大小与其他结局的关联（Table S2）
# -----------------------------------------------------------------------------

# 分析次要结局（如房颤、非CVD死亡等）
secondary_outcomes <- c("atrial_fib_incident", "non_CVD_death")

tableS2_results <- list()

for (outcome in secondary_outcomes) {
  
  event_var <- paste0(outcome, "_event")
  time_var <- paste0(outcome, "_time")
  
  model_formula <- as.formula(paste0(
    "Surv(", time_var, ", ", event_var, ") ~ LV_size_category + age + sex + ethnicity + ",
    "education + smoking_status + alcohol_status + BMI + SBP + DBP + ",
    "total_cholesterol + hypertension + diabetes + HbA1c + CRP + ",
    "hemoglobin + albumin + atrial_fibrillation + statin_use + ",
    "antihypertensive_use + LVEF"
  ))
  
  cox_model <- coxph(model_formula, data = data)
  
  coef_summary <- summary(cox_model)$conf.int
  small_row <- grep("Small", rownames(coef_summary))
  large_row <- grep("Large", rownames(coef_summary))
  
  results <- data.frame(
    Outcome = outcome,
    Group = c("Small LV", "Large LV"),
    HR = c(coef_summary[small_row, "exp(coef)"], 
           coef_summary[large_row, "exp(coef)"]),
    CI_lower = c(coef_summary[small_row, "lower .95"], 
                 coef_summary[large_row, "lower .95"]),
    CI_upper = c(coef_summary[small_row, "upper .95"], 
                 coef_summary[large_row, "upper .95"]),
    P_value = c(summary(cox_model)$coefficients[small_row, "Pr(>|z|)"],
                summary(cox_model)$coefficients[large_row, "Pr(>|z|)"])
  )
  
  tableS2_results[[outcome]] <- results
}

tableS2 <- do.call(rbind, tableS2_results)
write.csv(tableS2, "TableS2_Secondary_Outcomes.csv", row.names = FALSE)


# =============================================================================
# 第七部分：生成综合报告
# =============================================================================

cat("\n========== 生成分析报告 ==========\n")

# 创建报告文件夹
if (!dir.exists("Reports")) dir.create("Reports")

# 保存会话信息
sink("Reports/Session_Info.txt")
cat("R Session Information\n")
cat("======================\n\n")
print(sessionInfo())
sink()

# 保存所有结果对象
save.image("Results_Complete.RData")

cat("\n\n")
cat("=======================================================\n")
cat("         所有分析已完成！                              \n")
cat("=======================================================\n")
cat("\n生成的文件包括：\n")
cat("  - 表格文件：Table1-3, TableS1-S7 (CSV格式)\n")
cat("  - 图形文件：Figures文件夹中的所有PDF图形\n")
cat("  - 完整数据：Results_Complete.RData\n")
cat("  - 会话信息：Reports/Session_Info.txt\n")
cat("\n感谢使用本复现代码！\n")
cat("=======================================================\n")

################################################################################
# 代码结束
################################################################################

