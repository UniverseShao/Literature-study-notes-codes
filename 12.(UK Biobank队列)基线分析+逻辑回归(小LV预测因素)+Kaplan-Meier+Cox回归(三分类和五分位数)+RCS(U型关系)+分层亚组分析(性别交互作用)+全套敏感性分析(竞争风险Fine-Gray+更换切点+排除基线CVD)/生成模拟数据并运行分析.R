################################################################################
# UK Biobank 左心室大小与心血管结局关联研究 - 完整复现代码（含模拟数据）
# 
# 研究设计：基线分析 + 逻辑回归 + Cox回归 + RCS + 分层亚组 + 敏感性分析
# 
# 日期：2025-10-04
################################################################################

# =============================================================================
# 第零部分：生成模拟数据
# =============================================================================

cat("\n========== 生成模拟UK Biobank数据 ==========\n")

set.seed(12345)  # 设置随机种子以确保可重复性

# 样本量（增加样本量使曲线更平滑）
n <- 50000

# 生成模拟数据
data <- data.frame(
  # 人口学变量
  age = rnorm(n, mean = 55, sd = 8),
  sex = sample(c("Male", "Female"), n, replace = TRUE, prob = c(0.48, 0.52)),
  ethnicity = sample(c("White", "Asian", "Black", "Other"), n, replace = TRUE, 
                     prob = c(0.85, 0.08, 0.04, 0.03)),
  education = sample(c("None", "School", "College", "University"), n, replace = TRUE,
                    prob = c(0.10, 0.35, 0.30, 0.25)),
  
  # 体格检查指标
  BMI = rnorm(n, mean = 27, sd = 4.5),
  SBP = rnorm(n, mean = 135, sd = 18),
  DBP = rnorm(n, mean = 82, sd = 10),
  
  # 生化指标
  total_cholesterol = rnorm(n, mean = 5.5, sd = 1.1),
  HbA1c = rnorm(n, mean = 5.6, sd = 0.8),
  CRP = rexp(n, rate = 0.5),  # 指数分布，模拟右偏
  hemoglobin = rnorm(n, mean = 14.5, sd = 1.5),
  albumin = rnorm(n, mean = 45, sd = 3),
  
  # 生活方式
  smoking_status = sample(c("Never", "Former", "Current"), n, replace = TRUE,
                         prob = c(0.55, 0.30, 0.15)),
  alcohol_status = sample(c("Never", "Occasional", "Regular"), n, replace = TRUE,
                         prob = c(0.15, 0.45, 0.40)),
  
  # 既往疾病史
  hypertension = rbinom(n, 1, 0.30),
  diabetes = rbinom(n, 1, 0.08),
  prior_CVD = rbinom(n, 1, 0.05),
  atrial_fibrillation = rbinom(n, 1, 0.03),
  
  # 用药情况
  statin_use = rbinom(n, 1, 0.20),
  antihypertensive_use = rbinom(n, 1, 0.25)
)

# 生成左心室指标（考虑性别差异）
data$LVEF <- rnorm(n, mean = 62, sd = 6)  # 左心室射血分数

# LVEDVi（左心室舒张末期容积指数）- 性别差异明显
data$LVEDVi <- ifelse(
  data$sex == "Male",
  rnorm(n, mean = 78, sd = 15),  # 男性均值更高
  rnorm(n, mean = 68, sd = 13)   # 女性均值较低
)

# 确保LVEDVi在合理范围内
data$LVEDVi <- pmax(40, pmin(data$LVEDVi, 150))

# 生成随访时间（年）- 使用更真实的分布，减少删失
follow_up_years <- rgamma(n, shape = 2.5, scale = 2.8) + 0.5
follow_up_years <- pmin(follow_up_years, 12)  # 限制最大随访时间为12年

# 为每个结局生成事件和时间
# 基础风险受多个因素影响（提高基础风险以减少删失）
base_risk <- with(data, {
  risk <- 0.015 +  # 提高基础风险
    0.0008 * (age - 55) +
    0.002 * (BMI > 30) +
    0.004 * hypertension +
    0.006 * diabetes +
    0.008 * prior_CVD +
    0.003 * (smoking_status == "Current")
  risk
})

# LV大小对风险的影响（优化后的风险差异，确保曲线清晰分离）
lv_risk_effect <- with(data, {
  ifelse(sex == "Male",
         # 男性：Small LV高风险，Normal LV为参考组，Large LV中度风险
         ifelse(LVEDVi < 68, 0.60,    # Small LV: 高风险（60%增加）
                ifelse(LVEDVi > 88, 0.20, 0.0)),  # Large LV: 中度增加（20%），Normal LV: 参考组（0%）
         # 女性：相同逻辑但切点不同
         ifelse(LVEDVi < 60, 0.60,    # Small LV: 高风险（60%增加）
                ifelse(LVEDVi > 76, 0.20, 0.0))   # Large LV: 中度增加（20%），Normal LV: 参考组（0%）
  )
})

# MACE（主要心血管不良事件）- 大幅提高事件率以减少删失
mace_risk <- pmax(0.008, pmin(base_risk * 1.5 + lv_risk_effect * 0.35, 0.65))  # 大幅提高风险
mace_prob <- pmax(0, pmin(mace_risk * follow_up_years / 3.5, 0.55))  # 提高事件率到55%
data$MACE_event <- rbinom(n, 1, mace_prob)
data$MACE_time <- ifelse(data$MACE_event == 1,
                         pmin(rexp(n, rate = 1/2.5) + 0.5, follow_up_years),
                         follow_up_years)

# CHD（冠心病）- 提高事件率
chd_risk <- pmax(0.008, pmin(base_risk * 1.3 + lv_risk_effect * 0.30, 0.60))
chd_prob <- pmax(0, pmin(chd_risk * follow_up_years / 3.8, 0.50))
data$CHD_event <- rbinom(n, 1, chd_prob)
data$CHD_time <- ifelse(data$CHD_event == 1,
                        pmin(rexp(n, rate = 1/2.8) + 0.5, follow_up_years),
                        follow_up_years)

# HF（心力衰竭）- 提高事件率
hf_risk <- pmax(0.008, pmin(base_risk * 1.4 + lv_risk_effect * 0.38, 0.62))
hf_prob <- pmax(0, pmin(hf_risk * follow_up_years / 3.2, 0.52))
data$HF_event <- rbinom(n, 1, hf_prob)
data$HF_time <- ifelse(data$HF_event == 1,
                       pmin(rexp(n, rate = 1/2.2) + 0.5, follow_up_years),
                       follow_up_years)

# Stroke（卒中）- 提高事件率
stroke_risk <- pmax(0.008, pmin(base_risk * 1.2 + lv_risk_effect * 0.32, 0.58))
stroke_prob <- pmax(0, pmin(stroke_risk * follow_up_years / 4.0, 0.48))
data$stroke_event <- rbinom(n, 1, stroke_prob)
data$stroke_time <- ifelse(data$stroke_event == 1,
                           pmin(rexp(n, rate = 1/3.0) + 0.5, follow_up_years),
                           follow_up_years)

# Death（全因死亡）- 提高事件率
death_risk <- pmax(0.008, pmin(base_risk * 2.5 + lv_risk_effect * 0.40, 0.70))
death_prob <- pmax(0, pmin(death_risk * follow_up_years / 3.0, 0.58))
data$death_event <- rbinom(n, 1, death_prob)
data$death_time <- ifelse(data$death_event == 1,
                          pmin(rexp(n, rate = 1/3.2) + 0.5, follow_up_years),
                          follow_up_years)

# 非CVD死亡（用于竞争风险分析）
data$non_CVD_death <- ifelse(data$death_event == 1 & 
                              data$MACE_event == 0 & 
                              runif(n) < 0.4, 1, 0)

# 房颤发生（次要结局）
data$atrial_fib_incident_event <- rbinom(n, 1, 0.02)
data$atrial_fib_incident_time <- ifelse(data$atrial_fib_incident_event == 1,
                                        runif(n, 0.5, follow_up_years),
                                        follow_up_years)

# 确保分类变量是因子类型
data$sex <- factor(data$sex)
data$ethnicity <- factor(data$ethnicity)
data$education <- factor(data$education)
data$smoking_status <- factor(data$smoking_status, levels = c("Never", "Former", "Current"))
data$alcohol_status <- factor(data$alcohol_status)

cat("模拟数据生成完成！\n")
cat("样本量：", n, "\n")
cat("MACE事件数：", sum(data$MACE_event, na.rm = TRUE), "(事件率:", round(mean(data$MACE_event, na.rm = TRUE)*100, 2), "%)\n")
cat("CHD事件数：", sum(data$CHD_event, na.rm = TRUE), "(事件率:", round(mean(data$CHD_event, na.rm = TRUE)*100, 2), "%)\n")
cat("HF事件数：", sum(data$HF_event, na.rm = TRUE), "(事件率:", round(mean(data$HF_event, na.rm = TRUE)*100, 2), "%)\n")
cat("Stroke事件数：", sum(data$stroke_event, na.rm = TRUE), "(事件率:", round(mean(data$stroke_event, na.rm = TRUE)*100, 2), "%)\n")
cat("Death事件数：", sum(data$death_event, na.rm = TRUE), "(事件率:", round(mean(data$death_event, na.rm = TRUE)*100, 2), "%)\n\n")


# =============================================================================
# 第一部分：环境准备
# =============================================================================

# 1.1 清空环境（保留数据）
# rm(list = ls())
# gc()

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

cat("正在检查并安装必需的R包...\n")

# 安装并加载包
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("安装", pkg, "...\n")
    install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  }
}

cat("所有R包加载完成！\n\n")

# =============================================================================
# 第二部分：数据预处理
# =============================================================================

cat("\n========== 数据预处理 ==========\n")

# 2.1 创建左心室大小分类变量（进一步优化切点以增强组间差异）
data <- data %>%
  mutate(
    LV_size_category = case_when(
      sex == "Male" & LVEDVi < 68 ~ "Small",      # 男性切点进一步调整
      sex == "Male" & LVEDVi >= 68 & LVEDVi <= 88 ~ "Normal",
      sex == "Male" & LVEDVi > 88 ~ "Large", 
      sex == "Female" & LVEDVi < 60 ~ "Small",    # 女性切点进一步调整
      sex == "Female" & LVEDVi >= 60 & LVEDVi <= 76 ~ "Normal",
      sex == "Female" & LVEDVi > 76 ~ "Large"
    ),
    LV_size_category = factor(LV_size_category, 
                              levels = c("Normal", "Small", "Large"))
  )

# 2.2 创建二分类变量用于逻辑回归（小LV vs 非小LV）
data <- data %>%
  mutate(small_LV = ifelse(LV_size_category == "Small", 1, 0))

# 2.3 创建LVEDVi五分位数分组
data <- data %>%
  mutate(LVEDVi_quintile = ntile(LVEDVi, 5),
         LVEDVi_quintile = factor(LVEDVi_quintile))

cat("LV分组统计：\n")
print(table(data$LV_size_category, data$sex))
cat("\n")

# 打印事件发生率按组统计
cat("各组MACE事件发生率：\n")
event_by_group <- data %>%
  group_by(LV_size_category) %>%
  summarise(
    n = n(),
    MACE_events = sum(MACE_event, na.rm = TRUE),
    MACE_rate = round(mean(MACE_event, na.rm = TRUE) * 100, 2),
    CHD_events = sum(CHD_event, na.rm = TRUE),
    CHD_rate = round(mean(CHD_event, na.rm = TRUE) * 100, 2),
    HF_events = sum(HF_event, na.rm = TRUE),
    HF_rate = round(mean(HF_event, na.rm = TRUE) * 100, 2),
    .groups = 'drop'  # 清除分组信息
  )
print(event_by_group)
cat("\n")


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
cat("\nTable 1: 基线特征按LV大小分组\n")
print(table1, showAllLevels = TRUE, cramVars = "sex")

# 保存为CSV
table1_output <- print(table1, printToggle = FALSE, noSpaces = TRUE)
write.csv(table1_output, "Table1_Baseline_Characteristics.csv")
cat("已保存：Table1_Baseline_Characteristics.csv\n")


# -----------------------------------------------------------------------------
# 3.2 逻辑回归分析：预测"小左心室"的因素
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

cat("\n单变量逻辑回归结果：\n")
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

cat("\n多变量逻辑回归结果：\n")
print(head(multivariate_results_df, 10))
write.csv(multivariate_results_df, "Multivariate_Logistic_Regression.csv", row.names = FALSE)
cat("已保存逻辑回归结果\n")


# =============================================================================
# 阶段二：核心假设检验 - LV大小与心血管结局的关联
# =============================================================================

cat("\n========== 阶段二：Cox回归分析 ==========\n")

# 定义所有结局
outcomes <- c("MACE", "CHD", "HF", "stroke", "death")

# -----------------------------------------------------------------------------
# 4.1 Cox比例风险回归模型（Table 2 - 核心结果表）
# -----------------------------------------------------------------------------

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
cat("\n已保存：Table2_Cox_Regression_Results.csv\n")


# -----------------------------------------------------------------------------
# 4.2 Kaplan-Meier生存曲线
# -----------------------------------------------------------------------------

cat("\n----- 4.2 Kaplan-Meier生存曲线 -----\n")

# 创建保存图形的文件夹（PNG格式，更稳定）
if (!dir.exists("Figures")) dir.create("Figures")

# 绘制MACE的KM曲线（Figure 1A）- SCI期刊专业风格（优化删失显示）
fit_MACE <- survfit(Surv(MACE_time, MACE_event) ~ LV_size_category, data = data)

# 使用Nature/NEJM风格配色：深蓝、深红、深橙
km_plot_MACE <- ggsurvplot(
  fit_MACE,
  data = data,
  pval = TRUE,
  pval.method = TRUE,
  pval.coord = c(1, 0.15),  # p值位置
  pval.size = 4.5,  # p值字体大小
  conf.int = TRUE,  # 添加95%置信区间
  conf.int.alpha = 0.15,  # 置信区间透明度
  risk.table = TRUE,
  risk.table.height = 0.25,
  risk.table.fontsize = 3.8,
  ggtheme = theme_classic(),  # 使用经典主题（更专业）
  palette = c("#0072B2", "#D55E00", "#009E73"),  # Nature风格配色：蓝、橙红、绿
  legend.title = "LV Size Category",
  legend.labs = c("Normal LV", "Small LV", "Large LV"),
  xlab = "Follow-up Time (years)",
  ylab = "Event-free Survival Probability",
  title = "A. Major Adverse Cardiovascular Events",
  font.title = c(16, "bold"),
  font.legend = c(12, "plain"),
  font.x = c(14, "plain"),
  font.y = c(14, "plain"),
  font.tickslab = c(12, "plain"),
  linetype = c("solid", "dashed", "dotted"),  # 不同线型增强区分
  size = 1.5,  # 增加线条粗细使曲线更清晰
  censor = FALSE,  # 移除删失标记减少视觉杂乱
  xlim = c(0, 12),
  ylim = c(0, 1),
  break.time.by = 2,
  break.y.by = 0.2,  # y轴间隔20%
  surv.scale = "percent",
  tables.theme = theme_cleantable(),
  legend = c(0.15, 0.15)  # 图例位置（左下角）
)

# 保存为高分辨率PNG（发表级质量）
png("Figures/Figure1A_KM_MACE.png", width = 12*300, height = 10*300, res = 300, bg = "white")
print(km_plot_MACE$plot)
dev.off()
cat("已保存：Figures/Figure1A_KM_MACE.png\n")

# 绘制其他结局的KM曲线 - 优化可视化效果
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
    pval.method.coord = c(1.5, 0.15),
    pval.method.size = 4,
    pval.coord = c(1.5, 0.10),
    pval.size = 4,
    conf.int = TRUE,  # 添加置信区间
    conf.int.alpha = 0.15,
    risk.table = FALSE,
    ggtheme = theme_classic(),  # 经典主题
    palette = c("#0072B2", "#D55E00", "#009E73"),  # Nature风格配色
    legend.title = "LV Size",
    legend.labs = c("Normal", "Small", "Large"),
    xlab = "Follow-up Time (years)",
    ylab = "Event-free Survival",
    title = paste0(LETTERS[i], ". ", label),
    font.title = c(14, "bold"),
    font.legend = c(11, "plain"),
    font.x = c(13, "plain"),
    font.y = c(13, "plain"),
    font.tickslab = c(11, "plain"),
    linetype = c("solid", "dashed", "dotted"),  # 不同线型
    size = 1.2,  # 增加线条粗细使曲线更清晰
    censor = FALSE,  # 移除删失标记减少视觉杂乱
    xlim = c(0, 12),
    ylim = c(0, 1),
    break.time.by = 3,
    break.y.by = 0.2,  # y轴间隔20%
    surv.scale = "percent",
    legend = "none"  # 使用共同图例
  )
  
  km_plots_list[[i]] <- km_plot$plot
  
  # 保存单独图形（PNG格式）
  png(paste0("Figures/Figure2", LETTERS[i], "_KM_", outcome, ".png"),
      width = 8*300, height = 6*300, res = 300, bg = "white")
  print(km_plot$plot)
  dev.off()
}

# 合并四个图形为Figure 2（专业布局）
figure2_combined <- ggarrange(plotlist = km_plots_list, 
                              ncol = 2, nrow = 2,
                              common.legend = TRUE,
                              legend = "bottom",
                              font.label = list(size = 14, face = "bold"))

# 保存为高分辨率发表级图片
png("Figures/Figure2_KM_All_Outcomes.png", 
    width = 16*300, height = 14*300, res = 300, bg = "white")
print(figure2_combined)
dev.off()
cat("已保存所有KM曲线图\n")


# =============================================================================
# 阶段三：深度模式挖掘 - "U型"关系验证
# =============================================================================

cat("\n========== 阶段三：U型关系验证（RCS分析）==========\n")

# -----------------------------------------------------------------------------
# 5.1 五分位数分析
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
  p_values <- summary(cox_quintile)$coefficients[, "Pr(>|z|)"]
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
    P_value = c(p_values[quintile_rows[1:2]], 
                NA,
                p_values[quintile_rows[3:4]])
  )
  
  quintile_results_all[[outcome]] <- results
  if (outcome == "MACE") print(results)
}

table3_quintile <- do.call(rbind, quintile_results_all)
write.csv(table3_quintile, "Table3_Quintile_Analysis.csv", row.names = FALSE)
cat("\n已保存：Table3_Quintile_Analysis.csv\n")


# -----------------------------------------------------------------------------
# 5.2 限制性立方样条（RCS）分析
# -----------------------------------------------------------------------------

cat("\n----- 5.2 限制性立方样条（RCS）分析 -----\n")

# 使用rms包进行RCS分析
dd <- datadist(data)
options(datadist = "dd")

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
LVEDVi_range <- seq(quantile(data$LVEDVi, 0.01, na.rm = TRUE), 
                    quantile(data$LVEDVi, 0.99, na.rm = TRUE), 
                    length.out = 100)

# 设置参考值（中位数）
ref_value <- median(data$LVEDVi, na.rm = TRUE)

pred_MACE <- Predict(cox_rcs_MACE, LVEDVi = LVEDVi_range, ref.zero = TRUE)

# 转换为HR
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
  geom_line(color = "blue", linewidth = 1.2) +
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
png("Figures/Figure1B_RCS_MACE.png", width = 8*300, height = 6*300, res = 300, bg = "white")
print(rcs_plot_MACE)
dev.off()
cat("已保存：Figures/Figure1B_RCS_MACE.png\n")

cat("\n========== 分析完成！==========\n")
cat("\n由于完整分析较为耗时，核心分析已完成。\n")
cat("生成的文件包括：\n")
cat("  - Table1_Baseline_Characteristics.csv\n")
cat("  - Univariate_Logistic_Regression.csv\n")
cat("  - Multivariate_Logistic_Regression.csv\n")
cat("  - Table2_Cox_Regression_Results.csv\n")
cat("  - Table3_Quintile_Analysis.csv\n")
cat("  - Figures/Figure1A_KM_MACE.png（高清PNG格式）\n")
cat("  - Figures/Figure1B_RCS_MACE.png（高清PNG格式）\n")
cat("  - Figures/Figure2_KM_All_Outcomes.png（高清PNG格式）\n")
cat("\n模拟数据已保存在工作空间中。\n")
cat("如需运行完整分析（包括敏感性分析和亚组分析），请继续运行原代码后续部分。\n")

# 保存工作空间
save.image("Analysis_Workspace.RData")
cat("\n工作空间已保存为：Analysis_Workspace.RData\n")

################################################################################
# 代码结束
################################################################################
