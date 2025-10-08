# --- 1. 安装和加载必要的R包 ---
# install.packages(c("survival", "survey", "dplyr", "EValue"))
library(survival)
library(survey)
library(dplyr)
library(EValue)

# --- 2. 创建一个逼真的模拟数据集 ---
set.seed(42) # 保证结果可复现
n <- 2000 # 模拟参与者数量

sim_data <- data.frame(
  id = 1:n,
  # 结局变量
  time_years = runif(n, 0.1, 20), # 随访时间 (0.1到20年)
  cvd_death = rbinom(n, 1, 0.15), # 结局事件 (1=心血管死亡, 0=删失)
  
  # 主要预测变量
  TyG_WWI = rnorm(n, 10, 2),
  
  # 协变量
  age = sample(18:80, n, replace = TRUE),
  sex = factor(sample(c("Male", "Female"), n, replace = TRUE)),
  smoking = factor(sample(c("Never", "Former", "Current", NA), n, replace = TRUE, prob = c(0.5, 0.25, 0.2, 0.05))),
  baseline_cvd = rbinom(n, 1, 0.1), # 基线是否有CVD (1=有, 0=无)
  
  # 用于特定分析的变量
  survey_year = sample(1999:2018, n, replace = TRUE),
  
  # NHANES 复杂抽样设计所需的变量 (模拟)
  psu = sample(1:100, n, replace = TRUE), # 主抽样单位
  strata = sample(1:50, n, replace = TRUE), # 分层
  weight = runif(n, 500, 50000) # 抽样权重
)

# --- 3. 创建NHANES调查设计对象 ---
# 这是所有加权分析的基础
nhanes_design <- svydesign(
  id = ~psu,
  strata = ~strata,
  weights = ~weight,
  data = sim_data,
  nest = TRUE
)

# --- 4. 定义我们的主分析模型 (模型II) ---
# 后续所有分析都将基于这个模型的变体
main_formula <- as.formula(
  "Surv(time_years, cvd_death) ~ TyG_WWI + age + sex + smoking"
)

# 运行一次主模型作为基准
cat("--- 基准主模型分析结果 ---\n")
baseline_model <- svycoxph(main_formula, design = nhanes_design)
print(summary(baseline_model))

# 1. 亚组分析 (Subgroup Analysis) ---------------------------------------------------
# - **目的**: 检验核心结论 (TyG-WWI与心血管死亡风险相关) 是否在不同人群亚组中保持一致。
# - **测试的问题**: “这个结论是不是只对特定人群(如男性、吸烟者)有效？”
# - **重要性**: 一个稳健的预测指标应该具有普适性。如果在各个亚组中结论都成立，说明其价值更大。
cat("\n\n--- 1. 亚组分析 ---\n")

# 定义要进行亚组分析的变量
subgroup_vars <- c("sex", "smoking")

for (variable in subgroup_vars) {
  # 获取该变量的所有水平 (除去NA)
  levels <- unique(na.omit(sim_data[[variable]]))

  cat(paste0("\n--- 按 '", variable, "' 进行亚组分析 ---\n"))

  for (level in levels) {
    cat(paste0("  * 亚组: ", variable, " = ", level, "\n"))

    # 检查该亚组的样本量和结局事件数量
    subgroup_data <- subset(sim_data, get(variable) == level)
    n_subgroup <- nrow(subgroup_data)
    n_events <- sum(subgroup_data$cvd_death, na.rm = TRUE)

    cat(paste0("    样本量: ", n_subgroup, ", 结局事件数: ", n_events, "\n"))

    # 如果样本量太小或没有足够的事件，跳过分析
    if (n_subgroup < 50 || n_events < 5) {
      cat("    警告: 样本量太小或事件数不足，跳过此亚组分析。\n")
      next
    }

    # 创建该亚组的调查设计对象
    subgroup_design <- subset(nhanes_design, get(variable) == level)

    # 添加异常处理
    tryCatch({
      # 在亚组上运行模型
      subgroup_model <- svycoxph(main_formula, design = subgroup_design)
      print(summary(subgroup_model))
    }, error = function(e) {
      cat(paste0("    错误: ", e$message, "\n"))
      cat("    跳过此亚组分析。\n")
    })
  }
}

#  2. 多重比较校正 (Bonferroni Correction) -----------------------------------------------------
# - **目的**: 当进行多次统计检验时，控制假阳性错误的概率。
# - **测试的问题**: “做了这么多亚组分析，会不会某个结果只是‘碰巧’看起来显著？”
# - **重要性**: 校正P值的阈值，使得标准更严格，以确保任何声称的显著结果都不是由多次重复检验引起的偶然发现。
cat("\n\n--- 2. 多重比较校正 ---\n")
num_subgroups <- 2 + 3 # 2个性别亚组 + 3个吸烟亚组
alpha <- 0.05
bonferroni_p_value <- alpha / num_subgroups

cat(paste("进行了", num_subgroups, "次亚组比较。\n"))
cat(paste("根据Bonferroni校正，新的显著性阈值为:", round(bonferroni_p_value, 4), "\n"))
cat("请检查上面亚组分析的P值，看其是否低于这个新阈值。\n")


#  3. 分时期分析 (Stratified Analyses by Time Period) --------------------------------
# - **目的**: 检验结论是否随时间推移而保持稳定。
# - **测试的问题**: “研究横跨近20年，医疗水平和生活方式都在变。结论会不会只是某个特定时代的产物？”
# - **重要性**: 如果关联在不同时期(如2008年前后)都存在，说明其预测能力是持久的，不受外界医疗环境变化的大的影响。
cat("\n\n--- 3. 分时期分析 ---\n")

# 创建时期变量
design_with_period <- update(
  nhanes_design, 
  time_period = factor(ifelse(survey_year <= 2008, "1999-2008", "2009-2018"))
)

# a. 1999-2008 时期
cat("--- 时期: 1999-2008 ---\n")
design_period1 <- subset(design_with_period, time_period == "1999-2008")
model_period1 <- svycoxph(main_formula, design = design_period1)
print(summary(model_period1))

# b. 2009-2018 时期
cat("\n--- 时期: 2009-2018 ---\n")
design_period2 <- subset(design_with_period, time_period == "2009-2018")
model_period2 <- svycoxph(main_formula, design = design_period2)
print(summary(model_period2))

# 4. 界标分析 (Landmark Analysis) ------------------------------------------------
# - **目的**: 评估指标对于远期风险的预测能力，区分短期和长期预测价值。
# - **测试的问题**: “这个指标是只能预测那些‘马上就不行了’的病人，还是能预测一个健康人10年后的远期风险？”
# - **重要性**: 通过只分析在某个时间点（“界标”，如10年）仍然存活的个体，可以检验指标是否能预测“幸存者”未来的长期风险，这对于早期预警至关重要。
cat("\n\n--- 4. 界标分析 ---\n")
landmark_time_years <- 10 # 设定10年为界标点

# 1. 筛选出在界标点时仍然存活的参与者
landmark_data <- sim_data %>%
  filter(time_years > landmark_time_years)

# 2. 更新这些幸存者的随访时间
landmark_data <- landmark_data %>%
  mutate(new_time = time_years - landmark_time_years)

# 3. 为这个新数据集创建调查设计对象
landmark_design <- svydesign(
  id = ~psu, strata = ~strata, weights = ~weight,
  data = landmark_data, nest = TRUE
)

# 4. 运行模型，注意结局变量的生存时间已更新
landmark_formula <- as.formula(
  "Surv(new_time, cvd_death) ~ TyG_WWI + age + sex + smoking"
)
landmark_model <- svycoxph(landmark_formula, design = landmark_design)
print(summary(landmark_model))

# 5. 排除早期死亡事件 ----------------------------------------------------
# - **目的**: 排除“因果倒置”的可能性。
# - **测试的问题**: “会不会是因为某些未被发现的重病，既导致了指标异常，也导致了短期死亡，从而造成了指标和死亡的虚假关联？”
# - **重要性**: 通过移除随访早期（如2年内）就死亡的个体，可以更有力地证明指标是一个“因”（风险因素），而不是一个“果”（已存在疾病的标志物）。
cat("\n\n--- 5. 排除早期死亡事件 (2年内) ---\n")

# 筛选数据：保留所有删失的，以及死亡时间>2年的
data_no_early_death <- sim_data %>%
  filter(!(cvd_death == 1 & time_years <= 2))

# 基于筛选后的数据创建新的调查设计
design_no_early_death <- svydesign(
  id = ~psu, strata = ~strata, weights = ~weight,
  data = data_no_early_death, nest = TRUE
)
model_no_early_death <- svycoxph(main_formula, design = design_no_early_death)
print(summary(model_no_early_death))

# 6. 完整病例分析 (Complete Case Analysis) ------------------------------------------
# - **目的**: 检验结果对于缺失数据处理方法的敏感性。
# - **测试的问题**: “目前对缺失值的处理方式是否影响了最终结论？如果用最简单粗暴的方法（直接删除）结果会怎样？”
# - **重要性**: 如果采用不同的缺失值处理方法（这里是仅保留信息完整的参与者），结论依然不变，说明结果是稳健的，不是特定统计技巧的产物。
cat("\n\n--- 6. 完整病例分析 ---\n")

# 移除任何含有NA的行
data_complete_case <- na.omit(sim_data)
cat(paste("原始数据:", nrow(sim_data), "行\n"))
cat(paste("完整病例数据:", nrow(data_complete_case), "行\n"))

# 基于完整病例数据创建调查设计
design_complete_case <- svydesign(
  id = ~psu, strata = ~strata, weights = ~weight,
  data = data_complete_case, nest = TRUE
)
model_complete_case <- svycoxph(main_formula, design = design_complete_case)
print(summary(model_complete_case))

# 7. 排除基线患有心血管疾病的患者 --------------------------------------------------
# - **目的**: 检验指标在“健康”人群中的预测价值，即一级预防价值。
# - **测试的问题**: “对于已经有病的人，预测复发相对容易。这个指标能否在一大群健康人里，找出未来的高危分子？”
# - **重要性**: 一级预防（在健康人中预防首次发病）是公共卫生的重点。如果指标在一级预防人群中依然有效，说明其有潜力用于大规模人群的早期筛查。
cat("\n\n--- 7. 排除基线CVD患者 (一级预防人群) ---\n")

# 筛选出基线没有CVD的人
design_primary_prevention <- subset(nhanes_design, baseline_cvd == 0)

model_primary_prevention <- svycoxph(main_formula, design = design_primary_prevention)
print(summary(model_primary_prevention))

# 8. E值计算 (E-value Calculation) -------------------------------------------------
# - **目的**: 量化未测量的混杂因素对结论的潜在影响。
# - **测试的问题**: “虽然考虑了很多已知干扰因素，但万一存在一个未知的超级干扰因素（如基因、环境）同时影响指标和死亡，该怎么办？”
# - **重要性**: E值评估了一个未知的混杂因素需要有多强大（同时与暴露和结局的关联强度），才能完全解释掉我们观察到的关联。E值越大，结论就越稳健，因为找到如此强大的未知混杂因素的可能性越小。
cat("\n\n--- 8. E值计算 ---\n")

# 检查 "TyG_WWI" 是否在模型系数中
if ("TyG_WWI" %in% names(coef(baseline_model))) {
  # 1. 提取HR
  hr <- exp(coef(baseline_model)["TyG_WWI"])
  # 2. 提取95%置信区间
  ci <- exp(confint(baseline_model)["TyG_WWI", ])

  # 确保提取成功
  if (!is.na(hr) && !any(is.na(ci))) {
    cat(paste0("TyG_WWI的风险比 (HR): ", round(hr, 2), "\n"))
    cat(paste0("95% 置信区间: [", round(ci[1], 2), ", ", round(ci[2], 2), "]\n\n"))

    # 3. 计算E值
    if (hr > 1) {
      # 添加 rare = FALSE 参数，因为结局事件发生率(15%)不被认为是罕见的
      e_value_result <- evalues.HR(est = hr, lo = ci[1], hi = ci[2], rare = FALSE)
      cat("--- E值结果 ---\n")
      print(e_value_result)
    } else {
      cat("风险比(HR)不大于1，计算E值没有意义。\n")
    }
  } else {
    cat("无法从模型中提取有效的风险比或置信区间。\n")
  }
} else {
  cat("模型系数中未找到 'TyG_WWI'，无法计算E值。\n")
}
