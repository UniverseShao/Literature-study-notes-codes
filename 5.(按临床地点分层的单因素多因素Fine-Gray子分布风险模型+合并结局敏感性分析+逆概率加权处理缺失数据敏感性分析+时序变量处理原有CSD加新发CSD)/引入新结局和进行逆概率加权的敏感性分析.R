library(cmprsk)
library(survival)
library(riskRegression)
library(cmprsk)
library(prodlim) # FGR函数需要这个包
# 第一项敏感性分析：合并结局事件------------------------------------------------
# 1.第一步：创建模拟数据集------------------------------------------------------
# 设置随机数种子，以保证每次运行的结果都一样
set.seed(123)

# 定义患者数量
n_patients <- 6359 # 这是研究中基线无CSD的患者数量

# 创建模拟数据集
base_covariates <- data.frame(
  # 1. 生存时间和事件状态 (这是竞争风险模型的核心)
  # time: 从入组到发生事件或随访结束的时间（单位：月）
  time = runif(n_patients, 12, 48), 
  
  # status: 结局状态
  # 0 = 删失 (Censored): 随访结束时，既没有死亡，也没有发生CSD
  # 1 = 事件 (Event of Interest): 发生了新发CSD (这是我们关心的结局)
  # 2 = 竞争风险 (Competing Risk): 因其他原因死亡 (死亡阻止了CSD的发生)
  status = sample(c(0, 1, 2), n_patients, replace = TRUE, prob = c(0.70, 0.18, 0.12)),
  
  # 2. 核心自变量：治疗分组
  treatment = factor(sample(c("Standard", "Intensive"), n_patients, replace = TRUE)),
  
  # 3. 分层变量：临床中心
  clinical_site = factor(sample(paste0("Site_", 1:15), n_patients, replace = TRUE)),
  
  # 4. 需要校正的协变量 (按照方法学描述创建)
  age = rnorm(n_patients, mean = 65, sd = 5),
  sex = factor(sample(c("Male", "Female"), n_patients, replace = TRUE, prob = c(0.45, 0.55))),
  bmi = rnorm(n_patients, mean = 25.5, sd = 3),
  sbp_baseline = rnorm(n_patients, mean = 146, sd = 17),
  smoking = factor(sample(c("Yes", "No"), n_patients, replace = TRUE, prob = c(0.15, 0.85))),
  diabetes = factor(sample(c("Yes", "No"), n_patients, replace = TRUE, prob = c(0.20, 0.80))),
  hyperlipidemia = factor(sample(c("Yes", "No"), n_patients, replace = TRUE, prob = c(0.35, 0.65))),
  cvd_history = factor(sample(c("Yes", "No"), n_patients, replace = TRUE, prob = c(0.06, 0.94))),
  aspirin_use = factor(sample(c("Yes", "No"), n_patients, replace = TRUE, prob = c(0.10, 0.90))),
  statin_use = factor(sample(c("Yes", "No"), n_patients, replace = TRUE, prob = c(0.20, 0.80))),
  arb_use = factor(sample(c("Yes", "No"), n_patients, replace = TRUE, prob = c(0.44, 0.56))),
  egfr = rnorm(n_patients, mean = 109, sd = 24),
  total_chol = rnorm(n_patients, mean = 4.9, sd = 1.1),
  hdl_c = rnorm(n_patients, mean = 1.27, sd = 0.3)
)

# 1. 模拟"真实"的事件发生时间 (假设服从指数分布)
# 我们让强化治疗显著降低CSD风险 (HR < 1)，增加治疗效果差异
# 让年龄增加死亡和CSD风险 (HR > 1)
hazard_csd <- 0.08 * exp(-0.4 * (base_covariates$treatment == "Intensive") + 0.03 * (base_covariates$age - 65))
hazard_death <- 0.05 * exp(-0.2 * (base_covariates$treatment == "Intensive") + 0.04 * (base_covariates$age - 65))
time_to_csd <- rexp(n_patients, rate = hazard_csd)
time_to_death <- rexp(n_patients, rate = hazard_death)

# 2. 模拟“删失”时间
# 研究结束时间 (管理性删失)，假设最长随访48个月
study_end_time <- rep(48, n_patients)
# 真实失访时间 (Loss to follow-up)，假设一部分人会提前失访
time_to_loss <- runif(n_patients, 1, 48)
# 假设有10%的人会发生真实失访
is_lost <- rbinom(n_patients, 1, 0.1) == 1
# 结合两种删失情况
censoring_time <- ifelse(is_lost, time_to_loss, study_end_time)


# --- 确定最终的观察时间和结局状态 ---
# 最终观察时间是三个时间（CSD发生、死亡发生、被删失）中最早发生的那一个
final_time <- pmin(time_to_csd, time_to_death, censoring_time)

# 确定最终的结局状态
# status = 0 (删失), 1 (CSD), 2 (死亡)
final_status <- ifelse(final_time == censoring_time, 0, 
                       ifelse(final_time == time_to_csd, 1, 2))


# --- 组合成最终的数据集 ---
step_data_realistic <- cbind(base_covariates, time = final_time, status = final_status)
step_data_realistic$is_lost <- is_lost # 加入是否真实失访的标志

# 查看新数据集的结局分布
print("--- 更真实的数据集结局分布 ---")
print(table(status = step_data_realistic$status, is_lost = step_data_realistic$is_lost))

# 查看数据前几行，确保创建成功
head(step_data_realistic)

# --- 1. 在数据集中模拟"新植入起搏器"事件 ---
# 假设在随访期间，有大约3%的患者植入了新起搏器，增加复合事件发生率
# 我们只在那些原始状态为"删失"(status=0)的患者中模拟，因为如果已经发生CSD或死亡，就不会再观察到单纯的起搏器事件
n_censored <- sum(step_data_realistic$status == 0)
step_data_realistic$pacemaker_event <- 0
step_data_realistic$pacemaker_event[step_data_realistic$status == 0] <- sample(c(1, 0), 
                                                           size = n_censored, 
                                                           replace = TRUE, 
                                                           prob = c(0.03, 0.97))

# --- 2. 创建新的复合结局状态变量 status_composite ---
step_data_realistic$status_composite <- step_data_realistic$status

# 将那些发生了“起搏器事件”的删失患者，其结局状态从0(删失)改为1(复合事件)
step_data_realistic$status_composite[step_data_realistic$status == 0 & step_data_realistic$pacemaker_event == 1] <- 1

# 我们可以通过表格来验证我们的操作是否正确
print("--- 原始结局 vs. 复合结局对比 ---")
print(table(original_status = step_data_realistic$status, composite_status = step_data_realistic$status_composite))

# 2.第二步运行使用复合结局的多变量Fine-Gray模型---------------------------------
# 模型结构与我们之前的主分析完全相同，只是结局变量变了
sens_model_1 <- FGR(
  prodlim::Hist(time, status_composite) ~ treatment + age + sex + bmi + sbp_baseline + 
    smoking + diabetes + hyperlipidemia + cvd_history + 
    aspirin_use + statin_use + arb_use + egfr + 
    total_chol + hdl_c + strata(clinical_site),
  data = step_data_realistic,
  cause = 1
)

print("--- 敏感性分析 1 (复合结局) 的结果 ---")
summary(sens_model_1)

# 第二项敏感性分析：逆概率加权 (IPW)--------------------------------------------

# 加载必要的包
library(WeightIt)  # 用于计算权重
library(survival)  # 用于生存分析
library(etm)       # 用于加权的累积发病率函数

# 1. 创建倾向性评分模型来计算逆概率权重--------------------------------------------
# 我们需要模拟一些缺失数据的情况，因为IPW通常用于处理缺失数据或选择偏倚

# 1.1 模拟缺失数据机制
# 假设某些协变量存在缺失，缺失概率与其他变量相关
set.seed(456)

# 创建缺失指示变量 - 假设BMI数据有缺失
# 缺失概率与年龄、性别、治疗组相关（这是一个现实的缺失机制）
missing_prob <- plogis(-2 + 0.02 * step_data_realistic$age + 
                      0.3 * (step_data_realistic$sex == "Male") + 
                      0.2 * (step_data_realistic$treatment == "Intensive"))

# 生成缺失指示变量
step_data_realistic$bmi_missing <- rbinom(n_patients, 1, missing_prob)

# 创建包含缺失值的BMI变量
step_data_realistic$bmi_with_missing <- step_data_realistic$bmi
step_data_realistic$bmi_with_missing[step_data_realistic$bmi_missing == 1] <- NA

# 查看缺失情况
cat("BMI缺失情况:\n")
cat("缺失数量:", sum(is.na(step_data_realistic$bmi_with_missing)), "\n")
cat("缺失比例:", round(mean(is.na(step_data_realistic$bmi_with_missing)) * 100, 2), "%\n")

# 1.2 计算逆概率权重
# 我们使用logistic回归来预测观察到完整数据的概率
# 然后计算逆概率权重

# 创建完整观察指示变量（1=完整观察，0=有缺失）
step_data_realistic$complete_obs <- as.numeric(!is.na(step_data_realistic$bmi_with_missing))

# 拟合倾向性评分模型（预测完整观察的概率）
ps_model <- glm(complete_obs ~ treatment + age + sex + sbp_baseline + 
                smoking + diabetes + hyperlipidemia + cvd_history + 
                aspirin_use + statin_use + arb_use + egfr + 
                total_chol + hdl_c + clinical_site,
                data = step_data_realistic,
                family = binomial())

# 计算倾向性评分（完整观察的预测概率）
step_data_realistic$ps <- predict(ps_model, type = "response")

# 计算逆概率权重
# 对于完整观察的个体，权重 = 1/P(完整观察|协变量)
# 对于有缺失的个体，我们在这个分析中将其排除
step_data_realistic$ipw <- ifelse(step_data_realistic$complete_obs == 1, 
                                 1/step_data_realistic$ps, 
                                 0)

# 稳定化权重（可选，用于减少极端权重的影响）
mean_ps <- mean(step_data_realistic$ps)
step_data_realistic$stabilized_ipw <- ifelse(step_data_realistic$complete_obs == 1,
                                           mean_ps/step_data_realistic$ps,
                                           0)

# 查看权重分布
cat("\n逆概率权重统计:\n")
complete_data <- step_data_realistic[step_data_realistic$complete_obs == 1, ]
cat("权重范围:", round(range(complete_data$ipw), 3), "\n")
cat("权重均值:", round(mean(complete_data$ipw), 3), "\n")
cat("权重中位数:", round(median(complete_data$ipw), 3), "\n")

# 2. 使用加权数据进行竞争风险分析----------------------------------------------

# 2.1 方法1：使用加权的累积发病率估计
print("=== 方法1：使用加权的累积发病率估计 ===")

# 由于数据复制方法可能导致内存问题，我们改用更直接的方法
# 使用survey包的加权生存分析功能

library(survey)

# 准备完整观察的数据
complete_data <- step_data_realistic[step_data_realistic$complete_obs == 1, ]

# 创建加权设计对象
weighted_design <- svydesign(ids = ~1, weights = ~stabilized_ipw, data = complete_data)

# 按治疗组分别计算加权的累积发病率
print("计算标准治疗组的加权累积发病率...")
standard_data <- complete_data[complete_data$treatment == "Standard", ]
standard_design <- svydesign(ids = ~1, weights = ~stabilized_ipw, data = standard_data)

print("计算强化治疗组的加权累积发病率...")
intensive_data <- complete_data[complete_data$treatment == "Intensive", ]
intensive_design <- svydesign(ids = ~1, weights = ~stabilized_ipw, data = intensive_data)

# 使用加权的Kaplan-Meier方法计算CSD的累积发病率
# 创建CSD事件指示变量
standard_data$csd_event <- ifelse(standard_data$status == 1, 1, 0)
intensive_data$csd_event <- ifelse(intensive_data$status == 1, 1, 0)

# 首先检查数据的时间范围
print(paste("数据时间范围: ", round(min(complete_data$time), 2), " 到 ", round(max(complete_data$time), 2)))
print(paste("标准治疗组样本数: ", nrow(standard_data)))
print(paste("强化治疗组样本数: ", nrow(intensive_data)))

# 使用实际的随访时间范围来计算累积发病率
max_time <- max(complete_data$time)
time_point <- min(max_time, 60)  # 使用60个月或最大随访时间中的较小值

print(paste("使用的时间点: ", round(time_point, 2), " 个月"))

# 标准治疗组累积发病率 - 修正计算方法
# 计算在time_point时间内发生CSD事件的加权比例
standard_events_by_timepoint <- standard_data$time <= time_point & standard_data$csd_event == 1
standard_total_weight <- sum(standard_data$stabilized_ipw)
standard_event_weight <- sum(standard_data$stabilized_ipw[standard_events_by_timepoint])
standard_cif <- standard_event_weight / standard_total_weight

# 强化治疗组累积发病率 - 修正计算方法
intensive_events_by_timepoint <- intensive_data$time <= time_point & intensive_data$csd_event == 1
intensive_total_weight <- sum(intensive_data$stabilized_ipw)
intensive_event_weight <- sum(intensive_data$stabilized_ipw[intensive_events_by_timepoint])
intensive_cif <- intensive_event_weight / intensive_total_weight

print(paste("标准治疗组", round(time_point/12, 1), "年CSD累积发病率（加权）:", round(standard_cif * 100, 2), "%"))
print(paste("强化治疗组", round(time_point/12, 1), "年CSD累积发病率（加权）:", round(intensive_cif * 100, 2), "%"))
print(paste("风险差异（强化 vs 标准）:", round((intensive_cif - standard_cif) * 100, 2), "%"))

print("方法1完成")

# 2.2 方法2：使用加权的Fine-Gray模型（简化方法）
# 由于FGR没有直接的权重参数，我们使用一个简化的方法
# 通过对权重进行分层分析来近似加权效应

print("=== 方法2：使用分层Fine-Gray模型 ===")

# 将权重分为几个层次来进行分层分析
complete_data$weight_strata <- cut(complete_data$stabilized_ipw, 
                                  breaks = quantile(complete_data$stabilized_ipw, 
                                                   probs = c(0, 0.25, 0.5, 0.75, 1.0)),
                                  include.lowest = TRUE,
                                  labels = c("Low", "Medium-Low", "Medium-High", "High"))

print("权重分层情况:")
print(table(complete_data$weight_strata))

# 为每个权重层拟合Fine-Gray模型
strata_results <- list()
for(stratum in levels(complete_data$weight_strata)) {
  cat("\n拟合权重层:", stratum, "\n")
  stratum_data <- complete_data[complete_data$weight_strata == stratum, ]
  
  if(nrow(stratum_data) > 50) {  # 确保有足够的样本量
    tryCatch({
      strata_model <- FGR(
        prodlim::Hist(time, status) ~ treatment + age + sex + bmi + sbp_baseline + 
          smoking + diabetes + hyperlipidemia + cvd_history + 
          aspirin_use + statin_use + arb_use + egfr + 
          total_chol + hdl_c + strata(clinical_site),
        data = stratum_data,
        cause = 1
      )
      strata_results[[stratum]] <- strata_model
      cat("模型拟合成功，样本量:", nrow(stratum_data), "\n")
    }, error = function(e) {
      cat("模型拟合失败:", e$message, "\n")
    })
  } else {
    cat("样本量不足，跳过该层\n")
  }
}

# 计算加权平均的治疗效应
if(length(strata_results) > 0) {
  treatment_hrs <- sapply(strata_results, function(model) {
    if(!is.null(model) && !is.null(model$coef) && length(model$coef) > 0) {
      tryCatch({
        exp(model$coef[1])  # 假设treatment是第一个系数
      }, error = function(e) {
        cat("计算HR时出错:", e$message, "\n")
        return(NA)
      })
    } else {
      return(NA)
    }
  })
  
  # 计算各层的样本量权重
  strata_weights <- sapply(names(strata_results), function(stratum) {
    sum(complete_data$weight_strata == stratum)
  })
  
  # 加权平均HR
  valid_hrs <- !is.na(treatment_hrs)
  if(sum(valid_hrs) > 0) {
    weighted_hr <- weighted.mean(treatment_hrs[valid_hrs], 
                                strata_weights[valid_hrs])
    cat("\n分层加权平均治疗效应 (HR):", round(weighted_hr, 3), "\n")
  } else {
    cat("\n无法计算加权平均HR，所有分层模型都未能提供有效系数\n")
  }
}

# 3. 输出IPW敏感性分析结果----------------------------------------------------

cat("\n=== 敏感性分析 2: 逆概率加权 (IPW) 结果 ===\n")

cat("\n1. 倾向性评分模型摘要:\n")
cat("AIC:", round(AIC(ps_model), 2), "\n")
cat("完整观察预测准确性:", round(mean((ps_model$fitted.values > 0.5) == step_data_realistic$complete_obs), 3), "\n")

cat("\n2. 加权累积发病率结果:\n")
cat("标准治疗组", round(time_point, 1), "年CSD累积发病率（加权）:", round(standard_cif * 100, 2), "%\n")
cat("强化治疗组", round(time_point, 1), "年CSD累积发病率（加权）:", round(intensive_cif * 100, 2), "%\n")
cat("风险差异（强化 vs 标准）:", round((intensive_cif - standard_cif) * 100, 2), "%\n")

# 如果有分层模型结果，也输出
if(exists("strata_results") && length(strata_results) > 0) {
  cat("\n3. 分层Fine-Gray模型结果:\n")
  for(stratum in names(strata_results)) {
    if(!is.null(strata_results[[stratum]])) {
      model <- strata_results[[stratum]]
      tryCatch({
        if(!is.null(model$coef) && length(model$coef) > 0) {
          treatment_hr <- exp(model$coef[1])
          cat("权重层", stratum, "- 治疗效应HR:", round(treatment_hr, 3), "\n")
        } else {
          cat("权重层", stratum, "- 无有效系数\n")
        }
      }, error = function(e) {
        cat("权重层", stratum, "- 计算HR时出错:", e$message, "\n")
      })
    }
  }
  
  if(exists("weighted_hr")) {
    cat("分层加权平均治疗效应 (HR):", round(weighted_hr, 3), "\n")
  }
}

# 4. 与主分析结果比较--------------------------------------------------------

cat("\n=== 敏感性分析结果比较 ===\n")
if(exists("sens_model_1")) {
  tryCatch({
    if(!is.null(sens_model_1$coef) && length(sens_model_1$coef) > 0) {
      cat("1. 复合结局敏感性分析 - 治疗效应HR:", round(exp(sens_model_1$coef[1]), 3), "\n")
    } else {
      cat("1. 复合结局敏感性分析 - 无有效系数\n")
    }
  }, error = function(e) {
    cat("1. 复合结局敏感性分析 - 计算HR时出错:", e$message, "\n")
  })
} else {
  cat("1. 复合结局敏感性分析 - 模型不存在\n")
}

cat("2. IPW敏感性分析结果:\n")
cat("   - 累积发病率差异:", round((intensive_cif - standard_cif) * 100, 2), "%\n")
if(exists("weighted_hr")) {
  cat("   - 分层加权HR:", round(weighted_hr, 3), "\n")
} else {
  cat("   - 分层加权HR: 无法计算\n")
}

cat("\n=== 敏感性分析解释 ===\n")
cat("1. 复合结局敏感性分析：通过将'新发CSD'和'新植入起搏器'合并为复合结局，\n")
cat("   评估结局定义变化对结果的影响。\n")
cat("2. IPW敏感性分析：通过逆概率加权方法处理潜在的缺失数据偏倚，\n")
cat("   评估缺失数据对结果稳健性的影响。\n")
cat("3. 如果两种敏感性分析的结果与主分析一致，说明研究结论具有较好的稳健性。\n")

print("=== 敏感性分析完成 ===")

