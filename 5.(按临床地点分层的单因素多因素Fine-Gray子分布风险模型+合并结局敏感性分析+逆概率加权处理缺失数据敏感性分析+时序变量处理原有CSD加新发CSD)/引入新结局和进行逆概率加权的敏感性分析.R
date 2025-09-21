library(cmprsk)
library(survival)
library(riskRegression)
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
  status = sample(c(0, 1, 2), n_patients, replace = TRUE, prob = c(0.85, 0.06, 0.09)),
  
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

# 1. 模拟“真实”的事件发生时间 (假设服从指数分布)
# 我们让强化治疗稍微降低一点死亡风险 (HR < 1)
# 让年龄增加死亡和CSD风险 (HR > 1)
hazard_csd <- 0.02 * exp(0.1 * (base_covariates$treatment == "Intensive") + 0.05 * (base_covariates$age - 65))
hazard_death <- 0.03 * exp(-0.1 * (base_covariates$treatment == "Intensive") + 0.06 * (base_covariates$age - 65))
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

# --- 1. 在数据集中模拟“新植入起搏器”事件 ---
# 假设在随访期间，有大约1%的患者植入了新起搏器
# 我们只在那些原始状态为"删失"(status=0)的患者中模拟，因为如果已经发生CSD或死亡，就不会再观察到单纯的起搏器事件
n_censored <- sum(step_data_realistic$status == 0)
step_data_realistic$pacemaker_event <- 0
step_data_realistic$pacemaker_event[step_data_realistic$status == 0] <- sample(c(1, 0), 
                                                           size = n_censored, 
                                                           replace = TRUE, 
                                                           prob = c(0.01, 0.99))

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
