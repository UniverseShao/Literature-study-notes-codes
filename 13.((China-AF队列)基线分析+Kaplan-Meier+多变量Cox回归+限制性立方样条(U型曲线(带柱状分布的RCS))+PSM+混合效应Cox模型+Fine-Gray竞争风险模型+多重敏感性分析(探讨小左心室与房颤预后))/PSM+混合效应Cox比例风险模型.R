# -----------------------------------------------------------------------------
# 步骤 0: 加载必要的R包
# -----------------------------------------------------------------------------
library(tidyverse) # 用于数据处理和可视化
library(MatchIt)   # 用于倾向性评分匹配的核心包
library(survival)  # 用于生存分析 (Surv对象和标准Cox模型)
library(coxme)     # 用于混合效应Cox模型
library(tableone)  # 用于轻松创建基线特征表 ("Table 1")

# -----------------------------------------------------------------------------
# 步骤 1: 创建一个模拟的、不平衡的数据集
# -----------------------------------------------------------------------------
# 设定一个种子以便结果可复现
set.seed(123)

# 创建一个包含10000名患者的数据集
# 其中只有少数(n=500)是"小左心室"患者
n_total <- 10000
n_small_lv <- 500

patient_data <- tibble(
  id = 1:n_total,
  # 创建分组变量：'Small LV' vs 'Normal LV'
  lv_size = factor(c(rep("Small LV", n_small_lv), rep("Normal LV", n_total - n_small_lv))),
  
  # 创建一些混杂变量。我们故意让'Small LV'组的患者健康状况更差
  # 比如年龄更大，有心衰的比例更高
  age = ifelse(lv_size == "Small LV", rnorm(n_total, mean = 70, sd = 8), rnorm(n_total, mean = 62, sd = 10)),
  sex = factor(sample(c("Male", "Female"), n_total, replace = TRUE, prob = c(0.6, 0.4))),
  sbp = ifelse(lv_size == "Small LV", rnorm(n_total, mean = 145, sd = 15), rnorm(n_total, mean = 135, sd = 15)),
  has_hf = ifelse(lv_size == "Small LV", rbinom(n_total, 1, 0.40), rbinom(n_total, 1, 0.20)), # 心衰比例
  has_diabetes = ifelse(lv_size == "Small LV", rbinom(n_total, 1, 0.35), rbinom(n_total, 1, 0.25)), # 糖尿病比例
  
  # 根据患者情况模拟生存时间和结局
  # 让"Small LV"本身以及其他坏因素增加风险
  true_hazard = 0.5 * (lv_size == "Small LV") + 0.04 * (age - 65) + 0.3 * has_hf,
  time_to_event = rexp(n_total, rate = 0.1 * exp(true_hazard)),
  
  # 模拟一个随访删失时间 (有些人可能中途失访或研究结束时仍存活)
  censoring_time = runif(n_total, 1, 10),
  
  # 最终观察到的时间和事件状态
  observed_time = pmin(time_to_event, censoring_time),
  event_status = as.numeric(time_to_event <= censoring_time)
)

# 查看一下数据
glimpse(patient_data)

# -----------------------------------------------------------------------------
# 步骤 2: 检查匹配前的基线特征不平衡性
# -----------------------------------------------------------------------------
cat("\n\n===== 匹配前基线特征表 (Table 1) =====\n")
# 定义我们关心的变量
vars_to_check <- c("age", "sex", "sbp", "has_hf", "has_diabetes")
# 创建基线表
table_before <- CreateTableOne(vars = vars_to_check, strata = "lv_size", data = patient_data, test = TRUE)
print(table_before, smd = TRUE) # smd = TRUE 会显示"标准化平均差"，是衡量不平衡性的关键指标

cat("\n观察: 匹配前，两组在年龄、血压、心衰和糖尿病比例上存在巨大差异 (SMD > 0.1)。\n")
cat("P值也都非常小，说明差异是统计学显著的。直接比较这两组会得到有偏倚的结论。\n")


# -----------------------------------------------------------------------------
# 步骤 3: 执行1对2的倾向性评分匹配 (PSM)
# -----------------------------------------------------------------------------
cat("\n\n===== 正在执行1对2倾向性评分匹配... =====\n")

# 使用MatchIt包进行匹配
# 公式 lv_size ~ ... 的意思是：基于这些协变量，预测一个人属于'Small LV'组的概率
match_obj <- matchit(
  formula = lv_size ~ age + sex + sbp + has_hf + has_diabetes,
  data = patient_data,
  method = "nearest",  # 使用最近邻匹配法
  ratio = 2,           # 这就是关键的 1:2 匹配！
  caliper = 0.1,       # 设置卡钳值，确保匹配质量
  replace = FALSE      # 无放回匹配
)

# 显示匹配结果的摘要
summary(match_obj)


# -----------------------------------------------------------------------------
# 步骤 4: 检查匹配后的基线特征平衡性
# -----------------------------------------------------------------------------
cat("\n\n===== 匹配后基线特征表 (Table 1) =====\n")

# 从匹配结果中提取出匹配后的数据
matched_data <- match.data(match_obj)

# 在匹配后的数据上再次创建基线表
table_after <- CreateTableOne(vars = vars_to_check, strata = "lv_size", data = matched_data, test = TRUE)
print(table_after, smd = TRUE)

cat("\n观察: 匹配后，所有变量的SMD都远小于0.1，P值也变大了。\n")
cat("这说明匹配非常成功！现在两组患者的基线特征已经非常相似，具备了可比性。\n")


# -----------------------------------------------------------------------------
# 步骤 5 & 6: 运行生存模型并比较结果
# -----------------------------------------------------------------------------

### 模型A: 错误的方法 - 在原始不平衡数据上使用标准Cox模型
cat("\n\n===== 模型A: 在原始数据上运行标准Cox模型 (有偏倚的结果) =====\n")
model_naive <- coxph(Surv(observed_time, event_status) ~ lv_size, data = patient_data)
summary(model_naive)
cat(paste0("风险比 (Hazard Ratio) for Small LV: ", round(exp(coef(model_naive)), 2), "\n"))


### 模型B: 正确的方法 - 在匹配后的数据上使用混合效应Cox模型
cat("\n\n===== 模型B: 在匹配数据上运行混合效应Cox模型 (校正后的结果) =====\n")

# 这里的'subclass'是matchit()自动生成的变量，它代表了每一个匹配小组（1个Small LV + 2个Normal LV）
# (1 | subclass) 就是随机效应部分，它告诉模型要考虑来自同一个匹配小组内的相关性！
model_correct <- coxme(Surv(observed_time, event_status) ~ lv_size + (1 | subclass), data = matched_data)
summary(model_correct)
# 提取固定效应的系数
fixed_effects <- fixef(model_correct)
cat(paste0("风险比 (Hazard Ratio) for Small LV: ", round(exp(fixed_effects['lv_sizeSmall LV']), 2), "\n"))

# -----------------------------------------------------------------------------
# 最终结论
# -----------------------------------------------------------------------------
cat("\n\n===== 结论比较 =====\n")
cat("模型A (错误方法) 得到的HR: ", round(exp(coef(model_naive)), 2), "\n")
cat("模型B (正确方法) 得到的HR: ", round(exp(fixed_effects['lv_sizeSmall LV']), 2), "\n")
cat("\n你会发现，模型A高估了'Small LV'的风险，因为它错误地把年龄大、疾病多等因素的风险都归咎于'Small LV'了。\n")
cat("模型B通过匹配和混合效应模型，排除了这些混杂因素的干扰，得到了一个更真实、更可信的风险评估。\n")
cat("这完美地展示了这篇论文方法学的核心思想和价值！\n")