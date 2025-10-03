library(dagitty)
library(ggdag)
library(tidyverse)
library(survival)

# 使用 dagitty() 函数创建一个DAG对象
# 语法: 'Var1 -> Var2' 表示 Var1 是 Var2 的原因
# {Var1 Var2} -> Var3 表示 Var1 和 Var2 都是 Var3 的原因
# @exposure 指定暴露变量
# @outcome 指定结局变量

# 核心因果路径: SmallLV影响CVEvent
# 混杂因子: 年龄和性别影响SmallLV和CVEvent
# 中介因子: CKD介导SmallLV对CVEvent的影响
# 对撞因子: Medication受Hypertension和Diabetes影响
my_dag <- dagitty("dag {
  SmallLV -> CVEvent
  Age -> SmallLV
  Age -> CVEvent
  Sex -> SmallLV
  Sex -> CVEvent
  Age -> Hypertension
  Hypertension -> Diabetes
  Diabetes -> CVEvent
  CKD -> CVEvent
  SmallLV -> CKD
  Hypertension -> Medication
  Diabetes -> Medication
  Medication -> CVEvent
}")

# 设置暴露和结局变量
exposures(my_dag) <- "SmallLV"
outcomes(my_dag) <- "CVEvent"

# 打印DAG对象，查看基本信息
print(my_dag)

# 可视化DAG
ggdag(my_dag, text = TRUE) +
  theme_dag()

# 核心步骤：找到关闭所有后门路径的调整集
adjustment_sets <- adjustmentSets(my_dag, type = "minimal")
cat("--- 最小充分调整集 ---\n")
print(adjustment_sets)

# 可视化这个调整集
# 图中会高亮显示需要调整的变量
ggdag_adjustment_set(my_dag, shadow = TRUE) +
  theme_dag_blank()

# --- 模拟数据 (仅为演示) ---
set.seed(123)
n <- 2000
sim_data <- tibble(
  Age = rnorm(n, 65, 8),
  Sex = factor(rbinom(n, 1, 0.6), labels = c("Female", "Male")),
  Hypertension = factor(rbinom(n, 1, 0.5 + 0.01 * (Age - 65))),
  Diabetes = factor(rbinom(n, 1, 0.2 + 0.3 * (Hypertension == "1"))),
  Small_LV = factor(rbinom(n, 1, 0.1 + 0.2 * (Sex == "Female") - 0.01 * (Age - 65))),
  CKD = factor(rbinom(n, 1, 0.1 + 0.3 * (Small_LV == "1"))),
  Medication = factor(rbinom(n, 1, 0.5 + 0.2*(Hypertension=="1") - 0.2*(Diabetes=="1"))),
  
  # 结局模拟
  time = rexp(n, rate = 0.1 * exp(0.5*(Small_LV=="1") + 0.03*(Age-65) - 0.1*(Sex=="Male") + 0.2*(Diabetes=="1") + 0.3*(CKD=="1"))),
  status = rbinom(n, 1, 0.8)
)

# --- 运行模型 ---

# 模型1: DAG指导的、更科学的模型
# 只纳入 Small_LV 和 最小充分调整集中的变量
cat("\n\n--- 模型1: DAG指导的Cox回归模型 ---\n")
model_dag <- coxph(
  Surv(time, status) ~ Small_LV + Age + Sex + Hypertension + Diabetes,
  data = sim_data
)
summary(model_dag)


# 模型2: 传统“厨房水槽式”模型 (Kitchen Sink)
# 错误地纳入了所有看似相关的变量，包括中介和对撞因子
cat("\n\n--- 模型2: 错误的'全家桶'Cox回归模型 ---\n")
model_kitchen_sink <- coxph(
  Surv(time, status) ~ Small_LV + Age + Sex + Hypertension + Diabetes + CKD + Medication,
  data = sim_data
)
summary(model_kitchen_sink)

# --- 结果对比 ---
cat("\n\n--- 结果对比 (Small_LV的风险比) ---\n")
cat("DAG指导模型的HR:", exp(coef(model_dag)["Small_LV1"]), "\n")
cat("错误模型的HR:", exp(coef(model_kitchen_sink)["Small_LV1"]), "\n")