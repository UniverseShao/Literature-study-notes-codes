library(dagitty)
library(ggdag)
library(tidyverse)
library(survival)
library(patchwork)  # 用于图形拼接

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

# ===========================================================================
# 可视化部分：创建多个专业的DAG图形
# ===========================================================================

# 图1: 基础DAG结构图（带节点标签和颜色）
cat("\n正在生成图1: 基础DAG结构图...\n")

# 创建节点标签映射（中英文对照）
node_labels <- c(
  "SmallLV" = "Small Left Ventricle\n(小左心室)",
  "CVEvent" = "Cardiovascular Event\n(心血管事件)",
  "Age" = "Age\n(年龄)",
  "Sex" = "Sex\n(性别)",
  "Hypertension" = "Hypertension\n(高血压)",
  "Diabetes" = "Diabetes\n(糖尿病)",
  "CKD" = "CKD\n(慢性肾病)",
  "Medication" = "Medication\n(药物治疗)"
)

# 将DAG转换为tidy格式
dag_tidy <- tidy_dagitty(my_dag)

# 手动分类节点类型
dag_tidy$data <- dag_tidy$data %>%
  mutate(
    node_type = case_when(
      name == "SmallLV" ~ "Exposure",
      name == "CVEvent" ~ "Outcome",
      name %in% c("Age", "Sex") ~ "Confounder",
      name == "CKD" ~ "Mediator",
      name == "Medication" ~ "Collider",
      TRUE ~ "Other"
    ),
    label = node_labels[name]
  )

# 设置节点颜色
node_colors <- c(
  "Exposure" = "#E41A1C",      # 红色 - 暴露
  "Outcome" = "#377EB8",       # 蓝色 - 结局
  "Confounder" = "#4DAF4A",    # 绿色 - 混杂
  "Mediator" = "#FF7F00",      # 橙色 - 中介
  "Collider" = "#984EA3",      # 紫色 - 对撞
  "Other" = "#999999"          # 灰色 - 其他
)

p1 <- ggplot(dag_tidy, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(edge_color = "grey60", edge_width = 1, 
                 arrow_directed = grid::arrow(length = grid::unit(10, "pt"), type = "closed")) +
  geom_dag_point(aes(color = node_type), size = 20, alpha = 0.9) +
  geom_dag_text(aes(label = label), color = "white", size = 3.2, fontface = "bold", lineheight = 0.85) +
  scale_color_manual(
    values = node_colors,
    name = "Variable Type",
    labels = c("Confounder (混杂因子)", "Exposure (暴露)", 
               "Mediator (中介因子)", "Outcome (结局)", "Collider (对撞因子)")
  ) +
  labs(title = "Figure 1. Directed Acyclic Graph (DAG) for Causal Relationships",
       subtitle = "Exploring the relationship between small left ventricle and cardiovascular events") +
  theme_dag() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

print(p1)

# 保存图1
ggsave("DAG_Figure1_Structure.png", p1, width = 12, height = 8, dpi = 300, bg = "white")
ggsave("DAG_Figure1_Structure.pdf", p1, width = 12, height = 8, bg = "white")

# 图2: 调整集可视化（高亮显示需要调整的变量）
cat("\n正在生成图2: 调整集可视化图...\n")

# 核心步骤：找到关闭所有后门路径的调整集
adjustment_sets <- adjustmentSets(my_dag, type = "minimal")
cat("\n--- 最小充分调整集 ---\n")
print(adjustment_sets)

# 创建调整集数据
adj_dag <- dag_tidy
adj_dag$data <- adj_dag$data %>%
  mutate(
    adjusted = name %in% c("Age", "Sex"),
    node_status = case_when(
      name == "SmallLV" ~ "Exposure",
      name == "CVEvent" ~ "Outcome",
      name %in% c("Age", "Sex") ~ "Adjust for",
      name == "CKD" ~ "Do NOT adjust\n(Mediator)",
      name == "Medication" ~ "Do NOT adjust\n(Collider)",
      TRUE ~ "Other"
    )
  )

status_colors <- c(
  "Exposure" = "#E41A1C",
  "Outcome" = "#377EB8",
  "Adjust for" = "#4DAF4A",
  "Do NOT adjust\n(Mediator)" = "#FF7F00",
  "Do NOT adjust\n(Collider)" = "#984EA3",
  "Other" = "#FFED6F"
)

p2 <- ggplot(adj_dag, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(edge_color = "grey60", edge_width = 1,
                 arrow_directed = grid::arrow(length = grid::unit(10, "pt"), type = "closed")) +
  geom_dag_point(aes(color = node_status), size = 20, alpha = 0.9) +
  geom_dag_text(aes(label = label), color = "white", size = 3.2, fontface = "bold", lineheight = 0.85) +
  scale_color_manual(
    values = status_colors,
    name = "Adjustment Strategy"
  ) +
  labs(title = "Figure 2. Minimal Sufficient Adjustment Set",
       subtitle = "Variables to adjust: Age and Sex | Do NOT adjust: CKD (mediator) and Medication (collider)") +
  theme_dag() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

print(p2)

# 保存图2
ggsave("DAG_Figure2_AdjustmentSet.png", p2, width = 12, height = 8, dpi = 300, bg = "white")
ggsave("DAG_Figure2_AdjustmentSet.pdf", p2, width = 12, height = 8, bg = "white")

# 图3: 路径分析图（显示关键路径）
cat("\n正在生成图3: 因果路径分析图...\n")

# 创建路径标注
path_dag <- dag_tidy
path_dag$data <- path_dag$data %>%
  mutate(
    path_role = case_when(
      name %in% c("SmallLV", "CVEvent") ~ "Key variables",
      name == "CKD" ~ "Mediator",
      name %in% c("Age", "Sex") ~ "Confounder (adjust)",
      name %in% c("Hypertension", "Diabetes") ~ "Other confounders",
      name == "Medication" ~ "Collider (do not adjust)",
      TRUE ~ "Other"
    )
  )

path_colors <- c(
  "Key variables" = "#E41A1C",
  "Mediator" = "#FF7F00",
  "Confounder (adjust)" = "#4DAF4A",
  "Other confounders" = "#FFED6F",
  "Collider (do not adjust)" = "#984EA3",
  "Other" = "grey70"
)

p3 <- ggplot(path_dag, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(edge_color = "grey60", edge_width = 1.2,
                 arrow_directed = grid::arrow(length = grid::unit(10, "pt"), type = "closed")) +
  geom_dag_point(aes(color = path_role), size = 22, alpha = 0.9) +
  geom_dag_text(aes(label = label), color = "white", size = 3.2, fontface = "bold", lineheight = 0.85) +
  scale_color_manual(
    values = path_colors,
    name = "Variable Role in Pathway"
  ) +
  labs(title = "Figure 3. Complete Causal Pathway Analysis",
       subtitle = "Identifying all causal paths: direct, indirect (mediation), and backdoor (confounding)") +
  theme_dag() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

print(p3)

# 保存图3
ggsave("DAG_Figure3_PathwaysAnalysis.png", p3, width = 12, height = 8, dpi = 300, bg = "white")
ggsave("DAG_Figure3_PathwaysAnalysis.pdf", p3, width = 12, height = 8, bg = "white")

# 图4: 拼接所有图形
cat("\n正在生成图4: 综合拼接图...\n")

p_combined <- p1 / p2 / p3 +
  plot_annotation(
    title = "Comprehensive DAG Analysis for Small Left Ventricle and Cardiovascular Events",
    subtitle = "A complete causal inference framework using Directed Acyclic Graphs",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey40")
    )
  )

print(p_combined)

# 保存拼接图
ggsave("DAG_Combined_AllFigures.png", p_combined, width = 12, height = 20, dpi = 300, bg = "white")
ggsave("DAG_Combined_AllFigures.pdf", p_combined, width = 12, height = 20, bg = "white")

cat("\n========================================\n")
cat("所有DAG可视化图形已保存\n")
cat("========================================\n")
cat("- DAG_Figure1_Structure.png/pdf (基础结构图)\n")
cat("- DAG_Figure2_AdjustmentSet.png/pdf (调整集图)\n")
cat("- DAG_Figure3_PathwaysAnalysis.png/pdf (路径分析图)\n")
cat("- DAG_Combined_AllFigures.png/pdf (综合拼接图)\n")
cat("========================================\n\n")

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