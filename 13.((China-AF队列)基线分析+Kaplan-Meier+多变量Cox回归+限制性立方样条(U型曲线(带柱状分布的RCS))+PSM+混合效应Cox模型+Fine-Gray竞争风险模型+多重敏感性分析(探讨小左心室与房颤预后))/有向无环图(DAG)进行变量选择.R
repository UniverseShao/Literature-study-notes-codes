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

# 基于文献的完整DAG模型
# 包含左心室大小(LV_size)与房颤患者主要不良事件(Major_adverse_event)的因果关系
# 以及所有相关的混杂因子、中介因子和治疗变量
my_dag <- dagitty("dag {
  LV_size -> Major_adverse_event
  LV_size -> LVEF
  LV_size -> HF
  LV_size -> Ischemic_stroke
  
  Age -> LV_size
  Age -> Major_adverse_event
  Age -> CAD
  Age -> CKD
  Age -> HF
  Age -> AF_type
  
  Sex -> LV_size
  Sex -> Major_adverse_event
  Sex -> CAD
  Sex -> HF
  
  BMI -> LV_size
  BMI -> DM
  BMI -> Hypertension
  BMI -> HF
  
  Hypertension -> LV_size
  Hypertension -> CAD
  Hypertension -> CKD
  Hypertension -> HF
  Hypertension -> Ischemic_stroke
  Hypertension -> Major_adverse_event
  Hypertension -> ACEI_ARB
  
  DM -> CAD
  DM -> CKD
  DM -> PAD
  DM -> Major_adverse_event
  
  CAD -> HF
  CAD -> Major_adverse_event
  CAD -> Ischemic_stroke
  
  CKD -> HF
  CKD -> Major_adverse_event
  CKD -> LVEF
  
  PAD -> Major_adverse_event
  
  HF -> Major_adverse_event
  HF -> LVEF
  
  AF_type -> Major_adverse_event
  AF_type -> Ischemic_stroke
  AF_type -> OAC
  
  LVEF -> Major_adverse_event
  LVEF -> HF
  
  LAD -> Major_adverse_event
  LAD -> AF_type
  LAD -> Ischemic_stroke
  
  LVEDD -> LV_size
  LVEDD -> LVEF
  LVEDD -> HF
  
  IVS -> LV_size
  IVS -> Hypertension
  
  LVPW -> LV_size
  LVPW -> Hypertension
  
  Ischemic_stroke -> Major_adverse_event
  
  OAC -> Ischemic_stroke
  OAC -> Major_adverse_event
  
  AAD -> Major_adverse_event
  AAD -> AF_type
  
  ACEI_ARB -> HF
  ACEI_ARB -> Major_adverse_event
  ACEI_ARB -> LVEF
}")

# 设置暴露和结局变量
exposures(my_dag) <- "LV_size"
outcomes(my_dag) <- "Major_adverse_event"

# 打印DAG对象，查看基本信息
print(my_dag)

# ===========================================================================
# 可视化部分：创建多个专业的DAG图形
# ===========================================================================

# 图1: 基础DAG结构图（带节点标签和颜色）
cat("\n正在生成图1: 基础DAG结构图...\n")

# 创建节点标签映射（中英文对照）
node_labels <- c(
  "LV_size" = "LV size",
  "Major_adverse_event" = "Major\nadverse\nevent",
  "Age" = "Age",
  "Sex" = "Sex",
  "BMI" = "BMI",
  "Hypertension" = "Hyper-\ntension",
  "DM" = "DM",
  "CAD" = "CAD",
  "CKD" = "CKD",
  "PAD" = "PAD",
  "HF" = "HF",
  "AF_type" = "AF type",
  "LVEF" = "LVEF",
  "LAD" = "LAD",
  "LVEDD" = "LVEDD",
  "IVS" = "IVS",
  "LVPW" = "LVPW",
  "Ischemic_stroke" = "Ischemic\nstroke",
  "OAC" = "OAC",
  "AAD" = "AAD",
  "ACEI_ARB" = "ACEI/ARB"
)

# 将DAG转换为tidy格式
dag_tidy <- tidy_dagitty(my_dag)

# 手动分类节点类型
dag_tidy$data <- dag_tidy$data %>%
  mutate(
    node_type = case_when(
      name == "LV_size" ~ "Exposure",
      name == "Major_adverse_event" ~ "Outcome",
      name %in% c("Age", "Sex", "BMI") ~ "Baseline",
      name %in% c("Hypertension", "DM", "CAD", "CKD", "PAD") ~ "Comorbidity",
      name %in% c("HF", "LVEF", "Ischemic_stroke") ~ "Mediator",
      name %in% c("LAD", "LVEDD", "IVS", "LVPW", "AF_type") ~ "Cardiac",
      name %in% c("OAC", "AAD", "ACEI_ARB") ~ "Treatment",
      TRUE ~ "Other"
    ),
    label = node_labels[name]
  )

# 设置节点颜色
node_colors <- c(
  "Exposure" = "#90EE90",      # 浅绿色 - 暴露
  "Outcome" = "#87CEEB",       # 天蓝色 - 结局
  "Baseline" = "#FFB6C1",      # 浅粉色 - 基线变量
  "Comorbidity" = "#DDA0DD",   # 淡紫色 - 合并症
  "Mediator" = "#87CEEB",      # 蓝色 - 中介
  "Cardiac" = "#FFB6C1",       # 粉色 - 心脏参数
  "Treatment" = "#87CEEB",     # 蓝色 - 治疗
  "Other" = "#999999"          # 灰色 - 其他
)

p1 <- ggplot(dag_tidy, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(edge_color = "#B8336A", edge_width = 0.8, 
                 arrow_directed = grid::arrow(length = grid::unit(8, "pt"), type = "closed")) +
  geom_dag_point(aes(color = node_type), size = 12, alpha = 0.95) +
  geom_dag_text(aes(label = label), color = "black", size = 2.8, fontface = "bold", lineheight = 0.9) +
  scale_color_manual(
    values = node_colors,
    name = "Variable Type",
    labels = c("Baseline", "Cardiac parameters", "Comorbidities", 
               "Exposure (LV size)", "Mediators", "Outcome", "Treatments")
  ) +
  labs(title = "Supplemental Figure 1. Directed Acyclic Graph",
       subtitle = "Causal relationships between left ventricular size and major adverse events in atrial fibrillation") +
  theme_dag() +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 10, hjust = 0, color = "grey30"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(15, 15, 15, 15)
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
# 先获取调整集中的变量名称
adjustment_vars <- unlist(lapply(adjustment_sets, function(x) as.character(x)))
adj_dag$data <- adj_dag$data %>%
  mutate(
    adjusted = name %in% adjustment_vars,
    node_status = case_when(
      name == "LV_size" ~ "Exposure",
      name == "Major_adverse_event" ~ "Outcome",
      name %in% adjustment_vars ~ "Adjust for",
      name %in% c("HF", "LVEF", "Ischemic_stroke") ~ "Mediator\n(Do NOT adjust)",
      TRUE ~ "Other"
    )
  )

status_colors <- c(
  "Exposure" = "#90EE90",
  "Outcome" = "#87CEEB",
  "Adjust for" = "#FFD700",
  "Mediator\n(Do NOT adjust)" = "#FFA500",
  "Other" = "#DDA0DD"
)

p2 <- ggplot(adj_dag, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(edge_color = "#B8336A", edge_width = 0.8,
                 arrow_directed = grid::arrow(length = grid::unit(8, "pt"), type = "closed")) +
  geom_dag_point(aes(color = node_status), size = 12, alpha = 0.95) +
  geom_dag_text(aes(label = label), color = "black", size = 2.8, fontface = "bold", lineheight = 0.9) +
  scale_color_manual(
    values = status_colors,
    name = "Adjustment Strategy"
  ) +
  labs(title = "Figure 2. Minimal Sufficient Adjustment Set for LV Size Effect",
       subtitle = paste("Variables to adjust:", paste(adjustment_vars, collapse = ", "))) +
  theme_dag() +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 10, hjust = 0, color = "grey30"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(15, 15, 15, 15)
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
      name %in% c("LV_size", "Major_adverse_event") ~ "Exposure & Outcome",
      name %in% c("HF", "LVEF", "Ischemic_stroke") ~ "Mediator (causal pathway)",
      name %in% adjustment_vars ~ "Confounder (adjust)",
      name %in% c("Hypertension", "DM", "CAD", "CKD", "PAD") ~ "Comorbidities",
      name %in% c("OAC", "AAD", "ACEI_ARB") ~ "Treatments",
      name %in% c("LAD", "LVEDD", "IVS", "LVPW", "AF_type") ~ "Cardiac parameters",
      TRUE ~ "Other"
    )
  )

path_colors <- c(
  "Exposure & Outcome" = "#E74C3C",
  "Mediator (causal pathway)" = "#FFA500",
  "Confounder (adjust)" = "#4DAF4A",
  "Comorbidities" = "#9B59B6",
  "Treatments" = "#3498DB",
  "Cardiac parameters" = "#F39C12",
  "Other" = "grey70"
)

p3 <- ggplot(path_dag, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(edge_color = "#B8336A", edge_width = 0.8,
                 arrow_directed = grid::arrow(length = grid::unit(8, "pt"), type = "closed")) +
  geom_dag_point(aes(color = path_role), size = 12, alpha = 0.95) +
  geom_dag_text(aes(label = label), color = "black", size = 2.8, fontface = "bold", lineheight = 0.9) +
  scale_color_manual(
    values = path_colors,
    name = "Variable Role in Pathway"
  ) +
  labs(title = "Figure 3. Complete Causal Pathway Analysis",
       subtitle = "Direct effects, mediation pathways, and confounding relationships") +
  theme_dag() +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 10, hjust = 0, color = "grey30"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(15, 15, 15, 15)
  )

print(p3)

# 保存图3
ggsave("DAG_Figure3_PathwaysAnalysis.png", p3, width = 12, height = 8, dpi = 300, bg = "white")
ggsave("DAG_Figure3_PathwaysAnalysis.pdf", p3, width = 12, height = 8, bg = "white")

# 图4: 拼接所有图形
cat("\n正在生成图4: 综合拼接图...\n")

p_combined <- p1 / p2 / p3 +
  plot_annotation(
    title = "Comprehensive DAG Analysis: Left Ventricular Size and Major Adverse Events in Atrial Fibrillation",
    subtitle = "A complete causal inference framework for variable selection and confounding control",
    theme = theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
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

# ===========================================================================
# DAG分析结果总结
# ===========================================================================
cat("\n\n=== DAG Analysis Summary ===\n\n")

cat("1. EXPOSURE VARIABLE:\n")
cat("   - Left ventricular size (LV_size)\n\n")

cat("2. OUTCOME VARIABLE:\n")
cat("   - Major adverse event (Major_adverse_event)\n\n")

cat("3. MINIMAL SUFFICIENT ADJUSTMENT SET:\n")
cat("   Variables to adjust for:", paste(adjustment_vars, collapse = ", "), "\n")
cat("   These variables will block all backdoor paths and control confounding.\n\n")

cat("4. MEDIATORS (Do NOT adjust if studying total effect):\n")
cat("   - Heart failure (HF)\n")
cat("   - Left ventricular ejection fraction (LVEF)\n")
cat("   - Ischemic stroke\n")
cat("   Note: Adjusting for mediators will block causal pathways and bias the estimate.\n\n")

cat("5. CAUSAL PATHWAYS:\n")
cat("   Direct path: LV_size -> Major_adverse_event\n")
cat("   Indirect paths through mediators:\n")
cat("     - LV_size -> HF -> Major_adverse_event\n")
cat("     - LV_size -> LVEF -> Major_adverse_event\n")
cat("     - LV_size -> Ischemic_stroke -> Major_adverse_event\n\n")

cat("6. MODELING RECOMMENDATIONS:\n")
cat("   For Cox regression model:\n")
cat("   Model: Surv(time, event) ~ LV_size +", paste(adjustment_vars, collapse=" + "), "\n")
cat("   This model provides an unbiased estimate of the total causal effect.\n\n")

cat("========================================\n")
cat("Analysis complete!\n")
cat("========================================\n")