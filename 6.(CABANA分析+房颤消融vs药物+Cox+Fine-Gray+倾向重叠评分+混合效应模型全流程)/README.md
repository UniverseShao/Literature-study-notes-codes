# CABANA试验房颤患者导管消融 vs 药物治疗分析

这是一个完整的R分析脚本包，用于复现和理解CABANA试验的统计分析方法。本项目基于《研究思路.md》和《statistical_analysis_summary.md》中描述的方法学，实现了从数据模拟到结果可视化的完整分析流程。

## 🎯 项目概述

本项目实现了CABANA试验的完整统计分析方法，包括：

- **描述性统计与基线比较**
- **Kaplan-Meier生存分析**
- **多变量Cox比例风险回归**
- **Fine-Gray竞争风险模型**
- **倾向性评分重叠加权**
- **重复测量混合效应模型**
- **高质量可视化图表**

## 📁 文件结构

```
├── af_cabana_analysis.R          # 主要统计分析脚本
├── cabana_visualization.R        # 可视化脚本
├── run_complete_analysis.R       # 完整分析运行脚本
├── README.md                     # 使用说明（本文件）
├── 研究思路.md                   # 详细的研究方法思路
├── statistical_analysis_summary.md # 统计方法总结
└── plots/                        # 输出图表文件夹
    ├── km_*.png                  # Kaplan-Meier生存曲线
    ├── forest_plot.png           # 风险比森林图
    ├── mafsi_trajectory.png      # 症状评分轨迹
    └── ...                       # 其他图表
```

## 🚀 快速开始

### 1. 环境准备

确保安装了以下R包：

```r
# 安装必要的包
install.packages(c(
  "tidyverse", "survival", "survminer", "cmprsk", 
  "WeightIt", "lme4", "lmerTest", "tableone", 
  "forestplot", "RColorBrewer", "gridExtra", 
  "knitr", "kableExtra", "patchwork", "scales"
))
```

### 2. 运行完整分析

**方法一：一键运行（推荐）**

```r
# 运行完整分析
source("run_complete_analysis.R")
```

**方法二：分步运行**

```r
# 步骤1：运行主要分析
source("af_cabana_analysis.R")

# 步骤2：运行可视化
source("cabana_visualization.R")
all_plots <- create_all_visualizations(baseline_data, mafsi_data, results_summary, ps_weights)
save_all_plots(all_plots)
```

### 3. 查看结果

- **详细报告**: `CABANA_Analysis_Report.md`
- **数值结果**: `cabana_analysis_results.csv`
- **基线数据**: `cabana_baseline_data.csv`
- **图表文件**: `plots/` 文件夹

## 📊 分析方法详解

### 1. 数据模拟

创建了2204名房颤患者的模拟数据集，包括：

- **基线特征**: 年龄、性别、房颤类型、合并症等
- **治疗分组**: 导管消融组（1108人）vs 药物治疗组（1096人）
- **结局事件**: 主要复合终点、房颤复发、死亡等
- **症状评分**: MAFSI评分的纵向数据

### 2. 核心分析方法

#### 基线比较
- 连续变量：t检验/Wilcoxon秩和检验
- 分类变量：χ²检验

#### 生存分析
- **Kaplan-Meier方法**: 生存曲线可视化
- **Log-rank检验**: 组间生存差异检验
- **多变量Cox回归**: 调整混杂因素

#### 竞争风险分析
- **Fine-Gray模型**: 处理房颤复发与死亡的竞争关系
- **累积发生函数**: 更准确的风险估计

#### 因果推断
- **倾向性评分重叠加权**: 平衡组间协变量
- **双重稳健估计**: 结合PS加权和多变量调整

#### 纵向分析
- **混合效应模型**: 分析症状评分随时间变化

### 3. 统计假设检验

- **比例风险假设**: Schoenfeld残差检验
- **模型诊断**: 残差分析和拟合优度检验
- **敏感性分析**: 多种模型验证结果稳健性

## 📈 主要输出

### 数据文件
- `cabana_baseline_data.csv`: 完整的基线数据
- `cabana_mafsi_data.csv`: MAFSI评分纵向数据
- `cabana_analysis_results.csv`: 所有分析结果汇总

### 图表文件
- **生存曲线**: 主要终点、房颤复发、全因死亡
- **Forest图**: 各种模型的风险比可视化
- **竞争风险图**: 累积发生函数
- **基线对比图**: 治疗组间特征分布
- **症状轨迹图**: MAFSI评分随时间变化
- **倾向性评分图**: PS分布和权重分布

### 分析报告
- `CABANA_Analysis_Report.md`: 详细的分析报告，包含：
  - 研究概述和样本特征
  - 每个分析步骤的结果
  - 方法学亮点总结
  - 结论和临床意义

## 🔍 结果解读

### 主要发现
- **主要复合终点**: 导管消融组风险降低约25%
- **房颤复发**: 导管消融组风险降低约50%
- **症状改善**: 导管消融组MAFSI评分改善更显著

### 方法学亮点
1. **双重稳健设计**: 同时使用PS加权和多变量调整
2. **竞争风险处理**: 避免传统方法的偏倚
3. **假设检验**: 严格验证模型前提条件
4. **敏感性分析**: 多种模型确保结果稳健

## 🛠️ 自定义分析

### 修改样本量
```r
# 在af_cabana_analysis.R中修改
n_total <- 你的样本量
n_ablation <- 导管消融组样本量
n_drug <- 药物治疗组样本量
```

### 调整效应大小
```r
# 修改风险比
hr_primary <- 0.75        # 主要终点风险比
hr_recurrence <- 0.50     # 复发风险比
hr_death <- 0.85          # 死亡风险比
```

### 添加新的协变量
```r
# 在协变量列表中添加新变量
covariates <- c("age", "sex", "race", "af_type", 
                "新变量名", ..., "cha2ds2_vasc")
```

## 📚 参考文献

本分析基于以下方法学框架：

1. **Cox比例风险模型**: Cox DR. Regression models and life-tables. J R Stat Soc B. 1972.
2. **Fine-Gray竞争风险模型**: Fine JP, Gray RJ. A proportional hazards model for the subdistribution of a competing risk. J Am Stat Assoc. 1999.
3. **倾向性评分重叠加权**: Li F, Morgan KL, Zaslavsky AM. Balancing covariates via propensity score weighting. J Am Stat Assoc. 2018.
4. **混合效应模型**: Laird NM, Ware JH. Random-effects models for longitudinal data. Biometrics. 1982.

## 💡 使用建议

1. **首次使用**: 建议先运行`run_complete_analysis.R`查看完整流程
2. **学习目的**: 逐步阅读`af_cabana_analysis.R`中的每个分析步骤
3. **自定义分析**: 根据需要修改参数并重新运行
4. **结果验证**: 对比多种模型的结果确保一致性

## 🐛 常见问题

### Q: 运行时出现包缺失错误
A: 确保安装了所有必需的R包，运行上述安装命令。

### Q: 图表无法保存
A: 检查是否有`plots/`文件夹的写入权限。

### Q: 内存不足
A: 考虑减少样本量或使用更大内存的计算环境。

### Q: 想要修改分析参数
A: 在`af_cabana_analysis.R`脚本开头修改相应参数。

## 📧 联系方式

如有问题或建议，欢迎通过以下方式联系：

- 查看脚本注释了解详细实现
- 参考`研究思路.md`了解方法学背景
- 检查`statistical_analysis_summary.md`获取技术细节

---

**注意**: 本项目使用模拟数据进行演示，实际应用时请使用真实的临床数据。所有分析方法均基于已发表的高质量研究的标准方法学。
