# 安装必要的包 (如果尚未安装)
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("mice")

# 加载库
library(dplyr)
library(tidyr)
library(mice)

# 创建一个模拟数据集
# NA 代表缺失值
df <- data.frame(
  age = c(45, 62, 38, 55, 71, 49, 58, 67),
  sex = as.factor(c('male', 'female', 'female', 'male', 'female', 'male', 'male', 'female')),
  educational_attainment = as.factor(c(
    'college or above', 
    'high school or equivalent', 
    NA, # 缺失值
    'college or above', 
    'less than high school', 
    'high school or equivalent',
    NA, # 缺失值
    'college or above'
  )),
  smoking_status = as.factor(c(
    'former', 
    'never', 
    'current', 
    NA, # 缺失值
    'former', 
    'current', 
    'never', 
    NA # 缺失值
  )),
  hypertension = c(1, 0, 1, 1, 0, 0, 1, 1)
)

cat("--- 原始数据集 (Original DataFrame) ---\n")
print(df)
cat("\n--- 缺失值统计 ---\n")
print(sapply(df, function(x) sum(is.na(x))))

# 方法一：创建“未知 (Unknown)”类别 (主要分析)-------------------------------------------------

# 创建一个数据副本以进行操作
df_primary <- df

# 找出包含缺失值的分类变量列
categorical_cols_with_na <- c('educational_attainment', 'smoking_status')

# 使用 tidyr::replace_na() 填充缺失值
# 对于因子，我们首先需要添加 "Unknown" 作为一个新的因子水平
df_primary <- df_primary %>%
  mutate(across(all_of(categorical_cols_with_na), function(col) {
    # 添加新的因子水平
    col_with_new_level <- factor(col, levels = c(levels(col), "Unknown"))
    # 替换 NA
    replace_na(col_with_new_level, "Unknown")
  }))

cat("\n--- 方法一：创建 'Unknown' 类别 ---\n")
cat("处理后的数据集 (DataFrame after Method 1):\n")
print(df_primary)
cat("\n处理后的因子水平:\n")
print(levels(df_primary$smoking_status))
cat("\n处理后已无缺失值:\n")
print(sapply(df_primary, function(x) sum(is.na(x))))


# 方法二：多重插补法 (MI) (敏感性分析)-----------------------------------------------------
cat("\n--- 方法二：多重插补法 (MI) 使用 mice 包 ---\n")

# 执行多重插补
# m = 4 表示生成4个插补后的数据集，与论文一致
# seed = 123 确保结果可以复现
set.seed(123) # for reproducibility
imputed_mice <- mice(df, m = 4, method = 'polyreg', printFlag = FALSE)

# 'polyreg' (多项式回归) 是处理因子类型缺失值的常用方法

# 查看插补对象的摘要
cat("\n--- MICE 插补对象摘要 ---\n")
print(imputed_mice)

# 从插补对象中提取一个完整的数据集 (例如，第一个)
df_imputed_1 <- complete(imputed_mice, 1)

cat("\n--- 提取出的第一个完整数据集 (1st Imputed DataFrame) ---\n")
print(df_imputed_1)

# 你也可以提取所有4个数据集到一个列表中
# all_imputed_datasets <- complete(imputed_mice, "long", include = FALSE)



