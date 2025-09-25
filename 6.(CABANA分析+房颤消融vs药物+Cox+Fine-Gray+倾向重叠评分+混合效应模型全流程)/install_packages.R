# 安装CABANA分析所需的R包

cat("开始安装R包...\n")

# 定义需要的包
required_packages <- c(
  "tidyverse", "survival", "survminer", "cmprsk", 
  "WeightIt", "lme4", "lmerTest", "tableone", 
  "forestplot", "RColorBrewer", "gridExtra", 
  "knitr", "kableExtra", "patchwork", "scales"
)

# 检查并安装缺失的包
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("安装包:", pkg, "\n")
      tryCatch({
        install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org/")
        library(pkg, character.only = TRUE)
        cat("✓", pkg, "安装成功\n")
      }, error = function(e) {
        cat("✗", pkg, "安装失败:", e$message, "\n")
      })
    } else {
      cat("✓", pkg, "已安装\n")
    }
  }
}

install_if_missing(required_packages)

cat("\n包安装检查完成！\n")
