library(survival)
library(dplyr)
set.seed(42) # 设置随机种子以保证每次运行结果都相同
n_patients <- 1500 # 模拟1500名患者
# 1.1 创建基线数据集 (wide format) ---------------------------------------------
baseline_data <- tibble(
  id = 1:n_patients,
  # 随机分配治疗组: 0 = 标准治疗, 1 = 强化治疗
  treatment = sample(0:1, n_patients, replace = TRUE),
  age = round(rnorm(n_patients, mean = 68, sd = 6)),
  sex = sample(c("Male", "Female"), n_patients, replace = TRUE, prob = c(0.6, 0.4)),
  # 基线时的LVH指数 (Peguero-Lo Presti index, 单位 μV)
  lvh_baseline = rnorm(n_patients, mean = 2600, sd = 450)
)
# 1.2 创建包含随访信息的长格式数据集 (long format) -----------------------------
# 假设在基线(time=0), 第1, 2, 3年有LVH测量值
long_data <- baseline_data %>%
  slice(rep(1:n(), each = 4)) %>%
  group_by(id) %>%
  mutate(visit_time = c(0, 1, 2, 3)) %>%
  ungroup() %>%
  mutate(
    # 模拟LVH指数的动态变化：
    # 强化治疗组(treatment=1)的LVH指数每年下降得更多
    # 这就人为地创造了“路径 a” (Treatment -> LVH)
    lvh_index = case_when(
      visit_time == 0 ~ lvh_baseline,
      # 强化组每年平均多降100, 再加上一些随机波动
      TRUE ~ lvh_baseline - visit_time * (50 + 100 * treatment) + rnorm(n(), sd = 150)
    )
  )

# 1.3 模拟生存时间和事件 (心血管事件) ------------------------------------------
# 事件风险依赖于治疗组以及当前的LVH指数
event_data <- baseline_data %>%
  mutate(
    # 模拟事件发生时间。风险依赖于治疗和基线LVH
    # 强化治疗(treatment=1)降低风险, 高LVH增加风险
    # log(hazard) = -0.6 * treatment + 0.001 * (lvh_baseline - 2600)
    true_time = rexp(n(), rate = 0.04 * exp(-0.6 * treatment + 0.001 * (lvh_baseline - 2600))),
    # 模拟一个删失时间(0.1到4年)，代表失访或研究结束
    censor_time = runif(n(), min = 0.1, max = 4),
    # 最终观察时间是真实事件时间和删失时间中的较小者
    observed_time = pmin(true_time, censor_time),
    # 状态：1 = 发生事件, 0 = 删失
    status = as.numeric(true_time <= censor_time)
  ) %>%
  dplyr::select(id, observed_time, status)

# 1.4 将数据转换为适用于时变协变量Cox模型的格式 (start, stop, event) -----------

# 这是最关键的数据准备步骤，使用survival包的tmerge函数
final_analysis_data <- tmerge(
  data1 = baseline_data,
  data2 = event_data,
  id = id,
  event = event(observed_time, status)
)

# 将时变的LVH指数作为时变协变量(tdc)合并进来
final_analysis_data <- tmerge(
  data1 = final_analysis_data,
  data2 = long_data,
  id = id,
  lvh = tdc(visit_time, lvh_index) # tdc() 表示lvh是一个时变协变量
)

#  2.1 拟合模型 1: 估计总效应 (路径 c: Treatment -> Outcome) -------------------
cat("--- 模型 1: 治疗对心血管结局的总效应 (Total Effect) ---\n")
model_total_effect <- coxph(
  Surv(observed_time, status) ~ treatment + age + sex,
  data = left_join(event_data, baseline_data, by = "id")
)
print(summary(model_total_effect))
# 提取总效应系数
coef_c <- coef(model_total_effect)["treatment"]
coef_c

# 2.2 拟合模型 2: 估计路径 a (Treatment -> Mediator) ---------------------------
cat("\n--- 模型 2: 治疗对中介变量(LVH)的影响 (Path 'a') ---\n")
lvh_at_year_3 <- long_data %>% filter(visit_time == 3)
model_path_a <- lm(lvh_index ~ treatment + age + sex, data = lvh_at_year_3)
print(summary(model_path_a))

# 提取路径 a 的系数
coef_a <- coef(model_path_a)["treatment"]
coef_a

# 2.3 拟合模型 3: 估计路径 b 和直接效应 c' -------------------------------------
cat("\n--- 模型 3: 估计路径 'b' 和直接效应 'c_prime' ---\n")
model_paths_b_and_c_prime <- coxph(
  Surv(tstart, tstop, event) ~ treatment + lvh + age + sex,
  data = final_analysis_data
)
print(summary(model_paths_b_and_c_prime))

# 提取路径 b 的系数 (LVH 对结局的影响)
coef_b <- coef(model_paths_b_and_c_prime)["lvh"]
coef_b
# 提取直接效应 c' 的系数 (控制LVH后，治疗对结局的剩余影响)
coef_c_prime <- coef(model_paths_b_and_c_prime)["treatment"]
coef_c_prime


#  3.1 使用系数乘积法计算 ------------------------------------------------------
indirect_effect <- coef_a * coef_b
mediation_proportion <- indirect_effect / coef_c
mediation_proportion


# 第二种方法：手动bootstrap-----------------------------------------------------------------
# ==============================================================================
# 改进版: 使用手动Bootstrap对时变中介变量进行中介分析
# 采用"差值法"(c - c')计算中介比例，更适合Cox回归框架
# ==============================================================================

# *** 重要方法学说明 ***
# 传统的"乘积法"(a×b)在时变中介+Cox回归框架下存在以下问题：
# 1. HR不具备可加性：总效应c ≠ 直接效应c' + 间接效应a×b
# 2. 量纲不匹配：时变斜率(a) × 瞬时风险系数(b)缺乏明确的时间积分意义
# 3. 因果解释困难：在HR尺度下，乘积法的中介比例难以解释
#
# 本脚本采用的"差值法"优势：
# ✓ 基于同一Cox框架：总效应c和直接效应c'都来自Cox模型
# ✓ 逻辑一致性：间接效应 = c - c'，中介比例 = (c-c')/c
# ✓ 因果可解释：在RCT设定下，差值法具有较好的因果含义
# ✓ 稳健性：Bootstrap提供非参数置信区间，无需正态性假设
# ✓ 时变处理：使用cluster()调整个体内相关性，处理时变协变量更合理

# --- 步骤 0: 准备工作 - 加载额外包 ---
# lme4 包用于拟合线性混合效应模型
# 如果尚未安装，请取消下面代码的注释并运行
# install.packages("lme4")
library(lme4)

# 之前已加载 survival 和 dplyr

# --- 步骤 1: 设置Bootstrap参数 ---
n_boot <- 500 # Bootstrap重复次数。建议至少1000次，此处为演示设为500
set.seed(42) # 保证结果可重复

# 创建向量来存储每次Bootstrap的结果
total_effects_boot <- numeric(n_boot)      # 总效应 c
direct_effects_boot <- numeric(n_boot)     # 直接效应 c'
indirect_effects_boot <- numeric(n_boot)   # 间接效应 c - c'
mediation_proportions_boot <- numeric(n_boot) # 中介比例 (c - c')/c

# 记录成功的迭代次数
successful_iterations <- 0

# --- 步骤 2: 执行Bootstrap循环 ---
# 这个过程可能需要几分钟，取决于您的电脑性能和n_boot的大小
cat(paste0("开始执行 ", n_boot, " 次Bootstrap模拟...\n"))
cat("使用差值法 (c - c') 计算中介效应，更适合Cox回归框架\n\n")

for (i in 1:n_boot) {
  # 使用tryCatch处理可能的模型拟合失败
  tryCatch({
    # --- 2.1 对研究对象(id)进行有放回的重抽样 ---
    boot_ids <- sample(unique(baseline_data$id), size = n_patients, replace = TRUE)
    
    # 创建自助样本，需要从原始数据中提取对应id的所有记录
    boot_baseline_data <- baseline_data[match(boot_ids, baseline_data$id), ]
    # 需要处理重复ID的问题，为自助样本创建新的唯一ID
    boot_baseline_data$new_id <- 1:n_patients
    
    # 对应地创建自助长格式数据和生存数据
    boot_long_data <- long_data %>% filter(id %in% boot_ids) %>%
      left_join(dplyr::select(boot_baseline_data, id, new_id), by="id")
    
    boot_event_data <- event_data %>% filter(id %in% boot_ids) %>%
      left_join(dplyr::select(boot_baseline_data, id, new_id), by="id")
    
    # --- 2.2 估计总效应 c (不含时变LVH的Cox模型) ---
    boot_total_df <- dplyr::left_join(
      boot_event_data, 
      dplyr::select(boot_baseline_data, new_id, treatment, age, sex),
      by = "new_id"
    )
    
    model_total_boot <- coxph(Surv(observed_time, status) ~ treatment + age + sex,
                              data = boot_total_df)
    coef_c_boot <- coef(model_total_boot)["treatment"]
    
    # --- 2.3 估计直接效应 c' (含时变LVH的Cox模型) ---
    # 准备时变数据格式
    boot_final_data <- tmerge(
      data1 = dplyr::select(boot_baseline_data, new_id, treatment, age, sex),
      data2 = dplyr::select(boot_event_data, new_id, observed_time, status),
      id = new_id,
      event = event(observed_time, status)
    )
    boot_final_data <- tmerge(
      data1 = boot_final_data,
      data2 = dplyr::select(boot_long_data, new_id, visit_time, lvh_index),
      id = new_id,
      lvh = tdc(visit_time, lvh_index)
    )
    
    # 使用cluster(new_id)获得稳健标准误，因为每个个体被拆成多条记录
    model_direct_boot <- coxph(Surv(tstart, tstop, event) ~ treatment + lvh + age + sex + cluster(new_id), 
                               data = boot_final_data)
    coef_c_prime_boot <- coef(model_direct_boot)["treatment"]
    
    # --- 2.4 基于log-HR的"差值法"计算中介效应 ---
    indirect_boot <- as.numeric(coef_c_boot - coef_c_prime_boot)  # 间接效应 = c - c'
    
    # 计算中介比例，处理分母接近0的情况
    if (abs(coef_c_boot) > 1e-6) {  # 避免除以接近0的数
      pm_boot <- as.numeric(indirect_boot / coef_c_boot)  # 中介比例 = (c - c')/c
    } else {
      pm_boot <- NA  # 总效应接近0时，中介比例无意义
    }
    
    # --- 2.5 存储本次迭代的结果 ---
    total_effects_boot[i] <- as.numeric(coef_c_boot)
    direct_effects_boot[i] <- as.numeric(coef_c_prime_boot)
    indirect_effects_boot[i] <- indirect_boot
    mediation_proportions_boot[i] <- pm_boot
    
    successful_iterations <- successful_iterations + 1
    
  }, error = function(e) {
    # 如果模型拟合失败，记录为NA
    total_effects_boot[i] <- NA
    direct_effects_boot[i] <- NA
    indirect_effects_boot[i] <- NA
    mediation_proportions_boot[i] <- NA
    
    # 可选：打印错误信息（用于调试）
    # cat(paste0("迭代 ", i, " 失败: ", e$message, "\n"))
  })
  
  # 打印进度
  if (i %% 50 == 0) {
    cat(paste0("已完成 ", i, "/", n_boot, " 次迭代... (成功: ", successful_iterations, ")\n"))
  }
}

cat("Bootstrap模拟完成!\n")
cat(paste0("成功完成的迭代次数: ", successful_iterations, "/", n_boot, "\n\n"))

# --- 步骤 3: 分析Bootstrap结果 ---
cat("\n=================================================================\n")
cat("          改进版Bootstrap中介分析 - 最终结果\n")
cat("          使用差值法 (c - c') 计算中介效应\n")
cat("=================================================================\n\n")

# 移除失败的迭代（NA值）
valid_indices <- !is.na(mediation_proportions_boot) & !is.na(indirect_effects_boot)
valid_total <- total_effects_boot[valid_indices]
valid_direct <- direct_effects_boot[valid_indices]
valid_indirect <- indirect_effects_boot[valid_indices]
valid_pm <- mediation_proportions_boot[valid_indices]

cat(paste0("有效的Bootstrap样本数: ", sum(valid_indices), "/", n_boot, "\n\n"))

# === 1. 总效应 (c) 的Bootstrap结果 ===
cat("--- 1. 总效应 (c: Treatment → Outcome) ---\n")
total_point <- median(valid_total, na.rm = TRUE)
total_ci_lower <- quantile(valid_total, 0.025, na.rm = TRUE)
total_ci_upper <- quantile(valid_total, 0.975, na.rm = TRUE)

cat(paste0("总效应 log-HR (中位数): ", round(total_point, 5), "\n"))
cat(paste0("95% Bootstrap置信区间: [", round(total_ci_lower, 5), ", ", round(total_ci_upper, 5), "]\n"))
cat(paste0("总效应 HR (中位数): ", round(exp(total_point), 3), "\n"))
cat(paste0("总效应 HR 95%CI: [", round(exp(total_ci_lower), 3), ", ", round(exp(total_ci_upper), 3), "]\n\n"))

# === 2. 直接效应 (c') 的Bootstrap结果 ===
cat("--- 2. 直接效应 (c': Treatment → Outcome | LVH) ---\n")
direct_point <- median(valid_direct, na.rm = TRUE)
direct_ci_lower <- quantile(valid_direct, 0.025, na.rm = TRUE)
direct_ci_upper <- quantile(valid_direct, 0.975, na.rm = TRUE)

cat(paste0("直接效应 log-HR (中位数): ", round(direct_point, 5), "\n"))
cat(paste0("95% Bootstrap置信区间: [", round(direct_ci_lower, 5), ", ", round(direct_ci_upper, 5), "]\n"))
cat(paste0("直接效应 HR (中位数): ", round(exp(direct_point), 3), "\n"))
cat(paste0("直接效应 HR 95%CI: [", round(exp(direct_ci_lower), 3), ", ", round(exp(direct_ci_upper), 3), "]\n\n"))

# === 3. 间接效应 (c - c') 的Bootstrap结果 ===
cat("--- 3. 间接效应 (c - c': 通过LVH的中介效应) ---\n")
indirect_point <- median(valid_indirect, na.rm = TRUE)
indirect_ci_lower <- quantile(valid_indirect, 0.025, na.rm = TRUE)
indirect_ci_upper <- quantile(valid_indirect, 0.975, na.rm = TRUE)

cat(paste0("间接效应 log-HR (中位数): ", round(indirect_point, 5), "\n"))
cat(paste0("95% Bootstrap置信区间: [", round(indirect_ci_lower, 5), ", ", round(indirect_ci_upper, 5), "]\n\n"))

# === 4. 中介比例 (PM) 的Bootstrap结果 ===
cat("--- 4. 中介比例 (PM = (c - c')/c) ---\n")
pm_point <- median(valid_pm, na.rm = TRUE)
pm_ci_lower <- quantile(valid_pm, 0.025, na.rm = TRUE)
pm_ci_upper <- quantile(valid_pm, 0.975, na.rm = TRUE)

cat(paste0("中介比例 (中位数): ", round(pm_point, 3), " (", round(pm_point * 100, 1), "%)\n"))
cat(paste0("95% Bootstrap置信区间: [", round(pm_ci_lower, 3), ", ", round(pm_ci_upper, 3), "]"))
cat(paste0(" ([", round(pm_ci_lower * 100, 1), "%, ", round(pm_ci_upper * 100, 1), "%])\n\n"))

# === 5. 统计学显著性判断 ===
cat("--- 5. 统计学显著性评估 ---\n")

# 间接效应的显著性
if (indirect_ci_lower * indirect_ci_upper > 0) {
  cat("✓ 间接效应具有统计学显著性 (95%CI不包含0)\n")
  indirect_significant <- TRUE
} else {
  cat("✗ 间接效应没有统计学显著性 (95%CI包含0)\n")
  indirect_significant <- FALSE
}

# 中介比例的显著性（通常看间接效应是否显著）
if (indirect_significant) {
  cat("✓ 中介效应具有统计学意义\n")
  cat(paste0("  治疗效应中约有 ", round(abs(pm_point) * 100, 1), "% 是通过LVH变化实现的\n"))
} else {
  cat("✗ 中介效应没有统计学意义\n")
  cat("  没有足够证据表明LVH在治疗效应中起中介作用\n")
}

cat("\n--- 6. 方法学说明 ---\n")
cat("本分析采用'差值法'计算中介效应，相比传统'乘积法'具有以下优势:\n")
cat("• 所有效应估计都基于同一Cox回归框架，保证了方法的一致性\n")
cat("• 避免了HR尺度下'乘积法'的理论问题（HR不具备可加性）\n")
cat("• 时变协变量的处理更加合理，使用cluster()调整了个体内相关性\n")
cat("• Bootstrap方法提供了稳健的置信区间，无需正态性假设\n")
cat("• 差值法 (c - c') 在RCT设定下具有较好的因果解释性\n\n")

# === 7. 结果解释指南 ===
cat("--- 7. 结果解释 ---\n")
if (indirect_significant) {
  if (pm_point > 0) {
    cat("解释: 强化治疗通过降低LVH指标，进而降低了心血管事件风险。\n")
    cat(paste0("LVH的改善解释了治疗总效应的 ", round(pm_point * 100, 1), "%。\n"))
  } else {
    cat("解释: 存在统计学显著的中介效应，但方向需要进一步解释。\n")
  }
} else {
  cat("解释: 虽然治疗可能对LVH和心血管结局都有影响，\n")
  cat("但没有足够证据表明LVH变化是治疗效应的重要中介路径。\n")
  cat("治疗可能通过其他机制（如血压、血管功能等）发挥作用。\n")
}

cat("\n=================================================================\n")
cat("                    分析完成\n")
cat("=================================================================\n")

# 第二种方法是更可靠的，第一种方法是非常不严谨的

# ==============================================================================
# 总结与建议
# ==============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("                    脚本总结与建议\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("本脚本实现了两种中介分析方法的对比：\n\n")

cat("【方法1: 传统系数乘积法】\n")
cat("• 使用 a×b 计算间接效应，其中：\n")
cat("  - a: 治疗对第3年LVH的影响（线性回归）\n") 
cat("  - b: 时变LVH对结局的影响（Cox回归）\n")
cat("• 问题：量纲不匹配，HR尺度下缺乏严格的因果解释\n")
cat("• 结论：仅作为探索性分析，不建议作为主要结果报告\n\n")

cat("【方法2: 改进的差值法 + Bootstrap】★ 推荐\n")
cat("• 使用 c - c' 计算间接效应，其中：\n")
cat("  - c: 总效应（不含LVH的Cox模型）\n")
cat("  - c': 直接效应（含时变LVH的Cox模型）\n")
cat("• 优势：方法一致性好，因果解释清晰，统计推断稳健\n")
cat("• 结论：适合作为主要分析结果，符合现代中介分析标准\n\n")

cat("【实际应用建议】\n")
cat("1. 论文报告：以方法2的结果为主，方法1可作为敏感性分析\n")
cat("2. Bootstrap次数：正式分析建议使用1000-5000次\n")
cat("3. 模型诊断：检查Cox模型的比例风险假设\n")
cat("4. 敏感性分析：\n")
cat("   - 尝试不同的LVH测量时点\n")
cat("   - 考虑非线性关系（样条函数）\n")
cat("   - 评估未测量混杂因子的影响\n\n")

cat("【进一步改进方向】\n")
cat("• G-formula方法：在特定时点计算风险差的自然间接效应\n")
cat("• 联合建模：使用JM包处理纵向数据与生存数据的联合分析\n")
cat("• 因果中介分析：结合倾向性评分或工具变量方法\n\n")

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("脚本执行完毕。建议重点关注方法2的Bootstrap中介分析结果。\n")
cat(paste(rep("=", 70), collapse = ""), "\n")


