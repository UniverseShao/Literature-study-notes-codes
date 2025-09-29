## 安装常用包
# install.packages(c("tidyverse","WeightIt","cobalt","survey","survival","cmprsk"))

library(tidyverse)
library(WeightIt)   # 直接给 OW 权重
library(cobalt)     # 平衡诊断
library(survival)   # KM/Cox
library(cmprsk)     # Fine-Gray
library(survey)     # 设计对象 + 稳健方差（可选）

## 数据：df，变量：
## treat: 1=治疗 0=对照
## time: 随访时间
## status: 事件指示 (生存分析 1=事件 0=删失)
## fstatus: 竞争风险 (0=删失,1=目标事件,2=竞争事件)
## 协变量：age, sex, x1, x2, ...

## 1) 倾向性评分 + 重叠权重
w.out <- weightit(treat ~ age + sex + x1 + x2,
                  data = df, method = "overlap")  # OW
df$w <- w.out$weights
df$ps <- w.out$ps

## 2) 配平诊断
bal.tab(w.out, un = TRUE, m.threshold = 0.1)   # SMD表
love.plot(w.out, abs = TRUE, thresholds = c(m = .1))

## 3) 加权KM（2年累计风险与RD）
sf <- survfit(Surv(time, status) ~ treat, data = df, weights = w)
# 提取2年时点的累计风险，可用 summary(sf, times=2*365.25)
km2 <- summary(sf, times = 2*365.25)
# 组别顺序: treat=0,1；计算 RD = F1 - F0
F0 <- km2$cumhaz[1]; F1 <- km2$cumhaz[2]  # 若使用 cumulative hazard 需转化；更稳妥请用 km2$surv
# 建议：用 bootstrap 计算 2年RD 及 95%CI（略）

## 4) 加权 Cox（HR）
fit.cox <- coxph(Surv(time, status) ~ treat,
                 data = df, weights = w, robust = TRUE)
summary(fit.cox)

## 5) 竞争风险 Fine-Gray（SHR）
fit.fg <- crr(ftime = df$time, fstatus = df$fstatus,
              cov1 = model.matrix(~ treat, df)[, -1, drop=FALSE],
              weights = df$w, robust = TRUE)
fit.fg

## 6)（可选）双重稳健（增广 OW）思路
# 在加权基础上再建一个 outcome model（如 Cox: ~ age+sex+...），
# 用 AIPW/augmented 公式或现成包（如 PSweight）得到 DR-ATO 估计与方差。
