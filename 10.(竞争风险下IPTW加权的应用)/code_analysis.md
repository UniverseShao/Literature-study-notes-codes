# IPTW竞争风险分析R代码详细解析

## 代码概述

这个R代码实现了基于逆概率治疗加权(IPTW)和竞争风险分析的统计方法，用于评估治疗效果。代码主要包含以下几个部分：

1. 数据准备和预处理
2. 限制性立方样条(RCS)变换
3. 加权Aalen-Johansen估计器
4. IPTW-IPCW估计器
5. Bootstrap置信区间估计

## 详细逐行解析

### 1. 库加载和初始化

```r
library(survival)      # 生存分析基础包
library(cmprsk)        # 竞争风险分析
library(riskRegression) # 风险回归分析
library(prodlim)       # 产品限估计器
library(rms)           # 回归建模策略

N.boot <- 2000         # Bootstrap迭代次数
```

### 2. 数据读取和准备

```r
# 定义变量列表和类型
zlist <- list(age=0,prevMI=0,prevCHF=0,HeartRate=0,SBP=0,
              creatinine=0,enzymes=0,stemi=0,pci=0,treat=0,
              event.type=0,T.event=0)

# 读取数据文件
effect.df <- data.frame(scan("cohort.dat",zlist))
```

**变量说明：**
- `age`: 年龄
- `prevMI`: 既往心肌梗死史
- `prevCHF`: 既往充血性心力衰竭史
- `HeartRate`: 心率
- `SBP`: 收缩压
- `creatinine`: 肌酐水平
- `enzymes`: 酶水平
- `stemi`: ST段抬高型心肌梗死
- `pci`: 经皮冠状动脉介入治疗
- `treat`: 治疗指示变量
- `event.type`: 事件类型
- `T.event`: 事件发生时间

### 3. 限制性立方样条(RCS)变换

```r
# 为连续变量创建RCS项，用于非线性关系建模
age.mat <- rcspline.eval(effect.df$age,nk=5,inclx=F,norm=0)
HeartRate.mat <- rcspline.eval(effect.df$HeartRate,nk=5,inclx=F,norm=0)
SBP.mat <- rcspline.eval(effect.df$SBP,nk=5,inclx=F,norm=0)
creatinine.mat <- rcspline.eval(effect.df$creatinine,nk=5,inclx=F,norm=0)

# 设置列名
colnames(age.mat) <- c("age2","age3","age4")
colnames(HeartRate.mat) <- c("HeartRate2","HeartRate3","HeartRate4")
...

# 合并到主数据集
effect.df <- cbind(effect.df,age.mat.df,HeartRate.mat.df,SBP.mat.df,
                   creatinine.mat.df)
```

**RCS参数说明：**
- `nk=5`: 节点数，创建4个样条项
- `inclx=F`: 不包含原始变量
- `norm=0`: 不使用标准化

### 4. 倾向评分模型和IPTW权重计算

```r
# 构建倾向评分模型
psm <- glm(treat ~ age + prevMI + prevCHF + HeartRate + SBP + creatinine +
           enzymes + stemi + pci +  # 基础协变量
           age2 + age3 + age4 + HeartRate2 + HeartRate3 + HeartRate4 +  # RCS项
           SBP2 + SBP3 + SBP4 + creatinine2 + creatinine3 + creatinine4,
           family="binomial",data=effect.df)

# 计算倾向评分
ps <- psm$fitted
Z <- effect.df$treat

# 计算IPTW权重
effect.df$iptw <- (Z/ps) + ((1-Z)/(1-ps))
```

**IPTW权重公式：**
- 治疗组权重: 1/ps
- 对照组权重: 1/(1-ps)
- 稳定化权重: (Z/ps) + ((1-Z)/(1-ps))

### 5. 加权Aalen-Johansen估计器

```r
# 使用加权方法估计累积风险函数
cif.aj <- prodlim(data=effect.df,Hist(T.event,event.type) ~ treat,
                  type="risk",caseweights=effect.df$iptw)

# 预测1-5年的风险
risk.aj <- predict(cif.aj,times=365*(1:5),cause='1',
                   newdata=data.frame(treat=0:1))

# 计算风险差(RD)和相对风险(RR)
risk0.aj <- risk.aj$'treat=0'
risk1.aj <- risk.aj$'treat=1'
rd.aj <- risk1.aj - risk0.aj
rr.aj <- risk1.aj/risk0.aj
```

### 6. IPTW-IPCW估计器

```r
# 构建治疗模型
treat_mod <- glm(treat.factor ~ ..., family="binomial",data=effect.df)

# 构建竞争风险模型
cox_mod <- CSC(Hist(T.event,event.type) ~ ..., data=effect.df)

# 构建删失模型
censor_mod <- coxph(Surv(T.event,event.type==0) ~ 1,data=effect.df)

# 计算ATE估计
ate.iptw <- ate(event=cox_mod,treatment=treat_mod,censor=censor_mod,
                data=effect.df,estimator="IPTW",times=365*(1:5),cause=1,se=F)

# 计算相对风险
rr.iptw <- ate.iptw$diffRisk$estimate.B/ate.iptw$diffRisk$estimate.A
```

### 7. Bootstrap置信区间估计

```r
# 初始化存储向量
rd.aj.bs <- NULL
rr.aj.bs <- NULL
rr.iptw.bs <- NULL

# Bootstrap循环
for (iter in 1:N.boot){
    # 设置随机种子保证可重现性
    set.seed(iter)
    
    # 重采样
    sample.id <- sample(1:nrow(effect.df),size=nrow(effect.df),replace=T)
    boot.df <- effect.df[sample.id,]
    
    # 在bootstrap样本中重新执行所有分析步骤...
    # (重复上述所有步骤)
    
    # 存储结果
    rd.aj.bs <- rbind(rd.aj.bs,rd.aj)
    rr.aj.bs <- rbind(rr.aj.bs,rr.aj)
    rr.iptw.bs <- rbind(rr.iptw.bs,rr.iptw)
}

# 计算置信区间
rd.aj.ci.lower <- apply(rd.aj.bs,MARGIN=2,FUN=quantile,probs=0.025)
rd.aj.ci.upper <- apply(rd.aj.bs,MARGIN=2,FUN=quantile,probs=0.975)
...
```

## 统计方法总结

### 1. 逆概率治疗加权(IPTW)
- **目的**: 调整治疗选择偏倚
- **方法**: 基于倾向评分对观察单位进行加权
- **假设**: 无未测量的混杂因素

### 2. 竞争风险分析
- **目的**: 处理多个互斥终点事件
- **方法**: 使用Aalen-Johansen估计器
- **优势**: 避免传统生存分析中的独立删失假设

### 3. Bootstrap方法
- **目的**: 估计置信区间
- **方法**: 重采样技术
- **优势**: 不依赖正态分布假设

### 4. 限制性立方样条(RCS)
- **目的**: 建模连续协变量的非线性效应
- **方法**: 分段多项式函数
- **优势**: 灵活性和稳定性平衡

## 输出结果

代码生成以下输出文件：
- `boot.rd.aj.out`: Aalen-Johansen估计的风险差
- `boot.rr.aj.out`: Aalen-Johansen估计的相对风险
- `boot.rr.iptw.out`: IPTW-IPCW估计的相对风险
- `boot.iter`: Bootstrap迭代进度

每个输出文件包含点估计和95%置信区间。