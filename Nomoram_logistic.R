# 常用命令 ----------------------------------------------------------------
getwd()#显示默认路径
setwd ("C:/Users/95654/Documents") #设置默认路径，这是我的文档默认路径
rm(list=ls())#rm(list=ls())清除所有内存数据
gc()#垃圾回收数据
mydata <- as.data.frame(mydata)
write.table(mydata, file = "mydata.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
#数据打包 
dd = datadist(mydata_copy)
option <- options(datadist = "dd")
mydata <- mydata[,-7]
#删除第6列--------------------
# 下载所需要的R包 ----------------------------------------------------------------
# 检查 VIM 包是否已经安装
if (!require("VIM", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("VIM")
}

# 检查 rms 包是否已经安装
if (!require("rms", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("rms")
}

# 检查 nomogramFormula 包是否已经安装
if (!require("nomogramFormula", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("nomogramFormula")
}

# 检查 pROC 包是否已经安装
if (!require("pROC", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("pROC")
}

# 检查 rmda 包是否已经安装
if (!require("rmda", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("rmda")
}
install.packages("MASS")
install.packages( "ggplot2")

# 加载需要的R包 -----------------------------------------------------------------
library(readxl)#读取xlsx文件
library(VIM)#缺失值可视化
library(rms)#拟合模型
library(nomogramFormula)#计算列线图得分
library(pROC)#绘制ROC曲线，计算AUC和95%置信区间
library(rmda)#临床决策曲线和临床影响曲线
library(MASS)#逐步回归
library(regplot)#添加个案例列线图
library(DynNom)#动态列线图
library(ggplot2)
# 载入示例数据 ------------------------------------------------------------------
mydata <- read.csv("lasso_data.csv",row.names = 1,sep = ",",header = T)
# 或者文件是xlsx文件
mydata <- read_excel("501例分析确定数据7.xlsx",sheet="Sheet2")

mydata_copy <- mydata#数据备份
# 查看数据和简单数据清洗 -------------------------------------------------------------
# 查看数据中的缺失值
aggr(mydata, col=c('skyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(mydata), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Percentage"))
# 删除缺失值所在的行
mydata <- na.omit(mydata)
summary(mydata)
# 变量的转换与赋值 ------------------------------------------------------------
#生存变量的处理
table(mydata$stage)#检查生存状态
mydata$stage <- factor(mydata$stage,level=c(0,1),labels=c("Survival","Dead"))
#最终变量不宜过多，一般 变量数＜阳性时间/10-20

#二分类变量的转换（转换成因子，列线图能展示相应分组字符，计算评分时不要运行）---------------------------------------------------------------------------------
mydata$heart<- factor(mydata$heart ,level=c(0,1),labels=c("no disease","disease"))
mydata$gender <- factor(mydata$gender ,level=c(0,1),labels=c("male","female"))

#保留第一列，第3列，第5列到第七列
mydata <- mydata[,c(1,3,5:7)]
#删除第6列，
mydata <- mydata[,-6]

# 开始构建logistic回归模型 ------------------------------------------------------------
#数据打包 
dd = datadist(mydata)
option <- options(datadist = "dd")

#拟合模型
colnames(mydata)#查看列名，选择你要构建模型的变量
formula <- as.formula(stage ~  heart+gender + age + postoperativeMLR )
model <- lrm(formula, # 回归模型的公式,指定自变量和因变量
           data = mydata, # 包含所有变量的数据框
           x=TRUE, # logistic回归也称为"广义线性模型",这个参数指定响应变量的二分类
           y=TRUE) # 参数y也是指定因变量是二分类的
model#查看模型具体情况
#1. 模型概况:显示数据集obs数(152)、LR卡方值(29.32)和p值(0.0001),说明模型整体是显著的。R2值为0.239,表现一定的预测能力。
#2. 区分度指数:包括R2(决定系数)、Brier分数(预测精度)、gamma(区分度)、tau-a(区分度)等,用于评估模型的预测能力。值越大越好。
#3. 变量系数:每一列代表一个自变量,从左至右是:
#- Coef:变量的Logistic回归系数,代表变量的权重大小。正值意味着当变量增加时,1的对数几率增加,负值相反。
#- S.E.:系数的标准误,用于后续的Wald检验。
#- Wald Z:系数显著性的Wald统计量,Z越大绝对值越显著。
#- Pr(>|Z|):对应Z值的p值。小于0.05为显著变量。

OR <- exp(model$coefficients)#计算Logistic回归模型中每个自变量的比值比(Odds Ratio, OR)
OR
#逐步回归拟合模型-------------------------------------------------------------------
initial_model <- glm(formula,data = mydata,family=binomial())
summary(initial_model)
#查看模型具体情况，AIC信息准则即Akaike information criterion，是衡量统计模型拟合优良性的一种标准。AIC越小，模型越好，通常选择AIC最小的模型
forward_model <- stepAIC(initial_model,direction = "forward")#向前logistic逐步回归
summary(forward_model)
backward_model <- stepAIC(initial_model,direction ="backward")#向前logistic逐步回归
summary(backward_model)
both_model <- stepAIC(initial_model,direction = "both")#向后logistic逐步回归
summary(both_model)

# 结果可视化 -------------------------------------------------------------------
Nomogram_1 <- nomogram(model,
                       fun = function(x)1/(1+exp(-x)),
                       lp=F,
                       fun.at = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                       funlabel = "Risk")
plot(Nomogram_1)
# 答疑：
# 在Nomogram图的绘制中,线段长度同时取决于两个因素:
# 1) 变量的回归系数:系数绝对值越大的变量,其线段越长,其他条件相同。
# 2) 变量的取值范围:范围越大的变量,其线段也越长,其他条件相同。

#调整一些参数
plot(Nomogram_1,
     #模型的名称
     xfrac = .35,
     #变量与图形的占比（调整变量与坐标抽距离）
     cex.var = 1.6,
     #变量字体加粗
     cex.axis = 1.4,
     #数轴：字体的大小
     tcl = -0.5,
     #数轴：刻度的长度
     lmgp = 0.3,
     #数轴：文字与刻度的距离
     label.every = 1,
     #数轴：刻度下的文字，1=连续显示，2=隔一个显示一个
     col.grid = gray(c(0.8,0.95)))

#计算列线图评分(要计算评分，数据中不能含有因子变量)------------------------------------------------
options(option)
#设置options()函数相关参数
results <- formula_rd (nomogram = Nomogram_1)
#使用formula_rd()函数基于nomogram对象Nomogram_1生成公式formula和相应的数据框rd
mydata$points <- points_cal(formula = results$formula,rd=mydata)
mydata$points 
#使用points_cal()函数基于formula公式和原始数据框mydata计算points列,并添加到mydata
head(mydata) 
#查看mydata的数据框的前6行

#可以根据列线图评分（中位数、均值、四分位数）对所有样本进行分组


# 结果可视化 -------------------------------------------------------------------
Nomogram_1 <- nomogram(model,
                       fun = function(x)1/(1+exp(-x)),
                       lp=F,
                       fun.at = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                       funlabel = "Risk")
plot(Nomogram_1)
# 答疑：
# 在Nomogram图的绘制中,线段长度同时取决于两个因素:
# 1) 变量的回归系数:系数绝对值越大的变量,其线段越长,其他条件相同。
# 2) 变量的取值范围:范围越大的变量,其线段也越长,其他条件相同。

#调整一些参数
plot(Nomogram_1,
     #模型的名称
     xfrac = .35,
     #变量与图形的占比（调整变量与坐标抽距离）
     cex.var = 1.6,
     #变量字体加粗
     cex.axis = 1.0,
     #数轴：字体的大小
     tcl = -0.5,
     #数轴：刻度的长度
     lmgp = 0.3,
     #数轴：文字与刻度的距离
     label.every = 1,
     #数轴：刻度下的文字，1=连续显示，2=隔一个显示一个
     grid = c("#2178E1","#E5572A"))#数轴：刻度的长度
# 添加个案的列线图-------------------------------------------------------------------
Nomogram_2 <- glm(formula, family=binomial(link = logit),data=mydata)

#使用 regplot 函数创建一个逻辑回归图分解，其中包括密度图和箱线图---------------------------------
regplot(Nomogram_2,#Nomogram_2参数指定使用刚才创建的模型
        plots = c("density","boxes"),#plots参数指定要会制的图形类型
        observation = mydata[499,],#指定在图形中展示的样本，选择第几行
        center = T,#center参数指定是否将图形的中心点设置为零，对齐变
        points = TRUE,#points参数指定是否绘制每个观测值的点
        droplines = T,#droplines参数指定是否会制垂直于x轴的线
        title ="individual nomogram",#title参数指定图形的标题
        odds = T,#odds参数指定是否显示比值或者概率
        interval="confidence",#interval参数指定使用置信区间还是预测区
        rank = NULL,#rank参数指定是否对观测值进行排名
        clickable = F,#是否启用交互式模式，这个大家一定要换成T，玩一下
        dexcol ="skyblue" , #dencol参数指定密度图的颜色
        boxcol ="skyblue" ) #boxcol参数指定箱线图的颜色
# 动态列线图-------------------------------------------------------------------
Nomogram_3 <- glm(formula,data = mydata, family = binomial())

DynNom(Nomogram_3,clevel = 0.95,DNtitle ="Nomogram",DNxlab = "probability",DNylab = NULL,DNlimits = NULL)
DynNom(Nomogram_3,#Nomogram_3参数指定使用刚才创建的逻辑回归模型
       clevel = 0.95,#clevel参数指定要使用的置信水平
       DNtitle ="Nomogram",#DNxlab参数指定X轴的标签
       DNxlab = "probability",#DNxlab参数指定X 轴的标签
       DNylab = NULL,#DNylab参数指定 y轴的标签
       DNlimits = NULL)#DNlimits参数指定 X 轴的限制范围
warnings()







# 模型评估 -----------------------------------------------------------------
# ROC曲线 -------------------------------------------------------------------
#AUC计算
model_ROC <- glm(formula,data = mydata,family = binomial())
#lrm()函数属于Design包,专门用于拟合Logistic回归模型。
#glm()函数属于stats包,是用于拟合广义线性模型(Generalized Linear Models)的泛用函数,可以拟合Logistic回归、Poisson回归等多种模型。

#type=“response”给出具体的预测概率，而type=“class”按规定的阈值给出分类
mydata$predvalue <- predict(model_ROC,type="response")
# 计算AUC和95%CI
ROC <- roc(mydata$stage,mydata$predvalue)
auc(ROC)
ci(auc(ROC))
#绘制ROC曲线
plot(ROC, 
     #ROC为ROC对象,用于绘制ROC曲线图
     col="red",  
     #线条颜色为红色
     print.auc=F, 
     #在图中打印AUC值
     print.thres=T,
     #在图中打印最佳切割点的阈值
     lty=1,  
     #线型为实线
     xlab = "Specificity",
     #x轴标签为Specificity,特异度
     ylab = "Sensitivity", 
     #y轴标签为Sensitivity,敏感度
     lwd=2)  #线条宽度为2

legend(0.45,0.05,
       c("AUC of Nomogram:0.6769(95%CI 0.6275-0.7263)"),
       #图例位置及内容
       lty = c(1),   
       #图例线型为实线
       lwd = c(2),   
       #图例线条宽度为2
       col = c("red"), 
       #图例颜色为红色
       bty = "0")
plot(ROC,print.auc=TRUE,auc.polygon=T,max.auc.polygon =F,auc.polygon.col= "skyblue",print.thres= "best")
warnings()

# 假如你想绘制其中几个变量的ROC曲线，并在同一张图中展示 ---------------------------------------------------------------------
colnames(mydata)
roc1 <- roc(mydata$stage,mydata$postoperativeMLR)#模型的ROC
roc2 <- roc(mydata$stage,mydata$age)#年龄变量的ROC
roc3 <- roc(mydata$stage,mydata$gender)#KCNJ13变量ROC
roc4 <- roc(mydata$stage,mydata$heart)#KCNJ13变量ROC
# 计算AUC和95%CI
auc(roc1)#查看曲线下面积auc
ci(auc(roc1))#查看95%置信区间
auc(roc2)#查看曲线下面积auc
ci(auc(roc2))#查看95%置信区间
auc(roc3)#查看曲线下面积auc
ci(auc(roc3))#查看95%置信区间
auc(roc4)#查看曲线下面积auc
ci(auc(roc4))#查看95%置信区间

#方法一:绘制单条ROC曲线,灵敏度与特异度为坐标的总感觉有点问题
plot.new()#创建一个新的绘图窗口
plot(1-roc1$specificities,  
     #x轴为1-Specificity,特异度的补值
     roc1$sensitivities,type = "l", 
     #y轴为ROC$sensitivities,敏感度
     col="red", 
     #线条颜色为红色
     lty=1,
     #线型为实线
     xlab = "1-Specificity",
     #x轴标签为1-Specificity,特异度的补值
     ylab = "Sensitivity",
     #y轴标签为Sensitivity,敏感度
     lwd=2) #线条宽度为2
abline(0,1) #添加对角线
plot(roc1,print.auc=TRUE,auc.polygon=T,max.auc.polygon =F,auc.polygon.col= "skyblue",print.thres= "best")
warnings()

#方法二，显示最佳截断值
plot.new()#创建一个新的绘图窗口
plot(roc, 
     #roc1为ROC对象,用于绘制ROC曲线图
     col="red",  
     #线条颜色为红色
     print.auc=F, 
     #在图中打印AUC值
     print.thres=T,
     #在图中打印最佳切割点的阈值
     lty=1,  
     #线型为实线
     xlab = "Specificity",
     #x轴标签为Specificity,特异度
     ylab = "Sensitivity", 
     #y轴标签为Sensitivity,敏感度
     lwd=2)  
     #线条宽度为2
#增加其他变量
plot(roc1, add=TRUE, col="#D4DD98")
plot(roc2, add=TRUE, col="#99DD98")
plot(roc3, add=TRUE, col="#8F94E0")
plot(roc4, add=TRUE, col="#967345")



#简单计算另外两个变量的AUC和95%CI
auc(roc1)
ci(auc(roc1))
auc(roc2)
ci(auc(roc2))
auc(roc3)
ci(auc(roc3))
auc(roc4)
ci(auc(roc4))
#添加右下脚图例标注
legend("bottomright",legend=c("Nomogram:0.6769 (95% CI: 0.6275-0.7263)",
                              "postoperativeMLR:0.5952 (95% CI: 0.5368-0.6536)",
                              "age:0.6138 （95% CI: 0.5617-0.6659)",
                              "gender：0.56（95% CI: 0.5164-0.6036)",
                              "heart:0.5371 (95% CI: 0.5067-0.5675)"),
       col=c("red","#D4DD98","#99DD98","#8F94E0","#967345"),
       lty=1,
       lwd=3)
#其他方法使用plot.roc和lines.roc绘制组合图；
#percent=TRUE时以百分比的形式展示；
rocobj1<- plot.roc(mydata$术后双下肢静脉血管超声, mydata$术后NLR,
                   main= "Statistical comparison",
                   percent=TRUE, col= "#FF99CC90")
rocobj2<- lines.roc(mydata$术后双下肢静脉血管超声, mydata$术后SII,
                    percent=TRUE, col= "#99CC0090")

# 校准曲线 ---------------------------------------------------------------------
cal <- calibrate(model,method = "boot",B=1000)
# 使用calibrate()函数基于bootstrap方法对模型model进行校准,bootstrap重复次数为1000
plot(cal, 
     #绘制cal校准对象的校准曲线图
     xlim=c(0,1),  
     #x轴范围从0到1
     ylim=c(0,1), 
     #y轴范围从0到1
     xlab="Predicted Probability",
     #x轴标签为Predicted Probability,预测概率
     ylab="Observed Probability", 
     #y轴标签为Observed Probability,观察概率
     subtitles=TRUE) #不显示子标题



# 临床决策曲线和临床影响曲线 -----------------------------------------------------------
model_DCA <- decision_curve(formula, 
                            #公式formula为模型的公式 
                            data = mydata_copy, 
                            #data为模型拟合的数据集 
                            family = binomial(link = 'logit'), 
                            #模型属于二项式Logistic回归  
                            thresholds = seq(0,1,by=0.01),  
                            #阈值范围从0到1,间隔为0.01
                            confidence.intervals = 0.95, 
                            #置信区间为95%
                            study.design = 'case-control', 
                            #研究设计为病例对照研究 
                            population.prevalence = 0.3)  #人群基线发病概率为0.3
#DCA曲线
plot_decision_curve(model_DCA,  
                    #绘制model_DCA的决策曲线图 
                    curve.names = c('Nomogram'), 
                    #曲线名称为Nomogram 
                    xlim = c(0,1.0),   
                    #x轴范围0到1 
                    cost.benefit.axis = TRUE,  
                    #不显示成本效益坐标轴 
                    col = c('red'),   
                    #曲线颜色为红色
                    confidence.intervals = FALSE, 
                    #不显示置信区间 
                    standardize =FALSE ) #不进行标准化FALSE
#临床影响曲线
plot_clinical_impact(model_DCA, 
                     #绘制临床影响曲线图
                     population.size = 1000,  
                     #总人群规模为1000人 
                     cost.benefit.axis = T,    
                     #显示成本效益坐标轴 
                     n.cost.benefits = 8, 
                     #成本效益轴分为8段 
                     col = c('red','blue'), 
                     #红色代表病例,蓝色代表对照 
                     confidence.intervals = F) 




