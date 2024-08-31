
#联合诊断ROC-------------------
rm(list=ls())
options(stringsAsFactors = F)
setwd("F:/CRC/CRC/TCGA/")
load("TCGA_CMS.rda")

setwd("F:/CRC/GEO/GSE72970/")
load("GSE72970.rda")
samplefile <- read.table("TCGACMS4_123SampleFile.txt", header=TRUE, row.names=NULL, quote="", sep="\t", check.names = FALSE)
key <- read.table("F:/CRC/NOTCH3/33gene.txt", header=TRUE, row.names=NULL, quote="", sep="\t", check.names = FALSE)
CRC =new_exprSet#[,-1]
Inter_exp =as.data.frame(CRC)
Inter_exp = as.data.frame(Inter_exp[, match(samplefile$Sample_Name, colnames(Inter_exp))])

Inter_exp =Inter_exp[key$gene,]
Inter_exp =as.data.frame(t(Inter_exp))
combine<-cbind(samplefile,Inter_exp)
combine=combine[,-1]
colnames(combine)[1]="outcome"

write.table(combine, file = "TCGA_33ROCcombine.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = T)
#二分类变量改为1和0，读进来
combine <- read.table("TCGA_33ROCcombine.txt", header=TRUE, row.names=1, quote="", sep="\t", check.names = FALSE)
library(ggplot2) 
library(RColorBrewer) 
library(ggpubr)
library(preprocessCore)
library(e1071)
library(parallel)
library(data.table)
library(pROC)
library(glmnet)
#C1S+C1R+C3 TIMP2+AEBP1+PLAT+SPARC+IGFBP6+SERPINE1
fit1 <- glm(outcome ~NOTCH3+COL6A2+MRC2+COL5A1+COL1A2+COL6A1+
              THBS2+COL6A3+FBN1+FGFR1+HEG1+SPARC+C1R+FN1+AOC3+ISLR+SERPING1+
              COL15A1+FBLN5+DDR2+C1S+TAGLN+MFAP4+VIM+LMOD1+C3+
              GFPT2+CTSK+HSPB8+THBS1+IGFBP6+AEBP1+SERPINE1 ,
            data=combine,
            family = binomial())  
summary(fit1)
combine$prob <- predict(fit1, 
                        newdata=combine, 
                        type="response")
head(combine)
write.table(combine, file = "TCGA_33ROC.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = T)
roc1 <- roc(combine$outcome, combine$NOTCH3) 
roc2 <- roc(combine$outcome, combine$COL6A2)
roc3 <- roc(combine$outcome, combine$MRC2) 
roc4 <- roc(combine$outcome, combine$COL5A1)
roc5 <- roc(combine$outcome, combine$COL1A2) 
roc6 <- roc(combine$outcome, combine$COL6A1)
roc7 <- roc(combine$outcome, combine$THBS2)
roc8 <- roc(combine$outcome, combine$COL6A3) 
roc9 <- roc(combine$outcome, combine$FBN1)
roc10 <- roc(combine$outcome, combine$FGFR1) 
roc11 <- roc(combine$outcome, combine$HEG1)
roc12 <- roc(combine$outcome, combine$SPARC) 
roc13 <- roc(combine$outcome, combine$C1R)
roc14 <- roc(combine$outcome, combine$FN1)
roc15 <- roc(combine$outcome, combine$AOC3) 
roc16 <- roc(combine$outcome, combine$ISLR)
roc17 <- roc(combine$outcome, combine$SERPING1) 
roc18 <- roc(combine$outcome, combine$COL15A1) 
roc19 <- roc(combine$outcome, combine$FBLN5)
roc20 <- roc(combine$outcome, combine$DDR2) 
roc21 <- roc(combine$outcome, combine$C1S)
roc22 <- roc(combine$outcome, combine$TAGLN) 
roc23 <- roc(combine$outcome, combine$MFAP4)
roc24 <- roc(combine$outcome, combine$VIM)
roc25 <- roc(combine$outcome, combine$LMOD1) 
roc26 <- roc(combine$outcome, combine$C3)
roc27 <- roc(combine$outcome, combine$GFPT2) 
roc28 <- roc(combine$outcome, combine$CTSK)
roc29 <- roc(combine$outcome, combine$HSPB8) 
roc30 <- roc(combine$outcome, combine$THBS1)
roc31 <- roc(combine$outcome, combine$IGFBP6)
roc32 <- roc(combine$outcome, combine$AEBP1)
roc33 <- roc(combine$outcome, combine$SERPINE1) 
roc34 <- roc(combine$outcome, combine$prob)

roc1;roc2;roc3;roc4;roc5;roc6;roc7;roc8;roc9;roc10;roc11;roc12;roc13;roc14;roc15;roc16;roc17;roc18;roc19;roc20;roc21;roc22;roc23;roc24;roc25;roc26;roc27;roc28;roc29;roc30;roc31;roc32;roc33;roc34;
#tmp1 = as.data.frame(roc1$sensitivities) 
plot(roc1,  # 前面构建的ROC对象
     print.auc=TRUE, # 图像上输出AUC值
     print.auc.x=0.5, print.auc.y=0.85, # AUC值坐标为（x，y）
     auc.polygon=TRUE, # 将ROC曲线下面积转化为多边形
     auc.polygon.col="#DDDDDD",  # 设置多边形的填充颜色
     grid= FALSE, # 不显示网格背景线
     legacy.axes=TRUE)  # 使横轴从0到1，表示为1-特异度
plot.roc(roc2, add=TRUE,  # 增加曲线
         col = "#A500CC", # 设置曲线颜色
         print.auc=TRUE, # 图像上输出AUC
         print.auc.col = "#A500CC", # 设置AUC文本的颜色
         print.auc.x=0.5, print.auc.y=0.8)  # AUC的坐标为（x，y）
plot.roc(roc3, add=TRUE,  # 增加曲线
         col = "#FFFFBB", # 设置曲线颜色
         print.auc=TRUE, # 图像上输出AUC
         print.auc.col = "#FFFFBB", # 设置AUC文本的颜色
         print.auc.x=0.5, print.auc.y=0.75)  # AUC的坐标为（x，y）
plot.roc(roc4, add=TRUE,  # 增加曲线
         col = "#CCFF99", # 设置曲线颜色  #B8860B
         print.auc=TRUE, # 图像上输出AUC
         print.auc.col = "#CCFF99", # 设置AUC文本的颜色
         print.auc.x=0.5, print.auc.y=0.7)  # AUC的坐标为（x，y）
plot.roc(roc5, add=TRUE,  # 增加曲线
         col = "#00BBFF", # 设置曲线颜色
         print.auc=TRUE, # 图像上输出AUC
         print.auc.col = "#00BBFF", # 设置AUC文本的颜色
         print.auc.x=0.5, print.auc.y=0.65)  # AUC的坐标为（x，y）
plot.roc(roc6, add=TRUE,  # 增加曲线
         col = "#6A5ACD", # 设置曲线颜色
         print.auc=TRUE, # 图像上输出AUC
         print.auc.col = "#6A5ACD", # 设置AUC文本的颜色
         print.auc.x=0.5, print.auc.y=0.6)  # AUC的坐标为（x，y）
plot.roc(roc7, add=TRUE,  # 增加曲线
         col = "#99FFFF", # 设置曲线颜色
         print.auc=TRUE, # 图像上输出AUC
         print.auc.col = "#99FFFF", # 设置AUC文本的颜色
         print.auc.x=0.5, print.auc.y=0.55)  # AUC的坐标为（x，y）
plot.roc(roc8, add=TRUE,  # 增加曲线
         col = "#5599FF", # 设置曲线颜色
         print.auc=TRUE, # 图像上输出AUC
         print.auc.col = "#5599FF", # 设置AUC文本的颜色
         print.auc.x=0.5, print.auc.y=0.5)  # AUC的坐标为（x，y）
plot.roc(roc9, add=TRUE,  # 增加曲线
         col = "#32CD32", # 设置曲线颜色
         print.auc=TRUE, # 图像上输出AUC
         print.auc.col = "#32CD32", # 设置AUC文本的颜色
         print.auc.x=0.5, print.auc.y=0.45)  # AUC的坐标为（x，y）
plot.roc(roc10, add=TRUE,  # 增加曲线
         col = "#FFB3FF", # 设置曲线颜色  #B8860B
         print.auc=TRUE, # 图像上输出AUC
         print.auc.col = "#FFB3FF", # 设置AUC文本的颜色
         print.auc.x=0.5, print.auc.y=0.4)  # AUC的坐标为（x，y）