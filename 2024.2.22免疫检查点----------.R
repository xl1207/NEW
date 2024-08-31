#免疫检查点----------------
rm(list=ls())
options(stringsAsFactors = F)
setwd("F:/CRC/GEO/GSE40967/")

load("GSE40967.rda")
load("NOTCH3_45km.rda")#加载Risk group差异表达数据
#samplefile <- read.table(file = "NOTCH3_45Cluster_samplefile.txt", sep="\t", row.names = 1, header = T)#读取TCGA基因分类信息
key.markers <- read.table(file = "E:/生信学习资料/免疫浸润/Immunoinhibitor.txt",  header=TRUE, row.names=NULL, quote="", sep="\t", check.names = FALSE)
data_cp=na.omit(CRC[key.markers$gene,])


data_cp = as.data.frame(t(data_cp))
data_cp$Sample_Name = rownames(data_cp)
#load("step5.rda")# 加载样本分组信息数据
samplefile$Subtype = factor(samplefile$Subtype, levels = c("cluster1","cluster2"))
samplefile$Sample_Name = rownames(samplefile)
data_cp = merge(samplefile[,c("Sample_Name","Subtype")], data_cp, by="Sample_Name")

#箱线图


library(tidyverse)
library(reshape2)
library(ggsci)
library(ggsignif)
library(ggplot2)
library(ggpubr)
library(tidyr)
long_dat <- pivot_longer(data_cp, cols=3:ncol(data_cp), names_to="Gene", values_to = "Expression")
long_dat$Subtype = factor(long_dat$Subtype)


p=ggboxplot(long_dat, x="Gene", y="Expression", color = "Subtype", 
            ylab="Gene expression",
            xlab="",
            legend.title="Immunoinhibitor",
            palette = c('#58CDD9','#EF767A'),
            width=0.6, add = "none")
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Subtype),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")
outFile="GSE40967_Immunoinhibitor.pdf" 
#输出图片
pdf(file=outFile, width=10, height=5)
print(p1)
dev.off()