rm(list=ls())
setwd("F:/STAD/TCGA/")

#输入表达矩阵TCGA---------------------------------------
load("STAD.rda")
samplefile2=read.table("STAD_SampleFile.txt", header=T, sep="\t", check.names=F,row.names =NULL)
rownames(samplefile2)=samplefile2$Sample_Name
exp = as.data.frame(t(STAD))
#samplefile2 = SampleFile
#samplefile2 = samplefile2[samplefile2$Sample_Group =="CRC",]
#exp= as.data.frame(exp[samplefile2$Sample_Name, ]) #colnames(exp)
#CMS4表达矩阵
#samplefile <- read.table("TCGACMS4_SampleFile.txt", header=TRUE, row.names=NULL, quote="", sep="\t", check.names = FALSE)
exp=exp[samplefile2$Sample_Name,]

#exp=as.data.frame(t(exp))
cut.off = median(exp$ITIH3)   ##NOTCH3，RUNX2还有FJX1,MCTP1

#正常分组
samplefile2$Sample_Group = ifelse(exp$ITIH3 >= cut.off, 'high', 'low')
samplefile2$Sample_Group = factor(samplefile2$Sample_Group, levels = c("high", "low"))
table(samplefile2$Sample_Group)
exp = as.data.frame(t(exp))

write.table(exp, file="exp.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

inputFile="exp.txt"
gmtFile="F:/CRC/GSVA/c2.cp.kegg.v2023.1.Hs.symbols.gmt"                                         

#???ð?
library(GSVA)
library(limma)
library(GSEABase)


rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
geneSet=getGmt(gmtFile, 
               geneIdType=SymbolIdentifier())


ssgseaScore=gsva(rt, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
#??ssGSEA score???н???
ssgseaOut=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
write.table(ssgseaOut,file="ssgseaOut_test_kegg.txt",sep="\t",quote=F,col.names=F)

library(ggpubr)
library("data.table")
library("ggplot2")
library(ggplot2)
library(classInt) 
library(RColorBrewer)
reactom <- read.table("ssgseaOut_test_kegg.txt", header=TRUE, row.names=1, quote="", sep="\t", check.names = FALSE)
reactom = as.data.frame(t(reactom))



names <- rownames(reactom)   #把行名变为第一列
rownames(reactom) <- NULL
reactom <- cbind(names,reactom)
colnames(reactom)[1]="Sample_Name"
reactom = merge(samplefile2[,c("Sample_Name","Sample_Group")], reactom, by="Sample_Name")
compaired <- list(c("high", "low") )  #设置比较组别
reactom$Sample_Group <- factor(reactom$Sample_Group, levels = c("low", "high"))
#palette<-c(brewer.pal(7,"Set2")[c(1,2,4,5)])   #颜色调用 
palette<-c("#008ECB","#D14039", "#EA921D" ) 
colnames(reactom)
ggboxplot(reactom, 
          x = "Sample_Group", y = "KEGG_BUTANOATE_METABOLISM", 
          fill = "Sample_Group", palette = palette, 
          add = "jitter", size = 0.5)+   #添加抖动的散点 
  stat_compare_means(comparisons = compaired, 
                     method = "wilcox.test",   #设置统计方法 
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns"))) 
dev.off()