rm(list=ls())
options(stringsAsFactors = F)
setwd("F:/CRC/CRC/TCGA/") 

load("GSE48350.rda")#samplefile文件
CRC <- read.table("change_GSE40967-矩阵.txt", header=TRUE, row.names=1, quote="", sep="\t", check.names = FALSE)
SampleFile <- read.table("WDR72SampleFile.txt", header=TRUE, row.names=NULL, quote="", sep="\t", check.names = FALSE)
#SampleFile=samplefile
exprSet = CRC
exprSet = exprSet[,match(SampleFile$Sample_Name,colnames(exprSet))]#根据samplefile样本顺序排列矩阵顺序
write.table(exprSet, file="WDR72GSVA.txt", sep="\t", quote=FALSE, row.names=T, col.names=TRUE)

#CRC/HC CRC:HC=465:42
load("CRC.rda")
exprSet = exp_set_3[,match(SampleFile$Sample_Name,colnames(exp_set_3))]#根据samplefile样本顺序排列矩阵顺序
write.table(exprSet, file="GSVA507.txt", sep="\t", quote=FALSE, row.names=T, col.names=TRUE)


setwd("F:/CRC/GSVA/46GSVA/")
inputFile="WDR72GSVA.txt"
gmtFile="c5.all.v2023.1.Hs.symbols.gmt"

library(GSVA)
library(limma)
library(GSEABase)
rt <- read.table("WDR72GSVA.txt", header=TRUE, row.names=NULL, quote="", sep="\t", check.names = FALSE)
#rt=read.table(inputFile,sep="\t",header=T,check.names=F,rownames=1)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=normalizeBetweenArrays(mat)

c3gsc2=getGmt( gmtFile, 
               collectionType=BroadCollection(category="c5"), 
               geneIdType=SymbolIdentifier())
gsvaOut=gsva(mat, 
             c3gsc2, 
             min.sz=10, 
             max.sz=500, 
             verbose=TRUE,
             parallel.sz=1)
gsvaOut=rbind(id=colnames(gsvaOut),gsvaOut)
write.table(gsvaOut,file="WDR72GSVAOut.txt",sep="\t",quote=F,col.names=F)

library(limma)

rt=read.table("WDR72GSVAOut.txt",sep="\t",header=T,check.names=F,row.names=1)

logFCcutoff=0#0.2/0.3，实在没有可选0
adjPvalueCutoff=0.05
                                                      #Cluster 1
type=c( rep("low",283),rep("high",283) )    #分组信息及样本量，需修改 CRC:HC=465:42/46两型221：212
design=model.matrix(~ type)
colnames(design)=c("low", "high")
fit=lmFit(rt, design)
fit=eBayes(fit)

all=topTable(fit, coef="low", number=Inf,adjust.method="holm")
all=rbind(id=colnames(all),all)
write.table(all,file="WDR72h_all.txt",sep="\t",quote=F,col.names=F)

diff <- topTable(fit, coef="low", number=Inf,
                 p.value=adjPvalueCutoff, adjust="holm", lfc=logFCcutoff)
diffName=row.names(diff)
diff=rbind(id=colnames(diff),diff)
write.table(diff,file="WDR72h_diff.xls",sep="\t",quote=F,col.names=F)

hmExp=rt[diffName,]
hmExp=rbind(id=colnames(hmExp),hmExp)
write.table(hmExp,file="WDR72h_heatmap.txt",sep="\t",quote=F,col.names=F)

logFoldChange=0 #设置标准和上面一致                                        
adjustP=0.05
                 
rt=read.table("WDR72h_all.txt",sep="\t",header=T,check.names=F)                  

tiff(file="vol.tiff",
     width = 13,           
     height =13,            #ͼƬ?ĸ߶?
     units ="cm",
     compression="lzw",
     bg="white",
     res=300)  #分辨率300
yMax=10      #根据all.txt中adj.P.Val的最大值调整
xMax=1    #根据all.txt中logFC的最大值调整，在范围内即可
plot(rt$logFC,-log10(rt$adj.P.Val), ylab="-log10(adj.P.Val)",xlab="logFC",
     main="Volcanoh", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=1)  #cex为点的大小
diffSub=subset(rt, adj.P.Val<adjustP & logFC>logFoldChange)
points(diffSub$logFC,-log10(diffSub$adj.P.Val), pch=20, col="red",cex=2)
diffSub=subset(rt, adj.P.Val<adjustP & logFC<(-logFoldChange))
points(diffSub$logFC, -log10(diffSub$adj.P.Val), pch=20, col="green",cex=2)
abline(v=0,lty=2,lwd=3)
dev.off()

rt=read.table("WDR72h_heatmap.txt",sep="\t",header=T,row.names=1,check.names=F)

library(pheatmap)
Type=c(rep("low",283),rep("high",283))    
names(Type)=colnames(rt)
Type=as.data.frame(Type)

tiff(file="hheatmap0.1.tiff",
     width = 20,            #样本多调整宽度
     height =20,            #通路多调整高度
     units ="cm",
     compression="lzw",
     bg="white",
     res=600)
pheatmap(rt, 
         annotation=Type, 
         color = colorRampPalette(c("#2f5688", "#BBBBBB", "#CC0000"))(50),
         cluster_cols =F,
         fontsize = 6,
         fontsize_row=3,   #字体大小
         fontsize_col=5)
dev.off()

#聚类热图------------------------------------------------
load("nmf.rda")
load("STAD.rda")
#SampleFile <- read.table("ITIH3samplefile.txt", header=TRUE, row.names=NULL, quote="", sep="\t", check.names = FALSE)
KEGG <- read.table("ITIH3_goGSVAOut.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
#KEGG = KEGG[,-1]
KEGG1 = as.data.frame(t(KEGG))


exprSet=as.data.frame(STAD)
group <- read.table("ITIH3samplefile.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
#自定义分组颜色；
group_colors=list(type=c(clusterA = "#EF767A", clusterB = "#456990", clusterC = "#48C0AA"))
#添加分组颜色条；
library(pheatmap)
mycol1<-colorRampPalette(c("navy", "white", "orange"))(100)
exprSet = exprSet[apply(exprSet, 1, function(x) sd(x)!=0),] ## 删掉标准差为0的行
p <- pheatmap(exprSet,scale="row",
              fontsize = 6,
              fontsize_row = 6,
              fontsize_col = 7,
              cellwidth = 25,
              cellheight = 6,
              color = mycol1,
              border="white",
              treeheight_col=30,
              cluster_cols = F,
              gaps_col=c(3,6),
              cutree_rows=3,
              annotation_col=group,
              annotation_color=group_colors)

#提取热图的行方向（基因）的聚类树；
clu <- p$tree_row
#对聚类树进行分簇；
cluster <- factor(cutree(clu,4))
cluster
#转成数据框；
cut.df <- data.frame(cluster)
#绘制聚类树；
plot(clu,hang = -1,cex=0.6,axes=FALSE,ann=FALSE)

#绘制热图；
pheatmap(exprSet,scale="row",
         fontsize = 6,
         fontsize_row = 6,
         fontsize_col = 7,
         cellwidth = 1,
         cellheight = 6,
         color = mycol1,
         border="white",
         treeheight_col=30,
         cluster_cols = F,
         cutree_cols=3,
         cutree_rows=4,
         annotation_col=group,
         annotation_row =cut.df,
         annotation_color=group_colors)
dev.off()
#通路差异-----------------------------------------
rm(list=ls())
options(stringsAsFactors = F)
setwd("F:/CRC/CRC/TCGA/")
load("46step2_km433.rda")
load("TCGA_CRC.rda")
exprSet = CRC_ColsAllNa_removed[,match(samplefile$Sample_Name,colnames(CRC_ColsAllNa_removed))]#根据samplefile样本顺序排列矩阵顺序
write.table(exprSet, file="433exprSet.txt", sep="\t", quote=FALSE, row.names=T, col.names=TRUE)

inputFile="433exprSet.txt"
gmtFile="c2.cp.kegg.v2023.1.Hs.symbols.gmt"                                         

#???ð?
library(GSVA)
library(limma)
library(GSEABase)

rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
geneSet=getGmt(gmtFile, 
               geneIdType=SymbolIdentifier())

#ssgsea????
ssgseaScore=gsva(rt, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#????ssGSEA score????????
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
#??ssGSEA score???н???
ssgseaOut=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
write.table(ssgseaOut,file="433KEGGgseaOut_test.txt",sep="\t",quote=F,col.names=F)

library(ggpubr)
library(data.table)
library(ggplot2)
#reactom = fread("433KEGGgseaOut_test.txt", header = TRUE, data.table = FALSE)
reactom <- read.table("433KEGGgseaOut_test.txt", header=TRUE, row.names=1, quote="", sep="\t", check.names = F)

#compaired <- list(c("clusterA", "clusterB"),c("clusterA","clusterC"),c("clusterB","clusterC") )  #设置比较组别
compaired <- list("cluster1", "cluster2" )
palette<-c(brewer.pal(7,"Set2")[c(1,2,4,5)])   #颜色调用 
colnames(reactom)
ggboxplot(reactom, 
          x = "Subtype", y = "KEGG_N_GLYCAN_BIOSYNTHESIS", 
          fill = "Subtype", palette = palette, 
          add = "jitter", size = 0.5)+   #添加抖动的散点 
  stat_compare_means(comparisons = compaired, 
                     method = "wilcox.test",   #设置统计方法 
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns"))) 
dev.off()

