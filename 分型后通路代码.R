inputFile="combat_exprSet_test.txt"
gmtFile="c2.cp.reactome.v2023.1.Hs.symbols.gmt"                                         

#???ð?
library(GSVA)
library(limma)
library(GSEABase)

#??ȡ?????ļ????????????ļ?????
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
write.table(ssgseaOut,file="ssgseaOut_test_reactom.txt",sep="\t",quote=F,col.names=F)



library(ggplot2)
reactom = fread("reactom_test.txt", header = TRUE, data.table = FALSE)

compaired <- list(c("clusterA", "clusterB"),c("clusterA","clusterC"),c("clusterB","clusterC") )  #设置比较组别

palette<-c(brewer.pal(7,"Set2")[c(1,2,4,5)])   #颜色调用 
colnames(reactom)
ggboxplot(reactom, 
          x = "Subtype", y = "REACTOME_OTHER_INTERLEUKIN_SIGNALING", 
          fill = "Subtype", palette = palette, 
          add = "jitter", size = 0.5)+   #添加抖动的散点 
  stat_compare_means(comparisons = compaired, 
                     method = "wilcox.test",   #设置统计方法 
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns"))) 
dev.off()
