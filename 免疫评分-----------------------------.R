#免疫评分-----------------------------
library(utils)
rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
raw_data1 = exp
write.table(raw_data1, file="input_genes.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
filterCommonGenes(input.f="input_genes.txt", output.f="common_genes.gct", id="GeneSymbol")
estimateScore(input.ds="common_genes.gct", output.ds="estimate_score.gct", platform="illumina")
estimate_score <- read.table("estimate_score.gct", skip=2, header=TRUE)
rownames(estimate_score) <- estimate_score[ , 1]
estimate_score <- as.data.frame(t(estimate_score[ , 3:ncol(estimate_score)]))
rownames(estimate_score) = gsub("\\.", "\\-", rownames(estimate_score))
#save(estimate_score, file = "estimate_score.rda")
#分型分组1：
samplefile2=samplefile
samplefile2$Subtype = factor(samplefile2$Subtype, levels = c("cluster1","cluster2"))
#samplefile2$Subtype = factor(samplefile2$Subtype, levels = c("clusterA","clusterB","clusterC")) #,"clusterD"
estimate_score$Sample_Name = rownames(estimate_score)
estimate_score= merge(samplefile2[,c("Sample_Name","Subtype")], estimate_score, by="Sample_Name")

colnames(estimate_score)
library(ggpubr)
ggplot(data =estimate_score,aes(x = Subtype, y = ESTIMATEScore, fill = Subtype))+
  scale_fill_manual(values = c('#EF767A','#58CDD9')) +   #"#FAC795","#82C4E6","#EAA9C1"     "#FAC795","#82C4E6","#EAA9C1","#2BB69E" #"blue","red"
  geom_violin(alpha=0.4, position = position_dodge(width = .75),
              size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7)+ # 背景色透明化
  geom_point(shape = 21, size=2, 
             position = position_jitterdodge(), 
             color="black", alpha=1)+ # 边框线黑色
  theme_classic() +
  ylab(expression("ESTIMATEScore")) +
  xlab("Subtype")  +
  theme(axis.ticks = element_line(size=0.2,color="black"),
        axis.ticks.length = unit(0.2,"cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
   stat_compare_means(comparisons=list(c("cluster1","cluster2")), #,c("clusterA","clusterC"),c("clusterB","clusterC")  #,c("clusterB","clusterD"),c("clusterC","clusterD") 
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols=c("***", "**", "*", "ns")), 
                     label="p.signif")

dev.off()