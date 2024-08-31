#-----------------------------------------------------------------------------------------
rm(list = ls())
options(stringsAsFactors = F)
setwd("GSE42743/")
##1.使用GPL文件进行探针注释-------------------------------------------------------------------------------------------
exprSet <- read.table("C:/R/GSE/GSE41613/GSE41613_series_matrix.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
rownames(exprSet) = paste0("tmp", rownames(exprSet))
ids <- read.table("c:/R/GPL/GPL570.txt", header = TRUE, sep = "\t", row.names = NULL, check.names = FALSE)
ids$probe_id = paste0("tmp", ids$probe_id)
if(T){
  colnames(ids) = c("probe_id", "symbol")
  exprSet = exprSet[rownames(exprSet) %in% ids$probe_id, ]
  ids = ids[match(rownames(exprSet), ids$probe_id), ]
  
}
##2.使用最大值进行探针合并------------------------------------------------------------------------------------------
if(T){
  probeid2symbol <- function(expeSet, ids){
    tmp = by(exprSet, ids[ , 2], function(x) rownames(x)[which.max(rowMeans(x))]) 
    probes = as.character(tmp)
    print(dim(exprSet))
    
    new_exprSet = exprSet[rownames(exprSet) %in% probes, ]
    print(dim(new_exprSet))
    
    return(new_exprSet)
  }
  new_exprSet = probeid2symbol(expeSet, ids)
  rownames(new_exprSet) = ids[match(rownames(new_exprSet), ids$probe_id), 2]
  
}
new_exprSet <- as.data.frame(new_exprSet)
new_exprSet <- cbind(rownames(new_exprSet), new_exprSet)
colnames(new_exprSet)[1] <- "gene_symbol"
new_exprSet = new_exprSet[order(new_exprSet$gene_symbol), ]
write.table(new_exprSet, file = "GSE41613_gene_matrix.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
#分别导入matrix--------------------------------------------------------------------------------------------
rm(list = ls())
options(stringsAsFactors = F)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE) 
setwd("C:/R/GSE/")
#-------------------------------------------------------------------
library(data.table)
library(tidyverse)
library(pheatmap)
library(sva)
exprSet1 = fread("GSE42743log_gene_matrix.txt", header = TRUE, data.table = FALSE) %>% column_to_rownames("gene_symbol")
samplefile1 = fread("GSE42743Samplefile.txt", header = TRUE, data.table = FALSE)
table(samplefile1$Sample_Group)
exprSet1 = exprSet1[, match(samplefile1$Sample_Name, colnames(exprSet1))]

c <- cor(exprSet1, method = "spearman")
pheatmap(c)
tmp1 <- exprSet1[, rowMeans(c) > 0.7]
c1 <- cor(tmp1, method = "spearman")
pheatmap(c1)
setdiff(colnames(exprSet1), colnames(tmp1))
exprSet1 = exprSet1[,colnames(tmp1)]
samplefile1 = samplefile1[samplefile1$Sample_Name %in% colnames(tmp1),]
save(exprSet1, samplefile1, file="GSE65858.rda")
#-------------------------------------------------------------------
exprSet2 = fread("GSE42743log_gene_matrix.txt", header = TRUE, data.table = FALSE) %>% column_to_rownames("gene_symbol")
samplefile2 = fread("GSE42743Samplefile.txt", header = TRUE, data.table = FALSE)
table(samplefile2$Sample_Group)
exprSet2 = exprSet2[, match(samplefile2$Sample_Name, colnames(exprSet2))]

c <- cor(exprSet2, method = "spearman")
pheatmap(c)
tmp1 <- exprSet2[, rowMeans(c) > 0.7]
c1 <- cor(tmp1, method = "spearman")
pheatmap(c1)
setdiff(colnames(exprSet2), colnames(tmp1))
exprSet2 = exprSet2[,colnames(tmp1)]
samplefile2 = samplefile2[samplefile2$Sample_Name %in% colnames(tmp1),]
save(exprSet2, samplefile2, file="GSE65858.rda")
#-------------------------------------------------------------------
exprSet3 = fread("GSE42743log_gene_matrix.txt", header = TRUE, data.table = FALSE) %>% column_to_rownames("gene_symbol")
samplefile3 = fread("GSE42743Samplefile.txt", header = TRUE, data.table = FALSE)
table(samplefile3$Sample_Group)
exprSet3 = exprSet3[, match(samplefile3$Sample_Name, colnames(exprSet3))]

c <- cor(exprSet3, method = "spearman")
pheatmap(c)
tmp1 <- exprSet3[, rowMeans(c) > 0.7]
c1 <- cor(tmp1, method = "spearman")
pheatmap(c1)
setdiff(colnames(exprSet3), colnames(tmp1))
exprSet3 = exprSet3[,colnames(tmp1)]
samplefile3 = samplefile3[samplefile3$Sample_Name %in% colnames(tmp1),]
save(exprSet3, samplefile3, file="GSE65858.rda")
#-------------------------------------------------------------------
exprSet4 = fread("./GSE47404/GSE47404_gene_matrix.txt", header = TRUE, data.table = FALSE) %>% column_to_rownames("gene_symbol")
samplefile4 = fread("./GSE47404/SampleFile.txt", header = TRUE, data.table = FALSE)
table(samplefile4$Sample_Group)
exprSet4 = exprSet4[, match(samplefile4$Sample_Name, colnames(exprSet4))]

c <- cor(exprSet4, method = "spearman")
pheatmap(c)
tmp1 <- exprSet4[, rowMeans(c) > 0.7]
c1 <- cor(tmp1, method = "spearman")
pheatmap(c1)
save(exprSet4, samplefile4, file="4_GSE47404.rda")
#-------------------------------------------------------------------
exprSet5 = fread("./GSE53625/GSE53625_gene_matrix.txt", header = TRUE, data.table = FALSE) %>% column_to_rownames("gene_symbol")
samplefile5 = fread("./GSE53625/SampleFile.txt", header = TRUE, data.table = FALSE)
table(samplefile5$Sample_Group)
samplefile5 = samplefile5[samplefile5$Sample_Group != "normal",]
exprSet5 = exprSet5[, match(samplefile5$Sample_Name, colnames(exprSet5))]

write.table(samplefile5, file = "./GSE53625/SampleFile.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

c <- cor(exprSet5, method = "spearman")
pheatmap(c)
tmp1 <- exprSet5[, rowMeans(c) > 0.70]
c1 <- cor(tmp1, method = "spearman")
pheatmap(c1)
save(exprSet5, samplefile5, file="5_GSE53625.rda")
#-------------------------------------------------------------------
exprSet6 = fread("./GSE44021/GSE44021_gene_matrix.txt", header = TRUE, data.table = FALSE) %>% column_to_rownames("gene_symbol")
samplefile6 = fread("./GSE44021/SampleFile.txt", header = TRUE, data.table = FALSE)
table(samplefile6$Sample_Group)
samplefile6 = samplefile6[samplefile6$Sample_Group != "normal",]
exprSet6 = exprSet6[, match(samplefile6$Sample_Name, colnames(exprSet6))]

c <- cor(exprSet6, method = "spearman")
pheatmap(c)
tmp1 <- exprSet6[, rowMeans(c) > 0.80]
c1 <- cor(tmp1, method = "spearman")
pheatmap(c1)
save(exprSet6, samplefile6, file="6_GSE44021.rda")

load(GSE42743.rda)
#去除批间差---------------------------------------------------------
common_gene <- Reduce(intersect, list(rownames(exprSet1), rownames(exprSet2), 
                                      rownames(exprSet3), #rownames(exprSet4),
                                      rownames(exprSet5)))#, rownames(exprSet6)))

combined_exprSet <- cbind.data.frame(exprSet1[common_gene, ], exprSet2[common_gene, ],
                                     exprSet3[common_gene, ], #exprSet4[common_gene, ], 
                                     exprSet5[common_gene, ])#, exprSet6[common_gene, ])

combined_samplefile <- rbind.data.frame(cbind(samplefile1, Batch=rep("GSE69925_GPL570",nrow(samplefile1))),
                                        cbind(samplefile2, Batch=rep("GSE121931_GPL15207",nrow(samplefile2))),
                                        cbind(samplefile3, Batch=rep("GSE66258_GPL10558",nrow(samplefile3))),
                                        #cbind(samplefile4, Batch=rep("GSE47404_GPL6480",nrow(samplefile4))),
                                        cbind(samplefile5))
#,cbind(samplefile6, Batch=rep("GSE44021_GPL571",nrow(samplefile6))))
#combined_samplefile$Sample_Group = factor(combined_samplefile$Sample_Group, levels = c("HC", "SLE"))
#
source("./batchPCA.R")
batchPCA(indata = t(scale(t(combined_exprSet))), batch = combined_samplefile$Batch, fig.dir = ".",
         PCA.fig.title = "Raw PCA for combined expression profile",
         #cols = c(blue, yellow, green),
         showID = FALSE, cex = 0.7, showLegend = TRUE, pos="topright")
library(sva)
#mod = model.matrix(~as.factor(combined_samplefile$Sample_Group), data=combined_samplefile)
combat_exprSet <-  as.data.frame(ComBat(dat = as.matrix(combined_exprSet), batch = combined_samplefile$Batch))#, mod = mod

batchPCA(indata = t(scale(t(combat_exprSet))), batch = combined_samplefile$Batch, fig.dir = ".",
         PCA.fig.title = "Combat PCA for combined expression profile",
         #cols = c(blue, yellow, green),
         showID = FALSE, cex = 0.7, showLegend = TRUE, pos="topright")
#combat_samplefile = combined_samplefile[,c(1,3)]
#combat_samplefile = combat_samplefile[order(combat_samplefile$Sample_Group, decreasing = T), ]
combat_samplefile = combined_samplefile
combat_exprSet = combat_exprSet[, match(combat_samplefile$Sample_Name, colnames(combat_exprSet))]
save(combat_exprSet, combat_samplefile, file = "step1.rda")

#样本相关性
c <- cor(combat_exprSet, method = "spearman")
pheatmap(c)
tmp1 <- combat_exprSet[, rowMeans(c) > 0.70]
c1 <- cor(tmp1, method = "spearman")
pheatmap(c1)
combat_exprSet = combat_exprSet[,colnames(tmp1)]
combat_samplefile = combat_samplefile[combat_samplefile$Sample_Name %in% colnames(tmp1),]
save(combat_exprSet, combat_samplefile, file = "step2.rda")
#差异表达-----------------------------------------------------------
rm(list = ls())
options(stringsAsFactors = F)
setwd("F:/TB_T2D/TB/GSE19444/")
load("GSE19444.rda")
#new_exprSet = as.data.frame(normalizeBetweenArrays(as.matrix(new_exprSet), method = "quantile")) 
#new_samplefile = combat_samplefile
SampleFile$Sample_Group = factor(SampleFile$Sample_Group, levels = c("HC", "TB"))
library(limma)
design = model.matrix(~ 0 + SampleFile$Sample_Group)
colnames(design) = levels(SampleFile$Sample_Group)
contrast.matrix = makeContrasts("TB-HC", levels = design)
#contrast.matrix = makeContrasts("SLE-Normal", "SLE_LDG-Normal", "SLE_LDG-SLE", levels = design)

deg = function(exprSet, design, contrast.matrix, coef){
  fit = lmFit(exprSet, design)
  library(futile.logger)
  fit2 = contrasts.fit(fit, contrast.matrix)
  fit2 = eBayes(fit2)
  tepmOutput = topTable(fit2, coef = coef, adjust.method = "fdr", number = 50000)
  nrDEG = na.omit(tepmOutput)
  return(nrDEG)
  
}
nrDEG1 = deg(new_exprSet, design, contrast.matrix, coef = 1)
nrDEG1 <- as.data.frame(nrDEG1)
nrDEG1 <- cbind(rownames(nrDEG1), nrDEG1)
colnames(nrDEG1)[1] <- "gene_symbol"

nrDEG1_all <- nrDEG1[which(nrDEG1$adj.P.Val < 0.05 & abs(nrDEG1$logFC) >= 0.58), ]
nrDEG1_all[which(nrDEG1_all$logFC > 0), "group"] <- "up"
nrDEG1_all[which(nrDEG1_all$logFC < 0), "group"] <- "down"
table(nrDEG1_all$group)
#write.table(nrDEG1_all, file = "nrDEG_all.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
save(new_exprSet, new_samplefile, nrDEG1, nrDEG1_all, file = "step2.rda")

#可视化-----------------------------
library(ggplot2)
library(ggrepel)
data = nrDEG1
data$Group <- as.factor(ifelse(data$adj.P.Val < 0.01 & abs(data$logFC) >= 0.58, ifelse(data$logFC >= 0.58, 'up-regulated','down-regulated'),'not-significant'))
p <- ggplot(data = data, aes(x = logFC, y = -log10(adj.P.Val), color = Group))+ 
  geom_point(alpha = 0.8, size = 1.2)+  
  scale_color_manual(values = c("blue", "black", "red")) +  
  labs(title="Volcano plot", x="log2FoldChange", y="-log10(Agjust P-value)") +
  theme(plot.title = element_text(hjust = 0.4)) + 
  geom_hline(yintercept = -log10(0.01), lty = 4, lwd = 0.6, alpha = 0.8) + geom_vline(xintercept = c(0.58,-0.58), lty = 4, lwd = 0.6,alpha = 0.8) + 
  theme_bw() + 
  geom_text_repel(data = subset(data, abs(logFC) >= 1.5 & adj.P.Val < 1e-20), aes(label = gene_symbol),size = 3, segment.color = "black", show.legend = FALSE) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p
ggsave(p, filename = "volcano1.pdf", width = 20, height = 15, units = c("cm"))
dev.off()