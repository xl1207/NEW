setwd("F:/GSE/PP/GSE13355.1.5/")
#combat------------------------------------------------------------
rm(list = ls())
options(stringsAsFactors = F)
load('GSE13355.rda')
exprSet1 <- exprSet
samplefile1 <- SampleFile
load('GSE14905.rda')
exprSet2 <- exprSet
samplefile2 <- SampleFile
load('GSE78097.rda')
exprSet3 <- exprSet
samplefile3 <- SampleFile
load('GSE109248.rda')
exprSet4 <- exprSet
samplefile4 <- SampleFile
load('GSE53431.rda')
exprSet5 <- exprSet
samplefile5 <- SampleFile

common_gene <- Reduce(intersect, list(rownames(exprSet1), rownames(exprSet2), 
                                      rownames(exprSet3), rownames(exprSet4), 
                                      rownames(exprSet5)
                                      ))

combined_exprSet <- cbind.data.frame(exprSet1[common_gene, ], exprSet2[common_gene, ], 
                                     exprSet3[common_gene, ], exprSet4[common_gene, ], 
                                     exprSet5[common_gene, ]
                                     )

combined_samplefile <- rbind.data.frame(cbind(samplefile1, Batch=rep("GSE13355_GPL570",nrow(samplefile1))),
                                        cbind(samplefile2, Batch=rep("GSE14905_GPL570",nrow(samplefile2))),
                                        cbind(samplefile3, Batch=rep("GSE78097_GPL570",nrow(samplefile3))),
                                        cbind(samplefile4, Batch=rep("GSE109248_GPL10558",nrow(samplefile4))),
                                        cbind(samplefile5, Batch=rep("GSE53431_GPL10558",nrow(samplefile5)))
                                        )

colnames(combined_samplefile)[2] = "Sample_Group"
#combined_samplefile$Sample_Group = factor(combined_samplefile$Sample_Group, levels = c("NL" ,"LS"))

combat_samplefile = combined_samplefile
combat_exprSet = combined_exprSet
combat_exprSet = combat_exprSet[, match(combat_samplefile$Sample_Name, colnames(combat_exprSet))]
library(FactoMineR)
library(factoextra)
dist_mat <- dist(t(combat_exprSet))
combat_samplefile$Batch <- as.factor(combat_samplefile$Batch)
pre.pca <- PCA(t(combat_exprSet),graph = FALSE)
fviz_pca_ind(pre.pca,
             geom= "point",
             col.ind = combat_samplefile$Batch,
             addEllipses = TRUE,
             legend.title="Group")

#在校正的时候要保留生物学上的差异，不能矫正过枉。
model <- model.matrix(~as.factor(combat_samplefile$Sample_Group))
library(sva)
combat_exprSet <- ComBat(dat = combat_exprSet,batch = combat_samplefile$Batch,mod = model)
combat.pca <- PCA(t(combat_exprSet),graph = FALSE)
fviz_pca_ind(combat.pca,
             geom= "point",
             col.ind = combat_samplefile$Batch,
             addEllipses = TRUE,
             legend.title="Group"  )

combat_samplefile = combat_samplefile[,c(1,2)]
combat_samplefile = combat_samplefile[order(combat_samplefile$Sample_Group, decreasing = T), ]

save(combat_exprSet, combat_samplefile, file = "step1_combat.rda")
write.table(combat_samplefile, file = "combat_samplefile.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
