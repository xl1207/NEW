rm(list = ls())
options(stringsAsFactors = F)
setwd("F:/TB_T2D/T2D/GSE38642/")
load("GSE38642.rda")
new_exprSet = new_exprSet[,-1]
new_exprSet = exprSet
new_exprSet = as.data.frame(normalizeBetweenArrays(as.matrix(new_exprSet), method = "quantile")) 
new_samplefile = SampleFile
new_samplefile$Sample_Group = factor(new_samplefile$Sample_Group, levels = c("T2D", "HC"))
design = model.matrix(~ 0 + new_samplefile$Sample_Group)
colnames(design) = levels(new_samplefile$Sample_Group)
contrast.matrix = makeContrasts("T2D-HC", levels = design)
#contrast.matrix = makeContrasts("SLE-Normal", "SLE_LDG-Normal", "SLE_LDG-SLE", levels = design)

deg = function(new_exprSet, design, contrast.matrix, coef){
  fit = lmFit(new_exprSet, design)
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

nrDEG1_all <- nrDEG1[which(nrDEG1$adj.P.Val < 0.05 & abs(nrDEG1$logFC) >= 0.32), ]
nrDEG1_all[which(nrDEG1_all$logFC > 0), "group"] <- "up"
nrDEG1_all[which(nrDEG1_all$logFC < 0), "group"] <- "down"
table(nrDEG1_all$group)
write.table(nrDEG1_all, file = "GSE19444_1.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
save(new_exprSet, new_samplefile, nrDEG1, nrDEG1_all, file = "GSE19444_1.rda")
