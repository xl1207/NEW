rm(list = ls())
setwd("F:/GSE/GBM/")
#du <- read.table("GSE13355_series_matrix.txt", header = TRUE, sep = "\t", row.names = NULL, check.names = FALSE)
#1.找一个数据集做训练集-------------------
load('GSE13355.rda')
library(limma)
new_exprSet = as.data.frame(normalizeBetweenArrays(as.matrix(exprSet), method = "quantile")) 
SampleFile$Sample_Group = factor(SampleFile$Sample_Group, levels = c("tumor", "control"))
design = model.matrix(~ 0 + SampleFile$Sample_Group)
colnames(design) = levels(SampleFile$Sample_Group)
contrast.matrix = makeContrasts("tumor-control", levels = design)
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
#logFC卡差异基因
nrDEG1_all <- nrDEG1[which(nrDEG1$adj.P.Val < 0.05 & abs(nrDEG1$logFC) >= 1), ]#logFC至少为1
nrDEG1_all[which(nrDEG1_all$logFC > 0), "group"] <- "up"
nrDEG1_all[which(nrDEG1_all$logFC < 0), "group"] <- "down"
table(nrDEG1_all$group)#显示差异基因
write.table(nrDEG1_all, file = "差异基因1.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
save(exprSet,SampleFile, nrDEG1, nrDEG1_all, file = "step1_DEG_1.rda")

#可视化-火山图和热图-----------------------------
load("step1_DEG_1.rda")
#volcano2
library(ggpubr)
library(ggthemes)
data = nrDEG1
data$Group <- as.factor(ifelse(data$adj.P.Val < 0.05 & abs(data$logFC) >= 1, ifelse(data$logFC >= 1, 'up-regulated','down-regulated'),'not-significant'))
data$Label <- ""
data <- data[order(data$adj.P.Val), ]
up.genes <- head(data$gene_symbol[which(data$Group == "up-regulated")], 10)
down.genes <- head(data$gene_symbol[which(data$Group == "down-regulated")], 10)
top10.genes <- c(as.character(up.genes), as.character(down.genes))
data$Label[match(top10.genes, data$gene_symbol)] <- top10.genes
data$logP <- -log10(data$adj.P.Val)
p <-ggscatter(data, x = "logFC", y = "logP", color = "Group", 
              palette = c("#2f5688", "#BBBBBB", "#CC0000"), size = 1, label = data$Label, 
              font.label = 8, repel = TRUE, xlab = "log2FoldChange", ylab = "-log10(Agjust P-value)") + 
  theme_base() + 
  geom_hline(yintercept = 1.30, linetype = "dashed") + geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed")
p
ggsave(p, filename = "volcano2xing.pdf", width = 20, height = 15, units = c("cm"))
dev.off()
#heatmap
#heapmap1
library(pheatmap)                               #两型已修改
data <- as.matrix(exp)  #
#data <- as.matrix(exprSet)
x = nrDEG1_all$logFC
names(x) = nrDEG1_all$gene_symbol
choose_gene = c(names(head(sort(x, decreasing = TRUE), 30)), names(tail(sort(x, decreasing = TRUE), 30)))
choose_matrix = data[choose_gene, ]
#choose_matrix = exprSet[names(sort(apply(exprSet, 1, mad), decreasing = T)[1:500]), ]
choose_matrix = t(scale(t(choose_matrix)))
choose_matrix[choose_matrix > 2] = 2
choose_matrix[choose_matrix< -2] = -2
annotation_col = data.frame(Group = samplefile$Subtype)  #sample_name
#annotation_col = data.frame(Group = SampleFile$Sample_Group)
rownames(annotation_col) = colnames(data)
color <- colorRampPalette(c("#3300CC", "#3399FF", "white", "#FF3333", "#CC0000"))(50)
pdf("Heatmap2xing.pdf", width = 8, height = 8)
#pdf("Heatmap2.pdf", width = 8, height = 8)
pheatmap(choose_matrix, color = color, border_color = NA, annotation_col = annotation_col, 
         cluster_row = F, cluster_cols = F, show_colnames = F)
dev.off()

#heatmap2
#library(pheatmap)
#bk = unique(c(seq(-2, 2, length = 100)))
#color <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
#pdf("Heatmap.pdf", width = 10, height = 25)
#pheatmap(choose_matrix, annotation_col = annotation_col, labels_col = '', breaks = bk, 
#         color = color, cluster_rows = F, cluster_cols = F, show_rownames = TRUE, show_colnames = FALSE)
#dev.off()


#2.WGCNA分析 PP-------------------------------
rm(list = ls())
library(WGCNA)
library(RColorBrewer)
options(stringsAsFactor = FALSE)
allowWGCNAThreads()

#Step1
##表达谱数据准备
#expr <- read.table("GSE58095_gene_matrix.txt", header=TRUE, row.names=1, quote="", sep="\t", check.names = FALSE)
#SampleFile = read.table("SamplefileGSE56815.txt", header=TRUE, row.names=NULL, quote="", sep="\t", check.names = FALSE)
setwd("F:/GSE/PP/GSE13355U/")
load('GSE13355.rda')
expr = exprSet
dim(expr)
m.vars = apply(expr, 1, var)
expr1 = expr[which(m.vars > quantile(m.vars, probs = seq(0, 1, 0.25))[4]), ]#25%的基因做
#expr1 = expr[order(apply(expr,1,mad), decreasing = T)[1:5000],] #前5000个基因做

dim(expr1)

datExpr0 = as.data.frame(t(expr1))

##样本评估矩阵信息是否合格
gsg = goodSamplesGenes(datExpr0, verbose=3)
gsg$allOK

##对样本进行聚类分析(绘制样本聚类树，剔除离群样本)
sampleTree = hclust(dist(datExpr0), method="average")

#pdf(file="01.sample_cluster_tree.pdf", onefile=FALSE, paper="special", width=25, height=10, bg="white", pointsize=6)
par(mar = c(0, 4, 2, 0))
par(cex = 0.6)
plot(sampleTree, main="Sample clustering to detect outliers", sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)
cutHeight = 150 #若GSM离群，需要修改,介于一至二节点即可
abline(h=cutHeight, col="red")
dev.off()

clust = cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 10)
table(clust)

keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#表型数据(ClinicalTraits.txt), 重建样本聚类树(添加表型展示图)
#datTraits = read.table("~/Transcriptome/Data-SLE/trimmed_output/SampleFile.txt", sep='\t', header=TRUE, check.names=FALSE)
#traitRows = match(rownames(datExpr), datTraits$Sample_Name)
#datTraits = datTraits[traitRows, ]
#rownames(datTraits) = datTraits$Sample_Name
#datTraits$Subtype <- factor(datTraits$Subtype, levels = c("control", "case1", "case2"))
combat_samplefile1 = SampleFile

datTraits = combat_samplefile1
rownames(datTraits) = datTraits$Sample_Name
datTraits$Sample_Group = factor(datTraits$Sample_Group)
#datTraits = datTraits[,c(23, 9)]
#colnames(datTraits)[2] = "Tregs"
#datTraits$`Comment_anti-Ro` <- factor(datTraits$`Comment_anti-Ro`, levels = c("control", "none", "med", "high"))

#sampleTree2 = hclust(dist(datExpr), method="average")
#datTraitColors = numbers2colors(datTraits, signed = FALSE)#colorWhiteRed(50)
#pdf(file="01.sample_cluster_trait.pdf", onefile=FALSE, paper="special", width=25, height=10, bg="white", pointsize=6)
#plotDendroAndColors(sampleTree2, datTraitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")
#可选择参数:colorHeight=0.2, colorHeightBase=0.2, colorHeightMax=0.4, rowWidths=NULL, dendroLabels=NULL, addGuide=FALSE, guideAll=FALSE, guideCount=50, guideHang=0.2, addTextGuide=FALSE, cex.colorLables=0.8, cex.dendroLables=0.7, cex.rowText=0.8, marAll=c(1,5,3,1), saveMarTRUE)
#dev.off()

#Step2
##选择合适的软域值
powers = c(c(1:10), seq(from=12, to=20, by=2))#下一句 软阈值表
sft = pickSoftThreshold(datExpr, powerVector=powers, verbose=5, networkType="unsigned")

#pdf(file="02.sample_soft_threshold.pdf", onefile=FALSE, paper="special", width=25, height=10, bg="white", pointsize=6)
par(mfrow=c(1,2))
cex1 = 0.6
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2], xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n", main=paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2], labels=powers, cex=cex1, col="red")
abline(h=0.85, col="blue")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main=paste("Mean Connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels=powers, cex=cex1, col="red")
dev.off()

softPower =sft$powerEstimate#12
#softPower =20
#-------------------------------------------------------------------
#
#ADJ1 = abs(cor(datExpr, use="p"))^18#需要修改
#k = as.vector(apply(ADJ1, 2, sum, na.rm=TRUE))

#sizeGrWindow(25, 10)
#par(mfrow=c(1, 2))
#hist(k)
#scaleFreePlot(k, main="Check scale free topology\n")

#软域值的检验
#Connectivity = softConnectivity(datExpr, power=softPower)
#par(mfrow=c(1, 1))
#scaleFreePlot(Connectivity, main=paste("soft threshold, power=", softPower), truncated=TRUE)
#-------------------------------------------------------------------
#Step3
#邻接矩阵和拓扑矩阵，节点聚类树
adjacency = adjacency(datExpr, power=softPower, type="unsigned")
TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM

geneTree = hclust(as.dist(dissTOM), method="average")
#geneTree = flashClust(as.dist(dissTOM), method="average")

#Dynamic Tree Cut 1
#pdf(file="03.man_cluster_tree.pdf", onefile=FALSE, paper="special", width=25, height=10, bg="white", pointsize=6)
plot(geneTree, xlab="", sub="", main="Gene clustering on TOM-based dissimilarity", label=FALSE, hang=0.04)
dev.off()

minModuleSize = 25#每个模块至少n个基因，25/30/50/80
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

#Dynamic Tree Cut 2
#pdf(file="03.man_modules_color.pdf", onefile=FALSE, paper="special", width=25, height=10, bg="white", pointsize=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

##合并相似模块
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1 - cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average");

plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

MEDissThres = 0.25#该值以下模块合并
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

#Merged dynamic
#pdf(file="03.comparison_modules_color.pdf", onefile=FALSE, paper="special", width=25, height=10, bg="white", pointsize=6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleColors = mergedColors
table(moduleColors)
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder) - 1
table(moduleLabels)
MEs = mergedMEs

#Step4(无表型数据时可做部分)---------------------------------------
##ME计算
samplefile_design <- na.omit(datTraits[match(datTraits$Sample_Name,rownames(datExpr)),])
design = model.matrix(~0 + samplefile_design$Sample_Group)
colnames(design)=levels(samplefile_design$Sample_Group)


nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

#表型+模块相关性热图
moduleTraitCor = cor(MEs, design, use="p")#datTraits$Tregs
#moduleTraitCor = cor(MEs, datTraits$Tregs, use="p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

#pdf(file="04.modules_traits_relationships.pdf", onefile=FALSE, paper="special", width=25, height=10, bg="white", pointsize=6)
par(mar = c(6, 9, 3, 3))
#labeledHeatmap(Matrix = moduleTraitCor, xLabels = "Tregs", yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, xLabelsAngle=90, zlim = c(-1,1), main = paste("Module-trait relationships"))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(design), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, xLabelsAngle=90, zlim = c(-1,1), main = paste("Module-trait relationships"))
dev.off()
save(datExpr, moduleColors, design, file = "GSE13355_WGCNA.rda")

load("GSE13355_WGCNA.rda")
#提取疾病正模块基因
TOM = TOMsimilarityFromExpr(datExpr, power = softPower, networkType = "unsigned")

#TOM1 = TOM[1:round((nrow(TOM)/2)), 1:round((nrow(TOM)/2))]
#TOM2 = TOM[(round((nrow(TOM)/2)+1)):nrow(TOM), 1:round((nrow(TOM)/2))]
#TOM3 = TOM[1:round((nrow(TOM)/2)), (round((nrow(TOM)/2)+1)):nrow(TOM)]
#TOM4 = TOM[(round((nrow(TOM)/2)+1)):nrow(TOM), (round((nrow(TOM)/2)+1)):nrow(TOM)]

module = "brown" #megenta
inModule = (moduleColors==module)
probes = names(datExpr)
modProbes = probes[inModule]
modGenes = modProbes
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
##
#modTOM1 = TOM1[inModule[1:round((nrow(TOM)/2))], inModule[1:round((nrow(TOM)/2))]]
#modTOM2 = TOM2[inModule[(round((nrow(TOM)/2)+1)):nrow(TOM)], inModule[1:round((nrow(TOM)/2))]]
#modTOM3 = TOM3[inModule[1:round((nrow(TOM)/2))], inModule[(round((nrow(TOM)/2)+1)):nrow(TOM)]]
#modTOM4 = TOM4[inModule[(round((nrow(TOM)/2)+1)):nrow(TOM)], inModule[(round((nrow(TOM)/2)+1)):nrow(TOM)]]
#modTOM = rbind(cbind(modTOM1,modTOM3), cbind(modTOM2,modTOM4))
#dimnames(modTOM) = list(modProbes, modProbes)

#for(i in 1:(ncol(MEs)-1)){
#  module = labels2colors(i)
#  inModule = is.finite(match(moduleColors, module))
#  probes = names(datExpr)
#  modProbes = probes[inModule]
#  modGenes = modProbes
#  modTOM1 = TOM1[inModule[1:round((nrow(TOM)/2))], inModule[1:round((nrow(TOM)/2))]]
#  modTOM2 = TOM2[inModule[(round((nrow(TOM)/2)+1)):nrow(TOM)], inModule[1:round((nrow(TOM)/2))]]
#  modTOM3 = TOM3[inModule[1:round((nrow(TOM)/2))], inModule[(round((nrow(TOM)/2)+1)):nrow(TOM)]]
#  modTOM4 = TOM4[inModule[(round((nrow(TOM)/2)+1)):nrow(TOM)], inModule[(round((nrow(TOM)/2)+1)):nrow(TOM)]]
#  modTOM = rbind(cbind(modTOM1,modTOM3), cbind(modTOM2,modTOM4))
#  dimnames(modTOM) = list(modProbes, modProbes)
#}

#cyt = exportNetworkToCytoscape(modTOM, edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
#                               nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""), weighted = TRUE,
#                               threshold = 0.02, nodeNames = modProbes, altNodeNames = modGenes, nodeAttr = moduleColors[inModule])

##模块内连接度(节点基因与Hub基因)
module = "brown"
probes = names(datExpr)
inModule = (moduleColors==module)
modProbes = probes[inModule]
IMConn = softConnectivity(datExpr[, modProbes], power=softPower)

dat1 = datExpr[inModule]
datExpr_IMConn <- data.frame(IMConn, t(dat1))
#write.table(datExpr_IMConn, file=paste("Intramodule_connectivity_", module, ".txt"), sep="\t")

out = cbind(modProbes, IMConn)
colnames(out) = c("gene", "connectivity")
out = out[order(as.numeric(out[, 2]), decreasing=TRUE), ]
out = as.data.frame(out)
out$connectivity = as.numeric(out$connectivity)
write.table(out, file=paste(module, "-module-gene-GSE53431.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)


#xxx差异基因富集分析GSEA，用在线软件做-----------
load('GSE42632.rda')
load("step1_DEG_1.rda")
#raw_data1 = exprSet4
#raw_data1 = raw_data1[,samplefile4$Sample_Name]
#exprSet1 = raw_data1[nrDEG1_all[nrDEG1_all$group != "no sig",]$gene_symbol,]
#write.table(exprSet1, file = "差异基因表达谱.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

#xxx差异基因富集点图---------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(rvcheck)
library(ggplot2)

symbol2entrezid <- bitr(nrDEG1_all$gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
nrDEG_risk <- merge(nrDEG1_all, symbol2entrezid, by.x = "gene_symbol", by.y = "SYMBOL")
write.table(symbol2entrezid, file = "symbol2entrezid.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

#富集分析1
gene = nrDEG_risk[nrDEG_risk$group != "no sig",]$ENTREZID
##ALL
res_GO = enrichGO(gene = gene, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", pAdjustMethod = "BH",  pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE)
dotplot(res_GO, showCategory = 10, split = "ONTOLOGY") + facet_grid(ONTOLOGY~., scale = "free")

res_GO_BP = enrichGO(gene = gene, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",  pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE)
res_GO_BP_simp <- simplify(res_GO_BP, cutoff=0.7, by = "p.adjust", select_fun = min)                  
dotplot(res_GO_BP_simp, showCategory = 10)
write.table(as.data.frame(res_GO_BP_simp@result), file = "res_GO_BP_simp.txt", quote = FALSE, row.names=FALSE, col.names=TRUE)
BP = res_GO_BP_simp@result

res_GO_CC = enrichGO(gene = gene, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", pAdjustMethod = "BH",  pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE)
res_GO_CC_simp <- simplify(res_GO_CC, cutoff=0.7, by = "p.adjust", select_fun = min)                  
dotplot(res_GO_CC_simp, showCategory = 10)
write.table(as.data.frame(res_GO_CC_simp@result), file = "res_GO_CC_simp.txt", quote = FALSE, row.names=FALSE, col.names=TRUE)
CC = res_GO_CC_simp@result

res_GO_MF = enrichGO(gene = gene, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", pAdjustMethod = "BH",  pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE)
res_GO_MF_simp <- simplify(res_GO_MF, cutoff=0.7, by = "p.adjust", select_fun = min)                  
dotplot(res_GO_MF_simp, showCategory = 10)
write.table(as.data.frame(res_GO_MF_simp@result), file = "res_GO_MF_simp.txt", quote = FALSE, row.names=FALSE, col.names=TRUE)
MF = res_GO_MF_simp@result

res_KEGG <- enrichKEGG(gene = nrDEG_risk$ENTREZID, organism = "human", keyType = "kegg", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, use_internal_data = FALSE)
res_KEGG <- setReadable(res_KEGG, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
write.table(as.data.frame(res_KEGG@result), file = "res_KEGG_simp.txt", quote = FALSE, row.names=FALSE, col.names=TRUE)
dotplot(res_KEGG, showCategory = 10)
KEGG = res_KEGG@result


#
#三者取交集的基因-----------
Intersect <- read.table("Intersect.txt", header=TRUE, row.names=NULL, quote="", sep="\t", check.names = FALSE)
load("GSE13355.rda")

Inter_exp = exprSet[Intersect$Intersect,]
save(Inter_exp, SampleFile, file = "intersect.rda")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(rvcheck)
library(ggplot2)

symbol2entrezid <- bitr(Intersect$Intersect, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
nrDEG_risk <- merge(Intersect, symbol2entrezid, by.x = "Intersect", by.y = "SYMBOL")


#xxx三者交集基因富集分析1----
gene = nrDEG_risk$ENTREZID
##ALL
res_GO = enrichGO(gene = gene, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", pAdjustMethod = "BH",  pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE)
dotplot(res_GO, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY~., scale = "free")
#dev.off()

res_KEGG <- enrichKEGG(gene = nrDEG_risk$ENTREZID, organism = "human", keyType = "kegg", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, use_internal_data = FALSE)
res_KEGG <- setReadable(res_KEGG, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
write.table(as.data.frame(res_KEGG@result), file = "res_KEGG_simp.txt", quote = FALSE, row.names=FALSE, col.names=TRUE)
dotplot(res_KEGG, showCategory = 10)
KEGG = res_KEGG@result

#机器学习算法------------------
load("intersect.rda")
#特征选择
#SampleFile = samplefile4
SampleFile = SampleFile[order(SampleFile$Sample_Group,decreasing = F),]
#SampleFile$Sample_Group = gsub("SLE","1",SampleFile$Sample_Group)
#SampleFile$Sample_Group = gsub("HC","0",SampleFile$Sample_Group)
SampleFile$Sample_Group = as.factor(SampleFile$Sample_Group)
exprSet = Inter_exp
exprSet = exprSet[,match(SampleFile$Sample_Name,colnames(exprSet))]
exprSet = as.data.frame(t(exprSet))
exprSet$Sample_Group = SampleFile$Sample_Group; exprSet = exprSet[,c(ncol(exprSet),1:ncol(exprSet)-1)]
data = exprSet

#Boruta----------------------------------------------
library(Boruta)
set.seed(1234)
boruta.train <- Boruta(Sample_Group ~ ., data=data, doTrace=2, maxRuns=1000)
print(boruta.train)

final.boruta <- TentativeRoughFix(boruta.train)
print(final.boruta)

#
plot(final.boruta, xlab="", xaxt="n")
lz <- lapply(1:ncol(final.boruta$ImpHistory), function(i) final.boruta$ImpHistory[is.finite(final.boruta$ImpHistory[, i]), i])
names(lz) <- colnames(final.boruta$ImpHistory)
Labels <- sort(sapply(lz, median)) 
axis(side=1, las=2, labels=names(Labels), at=1:ncol(final.boruta$ImpHistory), cex.axis=0.5, font=4)
dev.off()

predBoruta_g = getSelectedAttributes(final.boruta, withTentative = F)
predBoruta_g
write.table(predBoruta_g, file = "predBoruta_g.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = T)
#LASSO&SVM-RFE-------------------------------------------------------
x <- as.matrix(data[,-1])
y <- data$Sample_Group

#' @LASSO-logistic-Algorithm
library(glmnet)
set.seed(1234)
fit = glmnet(x, y, family = "binomial",alpha = 1,lambda = NULL)

plot(fit, xvar = "dev", label = TRUE)

cvfit = cv.glmnet(x, y, family = "binomial", type.measure = "class")
plot(cvfit)
dev.off()

cvfit$lambda.min#二选一，下面也改
#cvfit$lambda.1se

#' @LASSO选出来的特征如下
myCoefs <- coef(cvfit, s="lambda.min");
lasso_fea<-myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
lasso_fea<-lasso_fea[-1];
lasso_fea
write.table(lasso_fea, "lasso.txt", sep = "\t", row.names=FALSE)

#' @SVM-RFE-Algorithm
library(parallel)
library(doParallel)
cl <- makeCluster(4)#detectCores()
registerDoParallel(cl)
library(caret)
library(randomForest)
set.seed(1234)
subsets = seq(1, ncol(x), by = 1)
caretFuncs$summary = defaultSummary
rfe_ctrl = rfeControl(functions=rfFuncs, method="repeatedcv", number=5, repeats=3, verbose=FALSE, returnResamp="final", allowParallel=TRUE)
#functions=rfFuncs/caretFuncs
system.time(
  rfe.train <- rfe(x, y, sizes=subsets, rfeControl=rfe_ctrl, method = "rf")#svmLinear
)
print(rfe.train)
stopCluster(cl)
registerDoSEQ()

plot(rfe.train, type=c("g", "o"), cex = 0.1, col = "firebrick3", lwd=3)
plot(rfe.train, type=c("g", "o"), cex = 0.1, col = "green", lwd=3)
#plot(rfe.train, type = c("o", "g"), cex = 1.0, col = 1:length(subsets))
dev.off()

#' @SVM-RFE选出来的特征如下
rfe_fea <- predictors(rfe.train)
rfe_fea
write.table(rfe_fea, "rfe_fea.txt", sep = "\t", row.names=FALSE)

com = intersect(lasso_fea, rfe_fea)
com
# "BCL2A1"  "MTHFD2"  "LYRM2"   "FASTKD3" "ACACA" 


#RF----
library(randomForest)
set.seed(1234)

#读取输入文件
data = Inter_exp
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))  
#View(group)

#随机森林树
rf=randomForest(as.factor(group)~., data=data, ntree=500)
#pdf(file="forest.pdf", width=6, height=6)
plot(rf, main="Random forest", lwd=2)
dev.off()

#找出误差最小的点
optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(as.factor(group)~., data=data, ntree=optionTrees)

#查看基因的重要性
importance=importance(x=rf2)

#调整RF的GiNi系数
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes=names(rfGenes[rfGenes>=1])     #挑选重要性评分≥1的基因
#rfGenes=names(rfGenes[1:30])         #挑选重要性评分最高的30个基因
write.table(rfGenes, file="RFGenes.txt", sep="\t", quote=F, col.names=F, row.names=F)

#绘制基因的重要性图
pdf(file="geneImportance.pdf", width=6, height=6)
varImpPlot(rf2, main="")
dev.off()


#特征基因的差异表达-------------------
setwd("F:/GSE/PP/GSE13355.1.5/")
library(ggplot2) 
library(RColorBrewer) 
library(ggpubr)
library(preprocessCore)
library(e1071)
library(parallel)
library(data.table)

load("intersect.rda")
tcga_tmp = Inter_exp
rownames(SampleFile)=SampleFile$Sample_Name

TSPYL2_exp = tcga_tmp[rownames(tcga_tmp) == "S100A9",]
TSPYL2_exp = as.data.frame(t(TSPYL2_exp))

TSPYL2_exp$Sample_Name = rownames(TSPYL2_exp)

TSPYL2_exp1 = as.data.frame(TSPYL2_exp[rownames(SampleFile),])


data2 = merge(TSPYL2_exp1, SampleFile, by="Sample_Name")
data2$Sample_Group = as.factor(data2$Sample_Group)
compaired <- list(c("PP", "HC"))  #设置比较组别
palette<-c(brewer.pal(7,"Set2")[c(1,2,4,5)])   #颜色调用 
colnames(data2)
ggboxplot(data2, 
          x = "GSE13355 Sample_Group", y = "S100A9", 
          fill = "Sample_Group", palette = palette, 
          add = "jitter", size = 0.5)+   #添加抖动的散点 
  stat_compare_means(comparisons = compaired, 
                     method = "wilcox.test",   #设置统计方法 
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns"))) 
dev.off()


library(pROC)
library(ggplot2)


BRCA=data2
res<-roc(Sample_Group~S100A9,data=BRCA,aur=TRUE,
         ci=TRUE, # 显示95%CI
         #percent=TRUE, # 是否需要以百分比显示
         smooth=F,# 是否平滑曲线
         levels=c('PP','HC'),direction=">" #设置分组方向
)

p<- ggroc(res, color ="blue",legacy.axes = TRUE)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype=4)+
  theme_bw() + # 设置背景
  ggtitle("GSE13355 ROC Curve")+
  theme(plot.title = element_text(hjust = 0.5,size = 16),
        axis.text=element_text(size=12,colour = "black"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

p+annotate("text",x=0.75,y=0.25,label=paste("AUC = ", round(res$auc,3)))+
  annotate("text",x=0.75,y=0.20,label=paste("95%CI: ", round(res$ci[1],3),'-',round(res$ci[3],3)))

dev.off()

# 验证基因的表达在其他数据集---------------------
load("F:/GSE/PP/GSE75343/GSE75343.rda")
library(ggplot2) 
library(RColorBrewer) 
library(ggpubr)
library(preprocessCore)
library(e1071)
library(parallel)
library(data.table)

#exprSet1 = read.table("GSE89408_GEO_count_matrix_rename.txt", header=T, sep="\t", check.names=F,row.names =1)
#samplefile1 = fread("samplefile_GSE89408.txt", header = TRUE, data.table = FALSE)
exprSet1 = exprSet[,SampleFile$Sample_Name]

tcga_tmp = exprSet1
rownames(SampleFile)=SampleFile$Sample_Name

TSPYL2_exp = tcga_tmp[rownames(tcga_tmp) == "S100A9",]#EIF4EBP1/TYMP
TSPYL2_exp = as.data.frame(t(TSPYL2_exp))

TSPYL2_exp$Sample_Name = rownames(TSPYL2_exp)

TSPYL2_exp1 = as.data.frame(TSPYL2_exp[rownames(SampleFile),])


data2 = merge(TSPYL2_exp1, SampleFile, by="Sample_Name")
data2$Sample_Group = as.factor(data2$Sample_Group)
compaired <- list(c("PP", "HC"))  #设置比较组别
palette<-c(brewer.pal(7,"Set2")[c(1,2,4,5)])   #颜色调用 
colnames(data2)
ggboxplot(data2, 
          x = "GSE14905 Sample_Group", y = "S100A9",         #Sample_Group
          fill = "Sample_Group", palette = palette, 
          add = "jitter", size = 0.5)+   #添加抖动的散点 
  stat_compare_means(comparisons = compaired, 
                     method = "wilcox.test",   #设置统计方法 
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns"))) 
dev.off()

library(pROC)
library(ggplot2)

#SampleFile = samplefile1
SampleFile = SampleFile[order(SampleFile$Sample_Group,decreasing = F),]
#SampleFile$Sample_Group = gsub("SLE","1",SampleFile$Sample_Group)
#SampleFile$Sample_Group = gsub("HC","0",SampleFile$Sample_Group)
SampleFile$Sample_Group = as.factor(SampleFile$Sample_Group)
exprSet = exprSet1
exprSet = exprSet[,match(SampleFile$Sample_Name,colnames(exprSet))]
exprSet = as.data.frame(t(exprSet))
exprSet$Sample_Group = SampleFile$Sample_Group; exprSet = exprSet[,c(ncol(exprSet),1:ncol(exprSet)-1)]
data = exprSet

BRCA=data2
res<-roc(Sample_Group~S100A9,data=BRCA,aur=TRUE,
         ci=TRUE, # 显示95%CI
         #percent=TRUE, # 是否需要以百分比显示
         smooth=F,# 是否平滑曲线
         levels=c('PP','HC'),direction=">" #设置分组方向
)

p<- ggroc(res, color ="blue",legacy.axes = TRUE)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype=4)+
  theme_bw() + # 设置背景
  ggtitle("GSE109248 ROC Curve")+     #S100A9
  theme(plot.title = element_text(hjust = 0.5,size = 16),
        axis.text=element_text(size=12,colour = "black"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

p+annotate("text",x=0.75,y=0.25,label=paste("AUC = ", round(res$auc,3)))+
  annotate("text",x=0.75,y=0.20,label=paste("95%CI: ", round(res$ci[1],3),'-',round(res$ci[3],3)))

dev.off()

#去批次后，分型--------------------------
#训练集尝试一   
setwd("F:/GSE/PP/GSE13355.1.5/")
load("step1_combat.rda")
Intersect <- read.table("Intersect.txt", header=TRUE, row.names=NULL, quote="", sep="\t", check.names = FALSE)

Inter_exp1 = combat_exprSet[,combat_samplefile$Sample_Group == "PP"]

Inter_exp = na.omit(Inter_exp1[Intersect$Intersect,])

#共识聚类
library(ConsensusClusterPlus)
ccp.input = as.matrix(Inter_exp)
ccp.input = sweep(as.matrix(ccp.input), 1, apply(as.matrix(ccp.input), 1, median, na.rm=TRUE))
#ccp.input = t(scale(t(ccp.input)))#Z

results = ConsensusClusterPlus(ccp.input,
                               maxK=6, reps=1000, pItem=0.8, pFeature=1, clusterAlg="hc",
                               innerLinkage = "ward.D", finalLinkage = "ward.D",
                               distance="pearson", seed=1234, plot="pdf", writeTable=TRUE)
#pam euclidean
#hc pearson/spearman
#km  euclidean

Kvec = 2:6
x1 = 0.1; x2 = 0.9
PAC = rep(NA,length(Kvec))
names(PAC) = paste("K=", Kvec, sep="")
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)}
rcc.ind = Kvec[which.min(PAC)]
rcc.ind
#group
group <- results[[rcc.ind]]$consensusClass; table(group)
samorder = sort(results[[rcc.ind]]$consensusClass)
#3
samplefile <- data.frame("Subtype" = ifelse(samorder == 1,"clusterA", ifelse(samorder == 2, "clusterB", "clusterC")), 
                         row.names = names(samorder), check.names = F, stringsAsFactors = F)
#2
samplefile <- data.frame("Subtype" = ifelse(samorder == 1,"cluster1", "cluster2"),
                         row.names = names(samorder), check.names = F, stringsAsFactors = F)

samplefile = as.data.frame(cbind(Sample_Name = rownames(samplefile), samplefile))

samplefile$Subtype = factor(samplefile$Subtype)
table(samplefile$Subtype)
exprSet = as.data.frame(ccp.input[, match(samplefile$Sample_Name, colnames(ccp.input))])
save(results, exprSet, combat_exprSet, samplefile,file = "step2_hc.rda")

load("step2_hc.rda")

#XXXNMF----
library(NMF)
library(Biobase)
library(doParallel)
nmf.input = as.matrix(Inter_exp)

cl <- makeCluster(detectCores()); registerDoParallel(cl)
system.time(res.nmf <- NMF::nmf(nmf.input, 2:6, nrun=10, method="brunet", .opt='p', .pbackend=NULL, seed=123456))
plot(res.nmf)
stopCluster(cl); registerDoSEQ(); rm(cl)
#save(res.nmf, file='res.nmf.RData')
#load(file='res.nmf.RData')

coph = res.nmf$measures$cophenetic
coph_diff <- NULL
for (i in 2:length(coph)) {
  coph_diff <- c(coph_diff, abs(coph[i-1]-coph[i]))
}
k.best <- which.max(coph_diff)+1
#k.best = 3

res.nmf1 <- nmf(nmf.input, rank=k.best, method="brunet", nrun=100, seed=123456)
index <- extractFeatures(res.nmf1)
sig.order <- unlist(index)
rownames(nmf.input)[sig.order]
group <- predict(res.nmf1); table(group)

samplefile = as.data.frame(cbind(names(group), group))
colnames(samplefile) = c("Sample_Name","Subtype")
rownames(samplefile) = NULL
samplefile = samplefile[order(samplefile$Subtype), ]
samplefile$Subtype = ifelse(samplefile$Subtype == 1,"clusterA", "clusterB")
samplefile$Subtype = factor(samplefile$Subtype)

group = factor(samplefile$Subtype)
names(group) = samplefile$Sample_Name

sample.order <- names(group[order(group)])
jco <- c("#2874C5","#EABF00","#C6524A","#868686")
#consensusmap
consensusmap(res.nmf1, labRow=NA, labCol=NA,
             #annCol=data.frame("cluster"=group[colnames(nmf.input)]),
             annColors=list(cluster=c("1"=jco[1],"2"=jco[2],"3"=jco[3],"4"=jco[4])))#tracks=NA,Rowv=NA
dev.off()

exprSet = as.data.frame(nmf.input[, match(samplefile$Sample_Name, colnames(nmf.input))])
save(exprSet, samplefile, Inter_exp, file = "step2_nmf.rda")


#PCA------------------------------------------------------
setwd("F:/GSE/PP/GSE13355.1.5/")
load("step2_hc.rda")
library(ggplot2)
library(plyr)
#library(devtools)
library(ggord)
data = exprSet
rownames(data) = gsub("-", " ", rownames(data))
rownames(data) = gsub(" ", "", rownames(data))
source('geom_ord_ellipse.R') #该文件位于当前文件夹
pca.results <- prcomp(t(data), center = TRUE, scale. = FALSE)
#install.packages('devtools') 
#定义足够多的颜色，用于展示分组
#mycol <- c("#223D6C","#D20A13","#088247","#FFD121","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
mycol <- c("#FFD121","#7A142C","#223D6C")

#用ggord画基本PCA图
ggord(pca.results, grp_in = samplefile$Subtype, repel=TRUE,
      alpha = 0.6, ellipse_pro = 0.95, size = 2,
      cols = mycol[1:length(unique(samplefile$Subtype))],
      arrow=0, vec_ext = 0, veccol="brown", txt=0) +
  theme(panel.grid =element_blank()) +
  geom_ord_ellipse(ellipse_pro = .95, color='darkgrey', size=0.5, lty=2 ) +
  geom_ord_ellipse(ellipse_pro = .98, size=0.5, lty=2)
dev.off()

#heatmap-----------------------------------------------------------------
library(pheatmap)
data1 = exprSet
#clin = new_samplefile[new_samplefile$Sample_Name %in% com_sample, ]
#rownames(clin) = clin$Sample_Name
annCol <- data.frame("Gene cluster" = samplefile$Subtype, row.names = samplefile$Sample_Name, check.names = F, stringsAsFactors = F)
#annCol <- cbind.data.frame(annCol,clin[rownames(annCol),])
#annCol <- annCol[,-2]
#annRow <- data.frame("Signature gene" = boruta.all$direct, row.names = boruta.all$gene, check.names = F, stringsAsFactors = F)
annColors <- list("cluster" = c("A" = "#008ECB", "B" = "#EA921D", "C" = "#D14039")) #, "Signature gene" = c("A" = "#D14039", "B" = "#008ECB"))
plotdata <- data1[, rownames(annCol)]
plotdata <- t(scale(t(plotdata)))
plotdata[plotdata > 1] <- 1
plotdata[plotdata < -1] <- -1

pheatmap(plotdata, color=colorRampPalette(c("#343493", "white", "#C24A45"))(80), border_color=NA,
         annotation_col=annCol, annotation_colors = annColors,
         #scale="row",
         cluster_row=, cluster_cols=F, show_colnames=F, show_rownames=T,
         clustering_distance_cols = "euclidean", clustering_method = "complete")
dev.off()

#免疫细胞浸润分析----------------
setwd("F:/GSE/PP/GSE13355.1.5/")
load("step2_hc.rda")
load("step1_combat.rda")
library(ggplot2) 
library(RColorBrewer) 
library(ggpubr)
library(preprocessCore)
library(e1071)
library(parallel)
source("CIBERSORT.R") 
combat_exprSet1 = combat_exprSet[ ,samplefile$Sample_Name]
sig_matrix <- "LM22.txt"
write.table(cbind(rownames(combat_exprSet1), combat_exprSet1), file="cibersort.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
mixture_file = 'cibersort.txt'
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE)
res_cibersort = read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F) 
res_cibersort = res_cibersort[, -23:-25]
res_cibersort = res_cibersort[,colSums(res_cibersort)>0]
save(res_cibersort, file = "res_cibersort_第一次分型hc.rda")

load("res_cibersort_第一次分型hc.rda")
data1 = as.data.frame(res_cibersort)
data1$Sample_Name = rownames(data1)
data2 = merge(data1, samplefile, by="Sample_Name")
write.table(data2, file="data2.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

samplefile$Subtype = factor(samplefile$Subtype, levels = c("cluster1","cluster2"))
res_cibersort$Sample_Name = rownames(res_cibersort)
samplefile$Sample_Name = rownames(samplefile)
res_cibersort = merge(samplefile[,c("Sample_Name","Subtype")], res_cibersort, by="Sample_Name")
library(ggpubr)
library(tidyr)
long_dat <- pivot_longer(res_cibersort, cols=3:ncol(res_cibersort), names_to="Gene", values_to = "Expression")

#箱线图-
ggboxplot(long_dat, x="Gene", y="Expression", fill="Subtype",
          xlab="",
          ylab="Gene expression",
          legend.title="Gene",
          width=0.8,
          palette = c("#ECE644", "#E6756F") )+
  rotate_x_text(50)+
  stat_compare_means(aes(group=Subtype),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
dev.off()

#免疫检查点----
library(ggplot2) 
library(RColorBrewer) 
library(ggpubr)
library(preprocessCore)
library(e1071)
library(parallel)
library(data.table)
rm(list = ls())
setwd("F:/GSE/PP/GSE13355.1.5/")
load("step2_hc.rda")
load("step1_combat.rda")
tcga_tmp = combat_exprSet
rownames(samplefile)=samplefile$Sample_Name

TSPYL2_exp = tcga_tmp[rownames(tcga_tmp) == "VTCN1",]
TSPYL2_exp = as.data.frame(t(TSPYL2_exp))#TSPYL2_exp = as.data.frame(t(TSPYL2_exp))
rownames(TSPYL2_exp) = "VTCN1"
TSPYL2_exp = as.data.frame(t(TSPYL2_exp))
TSPYL2_exp$Sample_Name = rownames(TSPYL2_exp)

TSPYL2_exp1 = as.data.frame(TSPYL2_exp[rownames(samplefile),])

data2 = merge(TSPYL2_exp1, samplefile, by="Sample_Name")
data2$Subtype = as.factor(data2$Subtype)
compaired <- list(c("cluster1", "cluster2"))  #设置比较组别
palette<-c(brewer.pal(7,"Set2")[c(1,2,4,5)])   #颜色调用 
colnames(data2)
ggboxplot(data2, 
          x = "Subtype", y = "VTCN1", 
          fill = "Subtype", palette = palette, 
          add = "jitter", size = 0.5)+   #添加抖动的散点 
  stat_compare_means(comparisons = compaired, 
                     method = "wilcox.test",   #设置统计方法 
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns"))) 
dev.off()



#XXXcell----
library(ConsensusClusterPlus)
library(xCell)
res_xcell <- xCellAnalysis(combat_exprSet1)
res_xcell <- as.data.frame(res_xcell)
table(rowSums(res_xcell)>0)
res_xcell1 <- res_xcell[-65:-67,]

data3 = t(as.data.frame(res_xcell1))
data3 <- as.data.frame(data3)
data3$Sample_Name = rownames(data3)
data4 = merge(data3, samplefile, by="Sample_Name")

write.table(data4, file = "data4.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
data5 = raw_data <- read.table(file = "data4.txt", header=T, sep="\t", check.names=F,row.names =1)

compaired <- list(c("clusterA", "clusterB"),c("clusterA","clusterC"),c("clusterB","clusterC") )  #设置比较组别
palette<-c(brewer.pal(7,"Set2")[c(1,2,4,5)])   #颜色调用 
colnames(data2)
ggboxplot(data2, 
          x = "Subtype", y = "Macrophages M2", 
          fill = "Subtype", palette = palette, 
          add = "jitter", size = 0.5)+   #添加抖动的散点 
  stat_compare_means(comparisons = compaired, 
                     method = "wilcox.test",   #设置统计方法 
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns"))) 
dev.off()

#两型之间的差异基因----
setwd("F:/GSE/PP/GSE13355.1.5/")
load("step2_hc.rda")
load("step1_combat.rda")
library(limma)
new_exprSet = as.data.frame(normalizeBetweenArrays(as.matrix(combat_exprSet), method = "quantile")) 
new_exprSet = new_exprSet[,samplefile$Sample_Name]
samplefile$Subtype = factor(samplefile$Subtype, levels = c("cluster1", "cluster2"))
design = model.matrix(~ 0 + samplefile$Subtype)
colnames(design) = levels(samplefile$Subtype)
contrast.matrix = makeContrasts("cluster2-cluster1", levels = design)
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

nrDEG1_all <- nrDEG1[which(nrDEG1$adj.P.Val < 0.05 & abs(nrDEG1$logFC) >= 1), ]
nrDEG1_all[which(nrDEG1_all$logFC > 0), "group"] <- "up"
nrDEG1_all[which(nrDEG1_all$logFC < 0), "group"] <- "down"
table(nrDEG1_all$group)
write.table(nrDEG1_all, file = "差异基因1hc.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
save(combat_exprSet, samplefile, nrDEG1, nrDEG1_all, file = "step3_DEG0_hc.rda")


#XXX富集分析----
library(msigdbr)
library(dplyr)
library(data.table)
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)

#查看msigdbr包里自带的物种H:
#hallmark gene sets
#C1: positional gene sets
#C2: curated gene sets
#C3: motif gene sets
#C4: computational gene sets
#C5: GO gene sets
#C6: oncogenic signatures
#C7: immunologic signatures
# kegg 在C2,GO 在C5

#msigdbr_show_species()#可查看物种名
h <- msigdbr(species = "Homo sapiens", category = "C7") 

# 示例数据表达矩阵的基因名是gene symbol，这里就选gene_symbol。
# 如果你的表达矩阵以ENTREZ ID作为基因名，就把下面这段的gene_symbol换成entrez_gene
h <- select(h, gs_name, gene_symbol) %>% #或entrez_gene
  as.data.frame %>% 
  split(., .$gs_name) %>% 
  lapply(., function(x)(x$gene_symbol)) #或entrez_gene

# 在每个geneset里面去掉重复的基因
gs <- lapply(h, unique)

# 接下来去掉那些在两个或更多个pathways里出现过的genes
count <- table(unlist(gs))
keep <- names(which(table(unlist(gs)) < 2))
gs <- lapply(gs, function(x) intersect(keep, x))

# 过滤之后，很多pathway一个gene都不剩了，去掉这些
gs <- gs[lapply(gs, length) > 0]
#save(gs, file = "C5.gs.RData")

raw_data1 = new_exprSet
raw_data1 = raw_data1[,rownames(samplefile)]
exprSet1 = raw_data1[nrDEG1_all[nrDEG1_all$group != "no sig",]$gene_symbol,]
gsym.expr <- exprSet1
head(gsym.expr)
gsva_es <- gsva(as.matrix(gsym.expr), gs)# 这一句就完成了GSVA分析

#通路的差异表达分析--
group_list <- samplefile
design <- model.matrix(~ 0 + factor(samplefile$Subtype))
colnames(design) <- levels(factor(samplefile$Subtype))
rownames(design) <- colnames(gsva_es)
head(design)

contrast.matrix <- makeContrasts(cluster1-cluster2, levels = design)

fit <- lmFit(gsva_es, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(x)

pathway <- str_replace(row.names(x), "C2", "")
df <- data.frame(ID = pathway, score = x$t)#输出t值，用做FigureYa39bar的输入数据
head(df)
#按照score的值分组
cutoff <- 1#这里可以修改
df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3))
#View(df)

#按照score排序
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)

write.table(sortdf, file = "sortdf_KEGG.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = T)
ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c('palegreen3', 'snow3', '#EA921D'), guide = FALSE) + #如果只有13两组，请将颜色“snow3"删掉
  
  #画2条虚线
  geom_hline(yintercept = c(-cutoff,cutoff), 
             color="white",
             linetype = 2, #画虚线
             size = 0.4) + #线的粗细
  
  #写label
  geom_text(data = subset(df, score < 0),
            aes(x=ID, y= 0, label= paste0(" ", ID), color = group),#bar跟坐标轴间留出间隙
            size = 3, #字的大小
            hjust = "inward" ) +  #字的对齐方式
  geom_text(data = subset(df, score > 0),
            aes(x=ID, y= -0.1, label=ID, color = group),
            size = 3, hjust = "outward") +  
  scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) + #如果只有13两组，请将颜色“snow3"删掉
  
  xlab("") +ylab("t value of GSEA score")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 0.6)) + #边框粗细
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) #去除y轴
dev.off()

#hub基因和免疫系细胞的相关性----
setwd("F:/GSE/PP/GSE13355.1.5/")
load("step2_hc.rda")
load("res_cibersort_第一次分型hc.rda")
library(ggplot2)
library(ggpubr)
library(ggExtra)

res_cibersort1 = as.data.frame(t(res_cibersort))
#空格不识别，需要写出去一个
write.table(res_cibersort1, file = "res_cibersort_hc.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = T)
#res_cibersort_km.txt
#CIBERSORT-Results
res_cibersort1 <- read.table("res_cibersort_hc.txt", header=TRUE, row.names=1, quote="", sep="\t", check.names = FALSE)
rownames(res_cibersort1)
gene1="S100A9"             #第一个基因名字 EIF4EBP1/TYMP
gene2="Neutrophils"  


x=as.numeric(exprSet[gene1,])
y=as.numeric(res_cibersort1[gene2,])

df1=as.data.frame(cbind(x,y))
corT=cor.test(x,y,method="spearman")
cor=corT$estimate
pValue=corT$p.value
p1=ggplot(df1, aes(x, y)) + 
  xlab(gene1)+ylab(gene2)+
  geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =x, y =y))
p1

p2=ggMarginal(p1, type = "density", xparams = list(fill = "green"),yparams = list(fill = "red"))

p2
dev.off()


#构建评分模型------------------
setwd("F:/GSE/PP/GSE13355.1.5/")
load("step2_hc.rda")
#brouta
exprSet = exprSet[, match(samplefile$Sample_Name, colnames(exprSet))]
rownames(exprSet) = gsub("-", "_", rownames(exprSet))
outTab <- NULL
#用pearson做基因与分组的相关性
for (i in rownames(exprSet)) 
{
  tmp <- as.numeric(exprSet[i, samplefile$Sample_Name])
  cor.res <- cor.test(tmp, as.numeric(samplefile$Subtype), method = "pearson")
  outTab <- rbind.data.frame(outTab,
                             data.frame(gene = i,
                                        r = cor.res$estimate,
                                        p = cor.res$p.value,
                                        stringsAsFactors = F),
                             stringsAsFactors = F)
}

# 按相关性正负来分类
#outTab = outTab[outTab$p < 0.05,]
outTab$direct <- ifelse(outTab$r > 0, "A","B") # 正相关标为A，否则标为B
outTab <- outTab[order(outTab$r, decreasing = T),]
table(outTab$direct)

#一起降维
model_dat = t(scale(t(exprSet[outTab$gene,])))#z标准化
model_dat = as.data.frame(t(model_dat))
Class = factor(samplefile$Subtype)
model_dat = cbind.data.frame(Class, model_dat)

library(Boruta)
set.seed(1234)
boruta.train <- Boruta(Class ~ ., data=model_dat, doTrace=2, maxRuns=1000)
print(boruta.train)

final.boruta <- TentativeRoughFix(boruta.train)
print(final.boruta)

#
plot(final.boruta, xlab="", xaxt="n")
lz <- lapply(1:ncol(final.boruta$ImpHistory), function(i) final.boruta$ImpHistory[is.finite(final.boruta$ImpHistory[, i]), i])
names(lz) <- colnames(final.boruta$ImpHistory)
Labels <- sort(sapply(lz, median))
axis(side=1, las=2, labels=names(Labels), at=1:ncol(final.boruta$ImpHistory), cex.axis=0.5, font=4)
dev.off()

predBoruta = getSelectedAttributes(final.boruta, withTentative = F)
#write.table(predBoruta, "predBoruta", sep = "\t", row.names=FALSE)

boruta.all = outTab[outTab$gene %in% predBoruta, ]
table(boruta.all$direct)
save(boruta.all, file = "boruta.all.rda")
#热图------------------------------------------------------
library(pheatmap)
#heatmap
data1 = exprSet
#rownames(data1)=gsub("-", "_", rownames(data1))
#clin = new_new_samplefile[new_new_samplefile$Sample_Name %in% com_sample, ]
#rownames(clin) = clin$Sample_Name
annCol <- data.frame("Gene cluster" = samplefile$Subtype, row.names = samplefile$Sample_Name, check.names = F, stringsAsFactors = F)
#annCol <- cbind.data.frame(annCol,clin[rownames(annCol),])
#annCol <- annCol[,-2]
annRow <- data.frame("Signature gene" = boruta.all$direct, row.names = boruta.all$gene, check.names = F, stringsAsFactors = F)
annColors <- list("Gene cluster" = c("cluster1" = "#008ECB", "cluster2" = "#EA921D"), "Signature gene" = c("A" = "#D14039", "B" = "#008ECB"))

plotdata <- t(scale(t(data1[rownames(annRow), rownames(annCol)])))
plotdata[plotdata > 2] <- 2
plotdata[plotdata < -2] <- -2

pheatmap(plotdata, color=colorRampPalette(c("#343493", "white", "#C24A45"))(64), border_color=NA,
         annotation_col=annCol, annotation_row = annRow, annotation_colors = annColors,
         #scale="row",
         cluster_row=T, cluster_cols=T, show_colnames=F, show_rownames=F,
         clustering_distance_cols = "euclidean", clustering_method = "complete")
dev.off()

#箱示图----
stat_compare_means(comparisons=list(c("gene_clusterA","gene_clusterB")),
                   symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1),
                                    symbols=c("***", "**", "*", "ns")),
                   label="p.signif")

boruta.A = boruta.all[boruta.all$direct=="A",]
expr.A <- t(scale(t(data1[boruta.A$gene, ])))
pca.A <- prcomp(t(expr.A), scale = F, center = F)
pca1.A <- pca.A$x[,1]

boruta.B = boruta.all[boruta.all$direct=="B",]
expr.B <- t(scale(t(data1[boruta.B$gene, ])))
pca.B <- prcomp(t(expr.B), scale = F, center = F)
pca1.B <- pca.B$x[,1]

gene.score <- pca1.A - pca1.B
gene.outtab <- data.frame(samID = rownames(annCol),
                          pca1.A = pca1.A[rownames(annCol)],
                          pca1.B = pca1.B[rownames(annCol)],
                          gene.score = gene.score[rownames(annCol)], Sample_Group = annCol$`Gene cluster`, stringsAsFactors = F)
save(gene.outtab,file = "gene_outtab.rda")

ggplot(data = gene.outtab,aes(x = Sample_Group, y = gene.score, fill = Sample_Group))+
  scale_fill_manual(values = c("#A8DADB", "#BEB3EA", "#66B98E")) +
  geom_violin(alpha=0.4, position = position_dodge(width = .75),
              size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = TRUE, outlier.size = -1,
               color="black", lwd=0.8, alpha = 0.7)+ # 背景色透明化
  geom_point(shape = 21, size=2,
             position = position_jitterdodge(),
             color="black", alpha=1)+ # 边框线黑色
  theme_classic() +
  ylab(expression("gene score")) +
  xlab("Sample_Group")  +
  theme(axis.ticks = element_line(size=0.2,color="black"),
        axis.ticks.length = unit(0.2,"cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  stat_compare_means(method = "kruskal.test", label.y = max(gene.score))


#ROC----

library(pROC)
library(ggsci)
roc.list <- pROC::roc(gene.outtab$Sample_Group, gene.outtab$gene.score, ci =T)
roc.list$auc
roc.list$ci
roc_result <- coords(roc.list, "best")
#roc_result$threshold
roc_result

g3 <- ggroc(roc.list, size = 1.2,alpha=.6)
g3 <- g3 + geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed")
g3 <- g3 + scale_fill_manual(values = "#EA921D") + theme_classic()
g3 <- g3 + scale_color_aaas() + theme_classic()
g3 #geom_hline(yintercept =0.9411765)

plot.roc(gene.outtab$Sample_Group, gene.outtab$gene.score, percent=TRUE, smooth=F, col="#e55709")
legend("bottomright", legend=c("AUC=99.33%"), col=c("#ffdd00"), lwd=2, cex=0.5, text.width=80)
dev.off()

library(pROC)
rocobj <- plot.roc(gene.outtab$Sample_Group, gene.outtab$gene.score,
                   main="Confidence intervals", percent=TRUE,
                   ci=TRUE, # compute AUC (of AUC by default)
                   print.auc=TRUE) # print the AUC (will contain the CI)
ciobj <- ci.se(rocobj, # CI of sensitivity
               specificities=seq(0, 100, 5)) # over a select set of specificities
plot(ciobj, type="shape", col="#fbb034", alpha = 0.5) # plot as a blue shape
plot(ci(rocobj, of="thresholds", thresholds="best")) # add one threshold

#gene_score和response是否有差异
#药物响应---------------------------------------------------
library(ggplot2)
library(ggpubr)
setwd("F:/GSE/YWXY/GSE11903药物响应/") #Ustekinumab/Etanercept
#load("F:/GSE/药物响应/Etanercept")
combat_samplefile = SampleFile
combat_exprSet = exprSet
#combat_exprSet = new_exprSet[,-1]
combat_exprSet = combat_exprSet[,combat_samplefile$Sample_Name]
combat_samplefile$Sample_Group = gsub("YES", "response", combat_samplefile$Sample_Group)
combat_samplefile$Sample_Group = gsub("NO", "non_response", combat_samplefile$Sample_Group)

#先手动加入boruta.all
length(intersect(boruta.all$gene, rownames(combat_exprSet)))

setdiff(boruta.all$gene, rownames(combat_exprSet))
combat_exprSet1 = combat_exprSet[intersect(boruta.all$gene, rownames(combat_exprSet)), ]

boruta.A = intersect(boruta.all[boruta.all$direct == "A", ]$gene, rownames(combat_exprSet))
expr.A <- t(scale(t(combat_exprSet[boruta.A, ])))
pca.A <- prcomp(t(expr.A), scale = F, center = F)
pca1.A <- pca.A$x[,1]

boruta.B = intersect(boruta.all[boruta.all$direct == "B", ]$gene, rownames(combat_exprSet))
expr.B <- t(scale(t(combat_exprSet[boruta.B, ])))
pca.B <- prcomp(t(expr.B), scale = F, center = F)
pca1.B <- pca.B$x[,1]
colnames(combat_samplefile)[2] = "Sample_Group"
gene.score <- pca1.B
gene.score <- pca1.A - pca1.B
gene.outtab1 <- data.frame(Sample_Name = combat_samplefile$Sample_Name,pca1.B = pca1.B[combat_samplefile$Sample_Name],#  pca1.A = pca1.A[combat_samplefile$Sample_Name],
                           gene.score = gene.score[combat_samplefile$Sample_Name],
                           gene.group = ifelse(gene.score > 
                                                 -0.8453973, "High", "Low"),#1239行 roc_result$threshold  -0.8453973  ,median(gene.outtab$gene.score)
                           subtype = combat_samplefile$Sample_Group,
                           stringsAsFactors = F)


gene.outtab1$subtype = factor(gene.outtab1$subtype, levels = c("response","non_response"))

ggplot(data = gene.outtab1,aes(x = subtype, y = gene.score, fill = subtype))+
  scale_fill_manual(values = c("#ffdd00", "#00a4e4", "#D14039")) +
  geom_violin(alpha=0.4, position = position_dodge(width = .75),
              size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = TRUE, outlier.size = -1,
               color="black", lwd=0.8, alpha = 0.7)+ # 背景色透明化
  geom_point(shape = 21, size=2,
             position = position_jitterdodge(),
             color="black", alpha=1)+ # 边框线黑色
  theme_classic() +
  ylab(expression("gene score")) +
  xlab("Subtype")  +
  theme(axis.ticks = element_line(size=0.2,color="black"),
        axis.ticks.length = unit(0.2,"cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  stat_compare_means(method = "wilcox.test", label.y = max(gene.score)+1) #kruskal.test/wilcox.test

#堆叠------------------------------------------------------
data = gene.outtab1[gene.outtab1$Sample_Name%in% combat_samplefile$Sample_Name,]
colnames(data)[1] = "Sample_Name"
data = merge(data, combat_samplefile, by = "Sample_Name")

ggboxplot(data, x="subtype", y="gene.score", color="subtype",
          xlab="subtype",
          ylab="gene.score",
          #legend.title=x,
          palette = c("orange","green"),#,"green"
          add = "jitter")+
  stat_compare_means(comparisons=list(c("response","non_response")), symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")#,c("sub1","sub3"),c("sub2","sub3")


rt = data[,c(4,5)]
colnames(rt)=c("gene.group", "subtype")
p = ggplot(rt, aes(gene.group)) +
  geom_bar(aes(fill=subtype), position="fill") +
  labs(x = 'gene.group', y = '', title="Methotrexate") + guides(fill = guide_legend(title ="")) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))
p

#构建诊断模型--------------------------------------
rm(list=ls())
options(stringsAsFactors = F)
setwd("F:/GSE/PP/GSE13355.1.5/")
Intersect <- read.table("Intersect.txt", header=TRUE, row.names=NULL, quote="", sep="\t", check.names = FALSE)
load("step1_combat.rda")
data=t(combat_exprSet[Intersect$Intersect,])
data1=as.data.frame(data)
data1$Sample_Name = rownames(data1)

samplefile <- read.table("combat_samplefile.txt", header=TRUE, row.names=NULL, quote="", sep="\t", check.names = FALSE)

data2 = merge(data1, samplefile, by="Sample_Name")
rownames(data2)= data2$Sample_Name
data2=data2[,-1]
#write.table(data2, file = "data2.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
#data2 <- read.table("data2.txt", header=TRUE, row.names=NULL, quote="", sep="\t", check.names = FALSE)
#data3=data2[,c(10,1,2,3,4,5,6,7,8,9)]
data3=data2[,c(37,1:36)] #最后一列，1：最后的基因列
colnames(data3)[1] = "group"



library(caret)
set.seed(1234)
train_dat = data3
inTrain = createDataPartition(data3$group, p = 5/10, list = FALSE)
train_dat = data3[inTrain, ]
test_dat = data3[-inTrain, ]

#支持向量机--------------------------------------------------
#set.seed(1234)
#library(e1071)
#svm.model <- svm(group ~ ., data=train_dat, type="C-classification", kernel="linear", probability=TRUE, scale=FALSE)
#print(svm_model)

set.seed(1234)
train_ctrl = trainControl(method="repeatedcv", number=10, repeats=10, classProbs=TRUE, summaryFunction=twoClassSummary, allowParallel=TRUE)
svm.model = train(group ~ ., data=train_dat, method="svmLinear", metric="ROC", trControl=train_ctrl)
print(svm.model)
#svm.model = train(group ~ ., data=train_dat, method="svmRadial", metric="ROC", trControl=train_ctrl)

#Accuracy
pred1 <- predict(svm.model, newdata=train_dat[, -1], by="class")
table(pred1, train_dat$group)
mean(pred1 == train_dat$group)
#0.787234

#Accuracy
pred1 <- predict(svm.model, newdata=test_dat[, -1], by="class")
table(pred1, test_dat$group)
mean(pred1 == test_dat$group)
#0.775

save(train_dat,test_dat,file = "ml_dat.rda")
save(svm.model, file = "svm.rda")

#随机森林---------------------------------------------------
#set.seed(1234)
#library(randomForest)
#ntree <- 1000
#rf.model <- randomForest(group ~ ., data=train_dat, ntree=1000, proximity=TRUE, importance=TRUE)
#print(rf.model)

set.seed(1234)
library(caret)
train_ctrl = trainControl(method="repeatedcv", number=10, repeats=10, classProbs=TRUE, summaryFunction=twoClassSummary, allowParallel=TRUE)
rf.model = train(group ~ ., data=train_dat, method="rf", metric="ROC", trControl=train_ctrl)
print(rf.model)

#Accuracy
pred1 <- predict(rf.model, newdata=train_dat[, -1], by="class")
table(pred1, train_dat$group)
mean(pred1 == train_dat$group)
#1
#Accuracy
pred1 <- predict(rf.model, newdata=test_dat[, -1], by="class")
table(pred1, test_dat$group)
mean(pred1 == test_dat$group)
#0.775

#0.7435065
save(rf.model, file = "rf.model.rda")



#逻辑回归---------------------------------------------------
#library(glmnet)
#set.seed(1234)
#glm.model <- glmnet(as.matrix(train_dat[, -1]), train_dat$group, family = "binomial")
#print(glm.model)


set.seed(1234)
train_ctrl = trainControl(method="repeatedcv", number=10, repeats=10, classProbs=TRUE, summaryFunction=twoClassSummary, allowParallel=TRUE)
glm.model = train(group ~ ., data=train_dat, method="glm", metric="ROC", trControl=train_ctrl)
print(glm.model)

#Accuracy
pred1 <- predict(glm.model, train_dat[, -1], by="response")
table(pred1, train_dat$group)
mean(pred1 == train_dat$group)
#0.7765957

#Accuracy
pred1 <- predict(glm.model, test_dat[, -1], by="response")
table(pred1, test_dat$group)
mean(pred1 == test_dat$group)
#0.775

pred2 <- predict(glm.model, newdata=test_dat[, -1], by="response")
table(pred2, test_dat$group)
mean(pred2 == test_dat$group)
#0.6850649
save(glm.model, file = "glm.model.rda")

#ROC曲线---------------------------------------------------
library(pROC)
svmL.probs = predict(svm.model, test_dat[, -1], type="prob")
svmL.ROC = roc(test_dat$group, svmL.probs[, 2], percent=TRUE)
svmL.ROC$auc


rf.probs = predict(rf.model, test_dat[, -1], type="prob")
rf.ROC = roc(test_dat$group, rf.probs[, 2], percent=TRUE)
rf.ROC$auc



glm.probs = predict(glm.model, test_dat[, -1], type="prob")
glm.ROC = roc(test_dat$group, glm.probs[, 2], percent=TRUE)
glm.ROC$auc


plot.roc(test_dat$group, rf.probs[, 2], percent=TRUE, smooth=T, col="black")
plot.roc(test_dat$group, svmL.probs[, 2], percent=TRUE, smooth=T, col="#c1d82f", add=TRUE)
plot.roc(test_dat$group, glm.probs[, 2], percent=TRUE, smooth=F, col="#00a4e4", add=TRUE)

legend("bottomright", legend=c("rf.model AUC=99.91%", "svm.model AUC=99.72%", "glm.model AUC=99.83%"), 
       col=c("black", "#c1d82f", "#00a4e4", "#840000"), lwd=2, cex=0.5, text.width=80)
dev.off()


library(pROC)
svmL.probs = predict(svm.model, train_dat[, -1], type="prob")
svmL.ROC = roc(train_dat$group, svmL.probs[, 2], percent=TRUE)
svmL.ROC$auc


rf.probs = predict(rf.model, train_dat[, -1], type="prob")
rf.ROC = roc(train_dat$group, rf.probs[, 2], percent=TRUE)
rf.ROC$auc


glm.probs = predict(glm.model, train_dat[, -1], type="prob")
glm.ROC = roc(train_dat$group, glm.probs[, 2], percent=TRUE)
glm.ROC$auc


plot.roc(train_dat$group, rf.probs[, 2], percent=TRUE, smooth=F, col="black")
plot.roc(train_dat$group, svmL.probs[, 2], percent=TRUE, smooth=F, col="#c1d82f", add=TRUE)
plot.roc(train_dat$group, glm.probs[, 2], percent=TRUE, smooth=F, col="#00a4e4", add=TRUE)

legend("bottomright", legend=c("rf.model AUC=100%", "svm.model AUC=100%", "glm.model AUC=100%"), 
       col=c("black", "#c1d82f", "#00a4e4", "#840000"), lwd=2, cex=0.5, text.width=80)
dev.off()


