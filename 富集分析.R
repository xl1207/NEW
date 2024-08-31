rm(list = ls())
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
setwd("F:/CRC/CRC/TCGA/")
load("ABnrDEG0.58.rda")
#rt=read.table("Deg_log1.txt", header=T, sep="\t", check.names=F)
rt=nrDEG1_all
library("DOSE")
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(GOplot)
#BiocManager::install('GOplot')
pvalueFilter=0.05
qvalueFilter=0.05 

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}


#基因名字转换为基因id
colnames(rt)[1]="Gene"
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO富集分析
kk=enrichGO(gene=gene,OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
write.table(GO,file="ABDEG0.58_GO.txt",sep="\t",quote=F,row.names = F)

#定义显示Term数目
showNum=10
if(nrow(GO)<30){
  showNum=nrow(GO)
}

#柱状图
pdf(file="AB0.58barplot.pdf", width=10, height=12)
bar=barplot(kk, drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

#气泡图
pdf(file="bubble.pdf", width=8, height=12)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

#获取GO信息
go=data.frame(Category=GO$ONTOLOGY, ID=GO$ID, Term=GO$Description, Genes = gsub("/", ", ", GO$geneID), adj_pval = GO$p.adjust)
#读取基因的logFC
genelist <- data.frame(ID=rt$Gene, logFC=rt$logFC)
row.names(genelist)=genelist[,1]
#设置圈图参数
circ <- circle_dat(go, genelist)
termNum =8       #显示GO数目
termNum=ifelse(nrow(go)<termNum,nrow(go),termNum)
geneNum=200      #限定基因数目
geneNum=ifelse(nrow(genelist)<geneNum, nrow(genelist), geneNum)
#绘制圈图
chord <- chord_dat(circ, genelist[1:geneNum,], go$Term[1:termNum])
pdf(file="GOcircos.pdf", width=10, height=10)
GOChord(chord, 
        space = 0.001,           #基因之间的距离
        gene.order = 'logFC',    #????logFCֵ?Ի???????
        gene.space = 0.25,       #????????ԲȦ?????Ծ???
        gene.size = 5,           #????????????С 
        border.size = 0.1,       #??????ϸ
        process.label = 6)       #GO??????С
dev.off()

#聚类图
pdf(file="GOcluster.pdf",width=12, height=10)
GOCluster(circ, 
          go$Term[1:termNum], 
          lfc.space = 0.2,        #logFC????֮???Ŀ?϶??С
          lfc.width = 1,          #logFC??ԲȦ????
          term.space = 0.2,       #logFC??GO֮????϶?Ĵ?С
          term.width = 1)         #GOԲȦ?Ŀ???
dev.off()          

#6.KEGG--------------------------------
R.utils::setOption("clusterProfiler.download.method",'auto')
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, )
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$Gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG, file="AB_KEGG.txt", sep="\t", quote=F, row.names = F)

showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

pdf(file="AB0.58barplotKEGG.pdf", width=10, height=10)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)
dev.off()

pdf(file="bubbleKEGG.pdf", width=8, height=8)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)
dev.off()

kegg=data.frame(Category="ALL", ID = KEGG$ID, Term=KEGG$Description, Genes = gsub("/", ", ", KEGG$geneID), adj_pval = KEGG$p.adjust)
genelist <- data.frame(ID=rt$Gene, logFC=rt$logFC)
row.names(genelist)=genelist[,1]
#????Ȧͼ????
circ <- circle_dat(kegg, genelist)
termNum =5       #??ʾͨ·????Ŀ
termNum=ifelse(nrow(kegg)<termNum,nrow(kegg),termNum)
geneNum=100      #??ʾ????????Ŀ
geneNum=ifelse(nrow(genelist)<geneNum, nrow(genelist), geneNum)
#????ͨ·??Ȧͼ
chord <- chord_dat(circ, genelist[1:geneNum,], kegg$Term[1:termNum])
pdf(file="KEGGcircos.pdf", width=10, height=10)
GOChord(chord, 
        space = 0.001,           #????֮???ļ???
        gene.order = 'logFC',    #????logFCֵ?Ի???????
        gene.space = 0.25,       #????????ԲȦ?????Ծ???
        gene.size = 5,           #????????????С 
        border.size = 0.1,       #??????ϸ
        process.label = 6)       #ͨ·??????С
dev.off()

#ͨ·?ľ???ͼ
pdf(file="KEGGcluster.pdf",width=12, height=10)
GOCluster(circ, 
          kegg$Term[1:termNum], 
          lfc.space = 0.2,        #logFC????֮???Ŀ?϶??С
          lfc.width = 1,          #logFC??ԲȦ????
          term.space = 0.2,       #logFC??ͨ·???Ŀ?϶??С
          term.width = 1)         #ͨ·ԲȦ?Ŀ???
dev.off()
#---------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(rvcheck)
library(ggplot2)
nrDEG_risk=nrDEG1_all
symbol2entrezid <- bitr(nrDEG_risk$gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
nrDEG_risk <- merge(nrDEG_risk, symbol2entrezid, by.x = "gene_symbol", by.y = "SYMBOL")
#write.table(symbol2entrezid, file = "symbol2entrezid.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

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
