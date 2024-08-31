rm(list = ls())
options(stringsAsFactors = F)
setwd("F:/GSE/YWXY/GSE11903药物响应")

exprSet <- read.table("GSE69967_series_matrix.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
SampleFile <- read.table("GSE69967SampleFile.txt", header = TRUE, sep="\t", row.names = NULL, check.names = FALSE)#comment.char = "!"
SampleFile$Sample_Group <- factor(SampleFile$Sample_Group, levels = c("NO", "YES"))
exprSet <- exprSet[ , match(SampleFile$Sample_Name, colnames(exprSet))]

#1.Boxplot 归一化------------------------------------------------------------------------------------------

library(limma)
palette(c("#6DC9F2","#f26f54", "#AABBCC"))
#pdf("Boxplot.pdf", width = 8, height = 6)
boxplot(exprSet, main = "Signal distribution, Not normalized", col = SampleFile$Sample_Group, names = 1:ncol(exprSet), xlab = "Array", ylab = "Intensity", cex = 0.8, las = 2)#outline = FALSE
#legend("topright", levels(SampleFile$Sample_Group), fill = palette(), border = 0, box.lwd = 0, bty = "n", cex = 0.5)
dev.off()

exprSet = normalizeBetweenArrays(as.matrix(exprSet), method = "quantile")

#pdf("Boxplot.pdf", width = 8, height = 6)
boxplot(exprSet, main = "Signal distribution, Quantile normalization",col = SampleFile$Sample_Group, names = 1:ncol(exprSet), xlab = "Array", ylab = "Intensity", cex = 0.8, las = 2)#outline = FALSE
#legend("topright", levels(SampleFile$Sample_Group), fill = palette(), border = 0, box.lwd = 0, bty = "n", cex = 0.5)
dev.off()
#2.pca------------------------------------------------------------------------------------------
library(ggfortify)
exp_data <- as.data.frame(t(exprSet))
exp_data$group <- SampleFile$Sample_Group
#exp_data<-na.omit(exp_data)
data.pca <- prcomp(exp_data[ , -ncol(exp_data)])

#pdf("PCA2.pdf", width = 6, height = 6)
autoplot(data.pca, data = exp_data, label = FALSE, colour = 'group')
dev.off()

#3.使用GPL文件进行探针注释-------------------------------------------------------
ids <- read.table("GPL570_matrix.txt", header = TRUE, sep = "\t", row.names = NULL, check.names = FALSE)
if(T){
  colnames(ids) = c("probe_id", "gene_symbol")#有时gene_name/symbol
  exprSet = exprSet[rownames(exprSet) %in% ids$probe_id, ]
  ids = ids[match(rownames(exprSet), ids$probe_id), ]
  
}


#2.使用最大值进行探针合并
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
write.table(new_exprSet, file = "GSE69967_gene_matrix.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
#存SampleFile与gene_matrix
rownames(new_exprSet) <- new_exprSet[,1]
new_exprSet <- new_exprSet[,-1]  #去掉第一列
exprSet <- new_exprSet
save(exprSet,SampleFile,file="GSE69967.rda")
rm(list = ls())
load("GSE69967.rda")


#FPKM转换log2----
rm(list = ls())
options(stringsAsFactors = F)
setwd("F:/GSE/GSE48350")
h <- read.table("F:/GSE/GSE48350/GSE48350_gene_matrix.txt", header=TRUE, row.names=NULL, quote="", sep="\t", check.names = FALSE)
h=h[,-1]
h <- as.data.frame(log2(h))

x2 <- data.table::fread("GSE48350_gene_matrix.txt",
                        select = c("gene_symbol"))
head(x2)
y=cbind(x2,h)
y<- as.data.frame(y)#data.table格式不能自定行名，因此我们先转换为数据框
rownames(y) <- x2$gene_symbol
write.table(y, file = "GSE48350log_gene_matrix.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
