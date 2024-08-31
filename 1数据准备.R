rm(list = ls())
options(stringsAsFactors = F)
setwd("F:/GPL")

#Step1 -idmap------------------------------------------------------------
devtools::install_local("F:/Rpackages/ggord-master/")
#remotes::install_github("GuangchuangYu/limma", build = FALSE,force = TRUE)
#devtools::install_local("F:/Rpackages/idmap2-master/")
#devtools::install_local("F:/Rpackages/ggrepel-master/")
library(idmap1)
library(idmap2)
library(idmap3)
ls('package:idmap1')
ls('package:idmap2')
ls('package:idmap3')
gpl = idmap1::getIDs('GPL571')
gpl = idmap2::get_soft_IDs('GPL6480')
gpl = idmap3::get_pipe_IDs('GPL19612')
#project <- (files = "GPL570.txt" expr = "SAMPLE", probeid = "ID_REF", other.columns = "Detection Pval")
#读入手动加工的GPL
gpl <- read.table("GPL19571.txt", header = TRUE, sep = "\t", row.names = NULL, check.names = FALSE)
#写出idmap中的GPL
library(data.table)
write.table(gpl, file = "GPL6480.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
gpl = fread("GPL6480.txt", header = T, data.table = F)


#probe2symbol = Table(gpl)[,c(1,3)]

probe2symbol = gpl[ , c(1, 2)]
colnames(probe2symbol) = c("probe_id", "gene_symbol")
dim(probe2symbol)
fil <- is.na(probe2symbol$gene_symbol) | probe2symbol$gene_symbol == "" 
probe2symbol <- probe2symbol[!fil, ]
dim(probe2symbol)
#若下两句结果为0，则第三句改为probe2symbol，hgnc merge时也是
probe2symbol2 <- probe2symbol[-grep("///", probe2symbol$gene_symbol), ]
dim(probe2symbol2)
#不是0，下一句与hgnc merge句使用probe2symbol2
write.table(probe2symbol, file = "GPL19571_matrix.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

#hgnc---------------------------------------------------------------------------
hgnc <- read.table("F:/GPL/hgnc_mRNA_matrix.txt", header=TRUE, row.names=NULL, quote="", sep="\t", check.names = FALSE)
colnames(hgnc) = c("gene_name", "gene_symbol")
probe2symbol1 <- merge(probe2symbol, hgnc, by="gene_symbol", all.x=TRUE)
probe2symbol1 = na.omit(probe2symbol1) 
probe2symbol1 = probe2symbol1[, -1]
dim(probe2symbol1)
write.table(probe2symbol1, file = "GPL19571_matrix.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


View(hgnc)

#Step2.1-oligo处理affymetix数据-----------------------------------------
rm(list = ls())
options(stringsAsFactors = F)
setwd("F:/GSE/YWXY/GSE69967药物响应")

#
library(oligo)
#BiocManager::install("oligo")
celFiles <- list.files(path = "GSE69967_RAW", pattern = ".CEL.gz", full.names = TRUE)
data.raw <- read.celfiles(filenames = celFiles)
data.eset <- oligo::rma(data.raw)
#
#探针过滤(P/A过滤)
calls <- paCalls(data.raw)
class(calls)
#---
#class(calls)为list
pmas <- calls$calls
absent <- rowSums(pmas == "A")
absent <- which(absent == ncol(pmas))
data.filter <- data.eset[-absent, ]
nrow(data.eset); nrow(data.filter)
data.exprs <- exprs(data.filter)
#---
#class(calls)为matrix
P <- apply(calls, 1, function(x) any(x < 0.05))
pids <- as.numeric(names(P[P]))
pinfo <- getProbeInfo(data.raw)
fids <- pinfo[pinfo$fid %in% pids, 2]
data.filter <- data.eset[rownames(data.eset) %in% fids, ]
nrow(data.eset); nrow(data.filter)
data.exprs <- exprs(data.filter)

#install.packages("tidyverse")
library("tidyverse")

#install.packages("stringi")
library("stringi")
colnames(data.exprs) = str_split(colnames(data.exprs), "\\_", simplify = TRUE)[ , 1]#\\_或\\.
data.exprs <- as.data.frame(data.exprs)
data.exprs <- cbind(rownames(data.exprs), data.exprs)
colnames(data.exprs)[1] <- "probe_id"
write.table(data.exprs, file = "GSE69967_series_matrix.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


###确认gse是否FPKM格式
h <- read.table("C:/R/GSE/GSE41613/GSE41613_gene_matrix.txt", header=TRUE, row.names=NULL, quote="", sep="\t", check.names = FALSE)
h=h[,-1]
range(h)
colSums(h)

#Step2.2-用limma处理illumina数据(non-normalized)---------------------------------
#colnames()[seq(1, 10, by = 6)]=paste0(,A$)

library(limma)
library(oligo)
setwd("C:/R/GSE/GSE65858")
#背景矫正
project <- read.ilmn(files = "GSE65858_non-normalized.txt", expr = "SAMPLE", probeid = "ID_REF", other.columns = "Detection Pval")
project.neqc <- neqc(project, detection.p="Detection Pval", offset = 16)

#探针过滤
isexpr <- rowSums(project.neqc$other$`Detection Pval` < 0.01) > 0
#table(isexpr)
project.neqc0 <- project.neqc[isexpr, ]
nrow(project.neqc):nrow(project.neqc0)

expdata <- project.neqc0$E

#数据写出
data.exprs<-cbind.data.frame(probe_id=rownames(expdata), expdata)
#expdata <- as.data.frame(expdata)
#expdata <- cbind(rownames(expdata), expdata)
#colnames(expdata)[1] <- "probe_id"
write.table(data.exprs, file = "GSE65858_series_matrix.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

View(expdata)
#Step2.3-aglient------------------------------------------------------------
rm(list = ls())
options(stringsAsFactors = F)
library(marray)
library(limma)
library(oligo)

setwd("F:/GSE/YWXY/GSE72246药物响应")
SDRF <- read.table("GSE72246.sdrf.txt", check.names = FALSE, header = TRUE, sep="\t", stringsAsFactors = FALSE)
#names(SDRF)
#SDRF <- readTargets("GSE85716.sdrf.txt", row.names="FileNames", sep="\t")
#FileNames; Array Data File
project <- read.maimages(files = SDRF[ , "Array Data File"], path = "GSE72246_RAW", source = "agilent", green.only = TRUE)
project.bgc <- backgroundCorrect(project, method = "normexp",offset=16)
project.nba <- normalizeBetweenArrays(project.bgc, method = "quantile")

#探针过滤
neg95 <- apply(project.nba$E[project.nba$genes$ControlType == -1, ], 2, function(x) quantile(x, p = 0.5))
cutoff <- matrix(1.1*neg95, nrow(project.nba), ncol(project.nba), byrow = TRUE)
isexpr <- rowSums(project.nba$E > cutoff) >= 4
table(isexpr)
y0 <- project.nba[project.nba$genes$ControlType == 0 & isexpr, ]

yavedata <- avereps(y0, ID = y0$genes$ProbeName)#y0$genes[,"SystematicName"]
colnames(yavedata$E) <- SDRF$FileNames

expdata <- as.data.frame(yavedata$E)

expdata <- cbind(rownames(expdata), expdata)
colnames(expdata)[1] <- "probe_id"
write.table(expdata, file = "GSE72246_series_matrix.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)


A = A[,-1]
