rm(list = ls())
options(stringsAsFactors = F)
setwd("E:/share/GPL/")

#Step1 -idmap------------------------------------------------------------
install_local("E:/share/soft/idmap2-master/")
library(devtools)
library(idmap1)
ls('package:idmap1')
library(idmap2)
ls('package:idmap2')
library(idmap3)
ls('package:idmap3')
gpl = idmap1::getIDs('GPL570')
gpl = idmap2::get_soft_IDs('GPL6244')
gpl = idmap3::get_pipe_IDs('GPL23080')


write.table(gpl, file = "GPL10558.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

#第二种方式-------------------------------
gpl <- read.table("GPL10558-50081.txt", header=TRUE, row.names=NULL, quote="", sep="\t", check.names = FALSE)


library(GEOquery)
library(Biobase)
gpl <- getGEO(filename = "GPL10558-50081.txt")
colnames(Table(gpl))

gpl <- read.delim("GPL10558-50081.txt", header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, comment.char = "#")
colnames(gpl)
#gpl = gpl[gpl$gene_assignment != "---" & gpl$gene_assignment != "", ]
#library(stringr)
#gpl$gene_symbol = str_split(gpl$gene_assignment, " // ", simplify = TRUE)[ , 2]

#--------------------------
probe2symbol = gpl[ , c(1,2)]
colnames(probe2symbol) = c("probe_id", "gene_symbol")
dim(probe2symbol)
fil <- is.na(probe2symbol$gene_symbol) | probe2symbol$gene_symbol == "" 
probe2symbol <- probe2symbol[!fil, ]
dim(probe2symbol)
probe2symbol2 <- probe2symbol[-grep("///", probe2symbol$gene_symbol), ]
#dim(probe2symbol2)
write.table(probe2symbol, file = "GPL13607_matrix.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

#hgnc---------------------------------------------------------------------------
hgnc <- read.table("E:/share/genome/hgnc_mRNA_matrix.txt", header=TRUE, row.names=NULL, quote="", sep="\t", check.names = FALSE)
colnames(hgnc) = c("gene_name", "gene_symbol")
probe2symbol1 <- merge(probe2symbol, hgnc, by="gene_symbol", all.x=TRUE)
probe2symbol1 = na.omit(probe2symbol1) 
probe2symbol1 = probe2symbol1[, -1]
dim(probe2symbol1)
write.table(probe2symbol1, file = "GPL6244_matrix.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

