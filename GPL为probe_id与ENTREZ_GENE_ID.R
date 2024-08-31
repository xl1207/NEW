rm(list = ls())
setwd("F:/GSE/PP/GSE79704/")
library(org.Hs.eg.db)
Entrez <- read.table('GPL19983.txt', header=TRUE, row.names=NULL, quote='',fill = TRUE, sep="\t", check.names = FALSE,colClasses = 'character')
k <- keys(org.Hs.eg.db,keytype = 'ENSEMBL') #ENTREZID
Entrez_Symbol <-  select(org.Hs.eg.db,keys = k ,columns = c('ENTREZID','SYMBOL'),keytype = 'ENSEMBL')
#list = select(org.Hs.eg.db,keys = k ,columns = c('ENTREZID','SYMBOL'),keytype = 'ENSEMBL')
Entrez <- Entrez[,c(1,2)] #这里要选择平台文件中 探针列和ENTREZID列
colnames(Entrez)[2] <- 'ENTREZID'
library(dplyr)
probe2symbol <- left_join(Entrez,Entrez_Symbol,by='ENTREZID')
dim(probe2symbol)
probe2symbol <- probe2symbol[,c('ID','SYMBOL')] %>% distinct(ID,.keep_all = T)
Entrez <- na.omit(Entrez)
dim(probe2symbol)
detach("package:dplyr", unload = TRUE)
probe2symbol <- na.omit(probe2symbol)

write.table(probe2symbol[,c(1,2)], file = 'GPL19983.txt', quote = FALSE, sep = '\t',col.names = TRUE, row.names = FALSE)
#write.table(probe2symbol[,c(1,2)], file = 'GPL19983_matrix.txt', quote = FALSE, col.names = TRUE, row.names = FALSE)
#write.table(probe2symbol[,c(1,2)], file = paste(GPL,'.txt',sep=''), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
