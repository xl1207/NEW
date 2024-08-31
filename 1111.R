rm(list = ls())
install.packages("ggord")
BiocManager::install("data.table")
BiocManager::install("classInt")

install.packages("pROC")
install.packages("nomogramFormula")
install.packages("DynNom")


library(ggpubr)
library(caret)
library(randomForest)

2 5 99 100 264 266 324 380
29logFC 122cutHeight 173软阈值 203模块至少基因数 222模块合并0.25 275309颜色 526GiNi系数
基因552 566 581 591 619 633 657 667
路径539 603
免疫检查点 891 893 905


library(org.Hs.eg.db)
g_id = mapIds(x = org.Hs.eg.db,#注释包
              keys = BIRC5, #需要转换的基因Symbol
              keytype = "SYMBOL", #需要转换的类型
              column = "ENTREZID") #需要转换为的类型
head(g_id)


exprSet <- read.table("GBM.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
SampleFile <- read.table("GBM_SampleFile.txt", header = TRUE, sep="\t", row.names = NULL, check.names = FALSE)#comment.char = "!"
exprSet <- new_exprSet
