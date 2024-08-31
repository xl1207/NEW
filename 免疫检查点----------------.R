#免疫检查点----------------
load("step2_km.rda") #分组信息
load("step1_combat.rda") #加载表达数据
cp_gene <- read.table(file = "E:/生信学习资料/免疫浸润/checkpoint.txt", sep="\t", row.names = 1, header = TRUE)#读取基因分类信息

raw_data_cp = as.data.frame(t(combat_exprSet))
raw_data_cp = raw_data_cp[samplefile$Sample_Name,]
#cp_gene = cp_gene[intersect(rownames(cp_gene),rownames(new_exp)),]
#nrDEG_risk = nrDEG_risk[rownames(cp_gene),]
#cp_gene = cbind(cp_gene,nrDEG_risk[,ncol(nrDEG_risk)])
#colnames(cp_gene)[3] = "DEG"

raw_data_cp = as.data.frame(t(raw_data_cp))
#这里做的是两组之间差异的cp
data_cp = raw_data_cp[rownames(raw_data_cp) , ]

#也可以选做常见的几个checkpoints,例如HLA系列，或常见的检查点
#data_cp = new_exp[rownames(new_exp)%in% rownames(cp_gene[cp_gene$Immune_checkpoint == "CD226",]),]
#tmp = c("PD1","PDL1","PDL2")
#tmp = c("CTLA4","CD28","CD86","CD80")
#tmp = c("BTLA","CD160","LAG3","TIM-3")
tmp = c("CD226","CD27","CD28","CD40","CD40LG","CEACAM1","CTLA4","HLA-B","HLA-DPB1","HLA-DRA","ICOS","SIRPA","TIGIT","TNFSF14","HLA-E","BTLA","HLA-A","HLA-DMB","HLA-DOA","HLA-DQA1","CD160","IDO1","BTN3A1","BTN2A2","BTN2A1")
#tmp = c("CCL4","CCL5","CCL19","CCL28","CCR5","CCR6","CCR7","CCR8","CCR9")
#tmp = c("CCL18","CCL21","CCR10")
data_cp = raw_data_cp[rownames(raw_data_cp)%in% tmp,]

data_cp = as.data.frame(t(data_cp))
data_cp$Sample_Name = rownames(data_cp)
#load("step5.rda")# 加载样本分组信息数据
samplefile$Subtype = factor(samplefile$Subtype, levels = c("cluster1","cluster2"))
samplefile$Sample_Name = rownames(samplefile)
data_cp = merge(samplefile[,c("Sample_Name","Subtype")], data_cp, by="Sample_Name")


#write.table(raw_data_cp, file="combat_raw_data.txt", sep="\t", quote=FALSE, row.names=T, col.names=TRUE)

#箱线图
library(ggpubr)
library(tidyr)
long_dat <- pivot_longer(data_cp, cols=3:ncol(data_cp), names_to="Gene", values_to = "Expression")
long_dat$Subtype = factor(long_dat$Subtype)
ggboxplot(long_dat, x="Gene", y="Expression", fill="Subtype",
          xlab="",
          ylab="Gene expression",
          legend.title="Gene",
          width=0.8,
          palette = c("#008ECB", "#EA921D") )+
  rotate_x_text(50)+
  stat_compare_means(aes(group=Subtype),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
dev.off()