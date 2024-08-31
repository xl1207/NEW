#单个基因和免疫细胞的相关性棒棒糖图----------------
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(dplyr)


# 打开两个想进行相关性分析的清洁数据
# 数据1为行名为样品名，列名为基因名
# 数据2为行名为样品名，列名为想做的和基因有相关性关系的，如免疫细胞

# 样品1:expr_data
# 样品2:immu_data
#先在训练集中作"FAM3C"   "BTLA"    "STRBP"   "RASGRP3" "COBLL1"  "CD79A"   "RRAS2"   "MS4A1"   "CXCR5"   "TCL1A" 

exprSet3 = exp       #表达矩阵
exprSet3 = t(exprSet3)

#进行相关性计算
gene <- "ITIH3"     #ITIH3
y <- as.numeric(exprSet3[,gene])

res_xcell2 = as.data.frame(t(res_xcell1))  #xcell结果：列名为样本
cor_data <- do.call(rbind,lapply(colnames(res_xcell2),function(x){
  dd <- cor.test(as.numeric(res_xcell2[,x]),y,method ="spearman",exact=FALSE)
  data.frame(cell=x,cor=dd$estimate,p.value=dd$p.value)
}))


# 画单基因与样本2中所有列名的相关性的棒棒糖图

#定义圆圈颜色的函数
p.col = c('gold','pink','orange','LimeGreen','darkgreen')
fcolor = function(x,p.col){
  color = ifelse(x>0.8,p.col[1],ifelse(x>0.6,p.col[2],ifelse(x>0.4,p.col[3],
                                                             ifelse(x>0.2,p.col[4], p.col[5])
  )))
  return(color)
}

#定义设置圆圈大小的函数
p.cex = seq(2.5, 5.5, length=5)
fcex = function(x){
  x=abs(x)
  cex = ifelse(x<0.1,p.cex[1],ifelse(x<0.2,p.cex[2],ifelse(x<0.3,p.cex[3],
                                                           ifelse(x<0.4,p.cex[4],p.cex[5]))))
  return(cex)
}

dat <- cor_data
gene <- 'ITIH3'
#根据pvalue定义圆圈的颜色
points.color = fcolor(x=dat$p.value,p.col=p.col)
dat$points.color = points.color

#根据相关系数定义圆圈的大小
points.cex = fcex(x=dat$cor)-0.8
dat$points.cex = points.cex
dat=dat[order(dat$cor),]


xlim = ceiling(max(abs(dat$cor))*10)/10         #x轴范围
pdf(file="GSE66229_ITIH3棒棒糖.pdf", width=12, height=12)      #输出图形
layout(mat=matrix(c(1,1,1,1,1,0,2,0,3,0),nc=2),width=c(8,2.2),heights=c(1,2,1,2,1))
par(bg="white",las=1,mar=c(5,18,2,4),cex.axis=1.5,cex.lab=2)
plot(1,type="n",xlim=c(-xlim,xlim),ylim=c(0.5,nrow(dat)+0.5),xlab="Correlation Coefficient",ylab="",yaxt="n",yaxs="i",axes=F)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col="#F5F5F5",border="#F5F5F5")
grid(ny=nrow(dat),col="white",lty=1,lwd=2)
#绘制图形的线段
segments(x0=dat$cor,y0=1:nrow(dat),x1=0,y1=1:nrow(dat),lwd=4)
#绘制图形的圆圈
points(x=dat$cor,y = 1:nrow(dat),col = dat$points.color,pch=16,cex=dat$points.cex)+
  scale_size_continuous(range =c(2,4))
#展示免疫细胞的名称
text(par('usr')[1],1:nrow(dat),dat$cell,adj=1,xpd=T,cex=1.5)
#展示pvalue
pvalue.text=ifelse(dat$p.value<0.001,'<0.001',sprintf("%.03f",dat$p.value))
redcutoff_cor=0
redcutoff_pvalue=0.05
text(par('usr')[2],1:nrow(dat),pvalue.text,adj=0,xpd=T,col=ifelse(abs(dat$cor)>redcutoff_cor & dat$p.value<redcutoff_pvalue,"red","black"),cex=1.5)
axis(1,tick=F)

#绘制圆圈大小的图例
par(mar=c(0,4,3,4))
plot(1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=c(0.1,0.2,0.3,0.4,0.5),col="black",pt.cex=p.cex,pch=16,bty="n",cex=2,title="abs(cor)")

#绘制圆圈颜色的图例
par(mar=c(0,6,4,6),cex.axis=1.5,cex.main=2)
barplot(rep(1,5),horiz=T,space=0,border=NA,col=p.col,xaxt="n",yaxt="n",xlab="",ylab="",main="pvalue")
axis(4,at=0:5,c(1,0.8,0.6,0.4,0.2,0),tick=F)
dev.off()

