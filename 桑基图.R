#桑基图：
rm(list=ls())
options(stringsAsFactors = F)
setwd("D:\\Rstudio\\test\\3.glaucoma\\Results\\Drug") 
#1.R包的引用
install.packages("ggalluvial")
library(ggalluvial)
library(ggplot2)
library(dplyr)
#对于桑基图的绘制，主要通过ggalluvial包来实现，包括geom_stratum()函数和geom_flow()函数。
#2.数据的准备与读取

#首先，是数据的准备情况；在前期的分析过程中，我们应该计算两组变量之间的相互性，并整理成表格。
#使用R读取数据；
#setwd()
rt<-read.table("dgidb_export(11 genes)_toop30.txt",header=TRUE, row.names=NULL, quote="", sep="\t", check.names = FALSE)
head(rt)
rt<-rt[,c(4,5,10)]
#对数据的流向进行定义
San<- to_lodes_form(rt[,c(1,2)],axes=2:1,id= "Group")
#对于图形的from和to分别进行了定义，即从第2列流向第1列；当修改为axes = 2:1时，则表示图形从第1列流向第2列。
#在绘制过程中，其主要根据抵达数据“to”中不同类别的总和，按降序排列后再行绘制。

#3.桑基图的绘制

ggplot(San,   #定义流向后的数据集
       aes(x = x, stratum = stratum, alluvium = Group,fill = stratum, label = stratum)) + 
  geom_flow(width = 1/8,aes.flow = "forward") +  
  #geom_flow()函数控制边的视觉通道映射设定，也就是线条的颜色；主要由alluvium和weight决定
  #forward表示线条的颜色和前面的一致，backward表示和后面一致
  geom_stratum(alpha=0.9,width=1/10)+
  #geom_stratum()控制节点的视觉通道映射设置，主要由stratum和weight决定
  geom_text(stat = "stratum",size=3.5,color="black")+
  #size=3.5代表名字大小
  scale_x_discrete(expand = c(0,0))+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())+   #去掉坐标轴
  theme(panel.grid = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank())+  #显示分组名字
  ggtitle("")+
  guides(fill=FALSE)

dev.off()
    
                             
                 