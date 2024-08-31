setwd("/Users/wangxiaoqi/Desktop/MR/mTOR/")
#1---------------------------------------------------------

dir = list.files()
library(data.table)
merge.data=fread("HIF1A.13089.6.3_chrom_22_meta_final_v1.tsv", header = T, data.table=F)

for( i in (2:length(dir))){
  print(i)
  new.data = fread(dir[i], header = T, data.table=F)
  merge.data = rbind(merge.data,new.data)
} 


merge.data$MR = paste(merge.data$chromosome, merge.data$position, sep="_")
merge.data$P <- 10^(merge.data$`log(P)`)
save(merge.data,file="HIF1A_merge.data.rda")
View(head(merge.data, 20))

load(file = "HIF1A_merge.data.rda")
load(file="DY_EIF4EBP2.rda")
merge.data1 = merge.data
merge.data1 <- subset(merge.data1, P<1e-5)

merge.data2 = merge(merge.data1, tmp2, by="MR",all.x = T)
table(duplicated(merge.data2$MR))
setdiff(merge.data1$MR, merge.data2$MR)
write.csv(merge.data2,file="HIF1A.csv")

merge.data2=read.csv("HIF1A.csv",header=T,sep=",",row.name=1)


colnames(merge.data2)
merge.data2 = merge.data2[, c(13,6,5,7,8,10)]
colnames(merge.data2) = c("SNP","effect_allele","other_allele","beta","se","p")
merge.data2$effect_allele = toupper(merge.data2$effect_allele)
merge.data2$other_allele = toupper(merge.data2$other_allele)
save(merge.data2,file="HIF1A_merge.data2_1e-5.rda")
write.csv(merge.data2,file="HIF1A.csv")



#2-----------------------------------------------
install.packages("remotes") 
remotes::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)


setwd("E:/MR/mTOR/data")
library(TwoSampleMR)
#第一步  读取exposure数据
mTOR<-read.table("EIF4G.csv",header = T,sep = ",",row.name=1)

##1查看暴露数据去掉eaf/maf值小于0.01的
mTOR$eaf.exposure
#mTOR = mTOR[-26,]#例如去除第26个[不对？]

#2计算F值与R2
mTOR$r2=2*(1-mTOR$Minor.Allele.Global.Frequency)*(mTOR$Minor.Allele.Global.Frequency)*(mTOR$Minor.Allele.Global.Frequency)^2#公式见PPT
mTOR$F=(mTOR$r2/(1-mTOR$r2))*(mTOR$samplesize.exposure-2)
mTOR=mTOR[mTOR$F>=10,]#确保F值大于10
#3计算总R2和总F值(总F>100)
  sum(mTOR$R2)#直接求和
F_all = (sum(mTOR$R2)/(1-sum(mTOR$R2)))*((unique(mTOR$samplesize.exposure)-length(mTOR$samplesize.exposure)-1)/length(mTOR$samplesize.exposure))
F_all

#4换列名

mTOR1 <- format_data(dat=mTOR,type = "exposure", beta_col = "beta",se_col = "se",
                           effect_allele_col = "effect_allele",other_allele_col = "other_allele",
                           pval_col = "p")
#View(mTOR1)


#5相关性设置，并将文件放到TwoSampleMR包所在文件夹
mTOR2<-subset(mTOR1,pval.exposure<5e-6)  #5e-8/1e-5

#6去除连锁不平衡
mTOR3<-clump_data(mTOR3,clump_r2=0.001, clump_kb = 10000, pop = "EUR") #r2可以调整
View(mTOR3)


##7去除混杂因素--phenoscannerv2
##SNP与混淆因素无关，一般包括年龄、性别、吸烟、饮酒、肥胖
mTOR3$SNP
#write.csv(mTOR3,file = "mTOR3.csv")#保存文件，复制SNP列生成txt，导入网站
#mTOR3 = mTOR3[,-10]#需要去除某SNP时，找到其位置，比如是10，去除

#第二步 从结果把工具变量选出来
t2d_out <- extract_outcome_data(snps=mTOR3$SNP,outcomes='finn-b-I9_VHD',proxies = FALSE, maf_threshold = 0.01, access_token = NULL)
#Harmonize data
dat <- harmonise_data(exposure_dat=mTOR3,outcome_dat=t2d_out,action= 2)
#加OR
generate_odds_ratios(mr_res = mr(dat))#转化成OR二分类
##若为连续型变量，另外计算出来OR值写进文章
##无论暴露是连续性变量还是二分类变量 都统一用beta值进行MR的计算(二分类OR，连续型β)
##所以如果是连续性变量的话 后面还要加上一行代码 把MR里面的OR值算出来，写到文章里面就行
#第三步 MR分析
mr_results <- mr(dat,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median"))
mr_results
mr_report(dat)

#mr_scatter_plot(mr_results = mr(dat,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),dat)
#mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
#mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat)) 
#清除图片

#dev.off()

#第四步 敏感性分析
mr_scatter_plot(mr_results = mr(dat,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),dat)
mr_funnel_plot(mr_singlesnp(dat))
##①异质性检测-
mr_heterogeneity(dat)
##②多效性检测-
mr_pleiotropy_test(dat)
##③Leave-one-out analysis
mr_leaveoneout_plot(mr_leaveoneout(dat))
dev.off()
#第五步 MR-PRESSO
#devtools::install_github("rondolab/MR-PRESSO")#安装MR-PRESSO
source("E:/MR/R/run_mr_presso.R")
run_mr_presso(dat,NbDistribution = 10000)

#第六步 根据实际情况
dat$SNP
#找出要去除的SNP，每个人的不一样
#dat1=dat[c(2,8,10)，]#知道是第几列后，直接删


#3-------------------------------------------------------
#设置工作环境
setwd("D:/languageR/MR")

#第一步 导入暴露工具变量
library(TwoSampleMR) #加载R包 
mTOR <-extract_instruments(outcomes= 'ebi-a-GCST005536', p1=5e-8, clump=TRUE, r2=0.001, kb=10000, access_token = NULL) #获取暴露数据

View(mTOR) #查看暴露数据

##1查看暴露数据去掉eaf/maf值小于0.01的
mTOR$eaf.exposure
ifelse(mTOR$eaf.exposure<0.01,'T','F')
#例如去除第26个[不对？]
mTOR = mTOR[-4,]
#2计算R2和F值
mTOR$r2 = 2*(1-mTOR$eaf.exposure)*(mTOR$eaf.exposure)*(mTOR$beta.exposure)^2
mTOR$F = (mTOR$r2/(1-mTOR$r2))*(mTOR$samplesize.exposure-2)
mTOR = mTOR[mTOR$F>=10, ]
View(mTOR)

#3总F

F_all = (sum(mTOR$r2)/(1-sum(mTOR$r2)))*(unique(mTOR$samplesize.exposure)-length(mTOR$samplesize.exposur)-1)/length(mTOR$samplesize.exposur)
save(mTOR,file="ieu-b-35.rda")
load(file = "ieu-b-35.rda")
##6去除混杂因素--phenoscannerv2
##SNP与混淆因素无关，一般包括年龄、性别、吸烟、饮酒、肥胖
mTOR$SNP
write.csv(mTOR,file = "mTOR.csv")#保存文件，复制SNP列生成txt，导入网站
mTOR =mTOR[-c(8,17,18,31,32),]
#需要去除某SNP时，找到其位置，比如是10，去除
#第二步 从结果把工具变量选出来
t2d_out <- extract_outcome_data(snps=mTOR$SNP,outcomes='ieu-a-977',proxies = FALSE, maf_threshold = 0.01, access_token = NULL)
#Harmonize data
dat <- harmonise_data(exposure_dat=mTOR,outcome_dat=t2d_out,action= 3)
#加OR
generate_odds_ratios(mr_res = mr(dat))#转化成OR二分类
##若为连续型变量，另外计算出来OR值写进文章
##无论暴露是连续性变量还是二分类变量 都统一用beta值进行MR的计算(二分类OR，连续型β)
##所以如果是连续性变量的话 后面还要加上一行代码 把MR里面的OR值算出来，写到文章里面就行
#第三步 MR分析
mr_results <- mr(dat,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median"))
mr_results
mr_report(dat)
###library("rmarkdown")
#mr_scatter_plot(mr_results = mr(dat,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),dat)
#mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
#mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat)) 
#清除图片
dev.off()

#第四步 敏感性分析
mr_scatter_plot(mr_results = mr(dat,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),dat)
mr_funnel_plot(mr_singlesnp(dat))
##①异质性检测-
mr_heterogeneity(dat)
##②多效性检测-
mr_pleiotropy_test(dat)
##③Leave-one-out analysis-
mr_leaveoneout_plot(mr_leaveoneout(dat))
dev.off()
#第五步 MR-PRESSO-
devtools::install_github("rondolab/MR-PRESSO")#安装MR-PRESSO
source("D:/languageR/MR/run_mr_presso.R")
run_mr_presso(dat,NbDistribution = 10000)

#第六步 根据实际情况-
dat$SNP
#找出要去除的SNP，每个人的不一样
dat1=dat[-c(2,8,10)，]#知道是第几列后，直接删
mTOR = mTOR[-26,]

