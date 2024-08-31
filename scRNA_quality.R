#对单细胞数据进行质控
setwd("F:/BC/GEO/HER2/")
options(stringsAsFactors = F)
#加载R包
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(survival)
library(RcppEigen)
library(rvcheck)
library(ggfun)
library(SingleR)
#.libPaths()
#.libPaths( c( .libPaths(), "/Library/Frameworks/R.framework/Versions/4.2/Resources/library") )


#读取数据
#1、批量读取单细胞的数据-----------------
dir_name=c('N1','N2','N3','N4','N6','T1','T2','T3','T5','T6')
datalist=list()
for (i in 1:length(dir_name)){
  dir.10x = paste0("her2/",dir_name[i])
  my.data <- Read10X(data.dir = dir.10x) 
  datalist[[i]]=CreateSeuratObject(counts = my.data, project = dir_name[i], 
                                   #每个基因至少在3个细胞中表达，每一个细胞至少有250个基因表达
                                   min.cells = 3, min.features = 250)
}
#修改名称
names(datalist)=dir_name
#2、细胞质控####---------------------
# 批量计算线粒体和rRNA占比
for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")# 计算线粒体占比
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")# 计算rRNA占比
  datalist[[i]] <- sce
  rm(sce)
}
#质控前的
violin=list()
for (i in 1:length(datalist)){
  violin[[i]] <- VlnPlot(datalist[[i]],
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                         pt.size = 0.1, 
                         ncol = 4)
}
#将每个图片进行合并
pearplot_befor <- CombinePlots(plots = violin , nrow=length(datalist), legend="none")
pearplot_befor
ggsave(filename = 'QC_before.pdf',plot = pearplot_befor,he=15,wi=15)
#样本合并
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
#统计每一个样本的个数
raw_count <- table(sce@meta.data$orig.ident)
table(sce@meta.data$orig.ident)

pearplot_befor1<-VlnPlot(sce,
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                        pt.size = 0.1, 
                        ncol = 4)
pearplot_befor1
ggsave(filename = 'QC_before1.pdf',plot = pearplot_befor1,he=7,wi=15)
rm(sce)

#过滤
datalist <- lapply(X = datalist, FUN = function(x) {
  x<-subset(x,subset = nFeature_RNA > 500 & #根据结果修改参数
              nFeature_RNA < 6000 & 
              quantile(percent.mt, 0.98) > percent.mt & percent.mt < 10 &
              quantile(percent.Ribo, 0.99) > percent.Ribo & percent.Ribo > quantile(percent.Ribo, 0.01) & 
              nCount_RNA < quantile(nCount_RNA, 0.97) & nCount_RNA > 1000 )
})
#合并数据
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
clean_count <- table(sce@meta.data$orig.ident)
table(sce@meta.data$orig.ident)

#过滤前后样本细胞数据的统计
summary_cells <- as.data.frame(cbind(raw_count,clean_count))
counts <- rbind(as.data.frame(cbind(summary_cells[,1],rep("raw",each = length(summary_cells[,1])))),
                as.data.frame(cbind(summary_cells[,2],rep("clean",each = length(summary_cells[,2])))))
counts$sample <- rep(rownames(summary_cells),times =2)
colnames(counts)<- c("count","Stat","sample")
counts[,1] <- as.numeric(counts[,1])
counts$Stat <- factor(counts$Stat, levels=c("raw", "clean"), ordered=TRUE)
fit_cell_count <- ggplot(data =counts, mapping = aes(x = sample, y=count))+ 
  geom_bar(aes(fill = Stat),stat = 'identity', position = 'dodge') + scale_fill_brewer(palette = "Set1") +
  theme(text=element_text(size=10),legend.title=element_blank(),
        panel.background = element_rect(fill = "white", colour = "black",size = 0.2),
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.background = (element_rect(colour= "white",fill = "white")))

fit_cell_count
ggsave(filename = 'fit_cell_count.pdf',plot = fit_cell_count,width = 9,height = 9)
#质控后的小提琴图
violin_after=list()
for (i in 1:length(datalist)){
  violin_after[[i]] <- VlnPlot(datalist[[i]],
                               features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"), 
                               pt.size = 0.1,
                               ncol = 4)
}
pearplot_after <- CombinePlots(plots = violin_after , nrow=length(datalist), legend="none")
pearplot_after
ggsave(filename = 'QC_after.pdf',plot = pearplot_after,he=15,wi=15)
pearplot_after1 <- VlnPlot(sce,
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"), 
                         pt.size = 0.1,
                         ncol = 4)
pearplot_after1
ggsave(filename = 'QC_after1.pdf',plot = pearplot_after1,he=7,wi=15)
#质控前后图片的合并
pearplot_befor1
pearplot_after1
qc_merge<- CombinePlots(plots = list(pearplot_befor1,pearplot_after1) , 
                        nrow=2, legend='none')
qc_merge
ggsave(filename = 'qc_merge.pdf',plot = qc_merge,he=9,wi=15)


#保存datalist文件
save(datalist,file = 'datalist.RData')

#3、数据预处理####--------------------
#合并数据
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
rm(datalist)
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", 
                            nfeatures = 2000,
                            mean.cutoff=c(0.0125,3),
                            dispersion.cutoff =c(1.5,Inf))
### 可视化前20个高变基因
top20 <- head(VariableFeatures(sce), 20)
plot1 <- VariableFeaturePlot(sce)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE, size=3.0)

feat_20 <- CombinePlots(plots = list(plot1, plot2),legend="bottom")
feat_20
ggsave(filename = 'feat_20.pdf',plot = feat_20,he=10,wi=15)
#ScaleData
scale.genes <-  rownames(sce)
sce <- ScaleData(sce, features = scale.genes)
#样本的分组
meta1<-data.frame(matrix(nrow=length(sce@meta.data$orig.ident), ncol=2)) 
colnames(meta1)=c('Sample','Group1')
meta1$Sample=sce@meta.data$orig.ident
unique(meta1$Sample)
### Group1 Tumor 为原发性肿瘤；Normal：正常
meta1[grep("N1",meta1$Sample),]$Group1="Normal"
meta1[grep("N2",meta1$Sample),]$Group1="Normal"
meta1[grep("N3",meta1$Sample),]$Group1="Normal"
meta1[grep("N4",meta1$Sample),]$Group1="Normal"
meta1[grep("N6",meta1$Sample),]$Group1="Normal"
#meta1[grep("N8",meta1$Sample),]$Group1="Normal"
meta1[grep("T1",meta1$Sample),]$Group1="Tumor"
meta1[grep("T2",meta1$Sample),]$Group1="Tumor"
meta1[grep("T3",meta1$Sample),]$Group1="Tumor"
meta1[grep("T5",meta1$Sample),]$Group1="Tumor"
meta1[grep("T6",meta1$Sample),]$Group1="Tumor"
#meta1[grep("T8",meta1$Sample),]$Group1="Tumor"

sce <- AddMetaData(sce, meta1$Sample,col.name = "Sample")
sce <- AddMetaData(sce, meta1$Group1,col.name = "Group1")
save(sce,file = 'sce.RData')



#4、umap_cluster--------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
#加载sce
load('sce.RData')

#PCA降维，选择合适的拐点
sce <- RunPCA(sce, features = VariableFeatures(sce)) 
dimplot1 <- DimPlot(sce, reduction = "pca") 
elbowplot1 <- ElbowPlot(sce, ndims=50, reduction="pca") 
sc_pca <- dimplot1+elbowplot1
sc_pca
ggsave(filename = 'sc_pca.pdf',plot = sc_pca,he=10,wi=15)
#可视化前2个PC的top20个基因
VizDimLoadings(sce, dims = 1:10, nfeatures = 20, reduction = "pca")

#前20个PC
pdf('1.pdf',he=15,wi=15)
DimHeatmap(sce, dims = 1:20, cells = 500, balanced = TRUE)
dev.off()

#pbmc <- JackStraw(sce,num.replicate = 100)

#pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
#JackStrawPlot(pbmc,dims = 1:20)

#选择锚点和分辨率
Dims <- 50
Resolution <- 0.1
sce <- FindNeighbors(object = sce, dims = 1:Dims)
sce <- FindClusters(object = sce, resolution = Resolution)
#颜色
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#CCCCFF","#000000","#7B68EE","#9400D3","#A0522D","#800080","#D2B48C","#D2691E",
            "#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A",
            "#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513",
            "#DEB887")
length(table(sce@active.ident))
mycolor = allcolour[1:length(table(sce@active.ident))]
#### 按cluster进行占比统计，每个cluster所占的比例
cluster.frequency.table <- sce@meta.data %>%
  dplyr::count(seurat_clusters) %>%
  dplyr::mutate(freq = n / sum(n)*100) %>%
  ungroup()%>%as.data.frame()

cluster.frequency.table
pdf('2.pdf',he = 8,width = 10)
pie(cluster.frequency.table$n, labels=round(cluster.frequency.table$freq,2),radius=1.0, main = "Percentage of Cluster", col=mycolor)   
legend("left",legend=unique(cluster.frequency.table$seurat_clusters),bty="n",fill=mycolor)
dev.off()
#统计每一个分组，每一个亚群所占的比例
cluster.frequency.sample=data.frame()
#sample
for (i in as.character(unique(sce@meta.data$Group1))){
  data1<-sce@meta.data[which(sce@meta.data$Group1==i),]
  dat1 <- data1 %>%
    dplyr::group_by(Group1) %>%
    dplyr::count(seurat_clusters) %>%
    dplyr::mutate(freq = n / sum(n)*100) %>%
    ungroup()%>%as.data.frame()
  cluster.frequency.sample=rbind(cluster.frequency.sample,dat1)
}
head(cluster.frequency.sample)
cluster.freq.sample<-tidyr::spread(data=cluster.frequency.sample[,c("Group1","seurat_clusters","freq")],
                                   key=Group1, value=freq)
cluster.freq.sample[is.na(cluster.freq.sample)]<-0
head(cluster.freq.sample)
#rownames(cluster.freq.sample)=cluster.freq.sample$seurat_clusters
#从内圈到外圈依次是 Normal 、Turmal
cluster.freq<-ggplot(data=cluster.frequency.sample, mapping=aes(x=Group1,y=freq,fill=seurat_clusters))+
  geom_bar(stat='identity',width=0.9)+coord_polar(theta="y",start = 0)+
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank())+
  scale_fill_manual(values=mycolor)
cluster.freq
pdf('cluster_freq.pdf',he=7,wi=9)
cluster.freq
dev.off()
write.csv(cluster.frequency.sample,file ="cluster.frequency.csv")

### UMAP
sce <- RunUMAP(sce, dims=1:Dims, reduction="pca")
###tsne 降维
sce <- RunTSNE(sce, 
               dims=1:Dims, 
               reduction="pca",
               perplexity=30,
               max_iter=1000)

#可视化
sc_umap = DimPlot(sce,cols=mycolor,
                  reduction="umap",
                  #reduction="tsne",
                  label = "T", 
                  pt.size = 0.2,
                  label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        axis.ticks = element_blank()
  ) 
sc_umap
ggsave('sc_umap_cluster.pdf',sc_umap,he=7,wi=7)

sc_tsne = DimPlot(sce,cols=mycolor,
                  #reduction="umap",
                  reduction="tsne",
                  label = "T", 
                  pt.size = 0.2,
                  label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        axis.ticks = element_blank()
  ) 
sc_tsne
ggsave('sc_tsne_cluster.pdf',sc_umap,he=7,wi=7)

# 细胞来自于那个样本umap可能改成tsne
sc_umap_group1 = DimPlot(sce,cols=mycolor,group.by='Sample',
                         reduction="tsne",
                         label = "T", 
                         pt.size = 0.2,
                         label.size = 0) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        axis.ticks = element_blank()
  ) 

sc_umap_group1
ggsave('sc_tsne_sample.pdf',sc_umap_group1,he=7,wi=7)

sc_umap_group2 = DimPlot(sce,cols=mycolor,group.by='Group1',
                         reduction="tsne",
                         label = "T", 
                         pt.size = 0.2,
                         label.size = 0) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        axis.ticks = element_blank()
  ) 

sc_umap_group2
ggsave('sc_tsne_group.pdf',sc_umap_group2,he=7,wi=7)

#5、marker基因的筛选-------------
#需要修改marker基因实际上就是一个差异基因
#dge.cluster <- FindMarkers(sce,ident.1 = 0,ident.2 = 1)
#sig_dge.cluster <- subset(dge.cluster, p_val_adj<0.01&abs(avg_log2FC)>1)
#寻找差异基因时的差异倍数
Logfc = 0.3
#差异基因时最小的表达比例（在某一个亚群里面）
Minpct = 0.25
DefaultAssay(sce) <- "RNA"
sce.markers <- FindAllMarkers(object = sce,logfc.threshold = Logfc, 
                              min.pct = Minpct,only.pos = T)
sce.markers["pct.diff"]=sce.markers$pct.1-sce.markers$pct.2
sce.markers <- sce.markers[sce.markers$p_val_adj<0.05,]
length(unique(sce.markers$gene))
head(sce.markers)

View(sce.markers)
write.table(sce.markers,'scRNA_marker_gene.txt',quote = F,row.names = F,sep='\t')
### 选择前5个marker基因
Top5 <- sce.markers %>% 
  group_by(cluster) %>% 
  slice_max(n =5, order_by = avg_log2FC)  
Top5 <- unique(Top5$gene)

sc_marker_dotplot <- DotPlot(object = sce, 
                             features = Top5,
                             cols=c("blue", "red"),
                             scale = T)+ 
  RotatedAxis()+ ggtitle("Top 5 Marker Genes")+ 
  theme(plot.title = element_text(hjust = 0.5)) 

sc_marker_dotplot
ggsave(filename = 'sc_marker_dotplot.pdf',
       plot = sc_marker_dotplot,
       height = 9,width = 25)
#热图展示
library(viridisLite)
sc_marker_heatmap<- DoHeatmap(object = sce,
                              features = Top5,
                              group.colors = mycolor,
                              label = F) + 
  ggtitle("Top 5 Marker Genes") + 
  theme(plot.title = element_text(hjust = 0.5)) 
sc_marker_heatmap
ggsave(filename = 'sc_marker_heatmap.pdf',
       plot = sc_marker_heatmap,
       width = 12,height = 12)

save(sce,file = 'sce2.RData')

load('sce2.RData')

#setwd("F:/BC/GEO/3/")
#library(devtools)
#devtools::install_local("SingleR.zip")

#6、cluser_anno进行注释---------------------
library(SingleR)
library(celldex)
library(BiocParallel)
library(Seurat)
#BiocManager::install("BiocParallel")
#下载
hpca.se <- HumanPrimaryCellAtlasData()
load('sce2.RData')
#获取基因的表达谱的count数据
testdata <- GetAssayData(sce, slot="data")
#获取聚类的亚群
clusters <- sce@meta.data$seurat_clusters
View(clusters)
pred.sce <- SingleR(test =  testdata, 
                    ref = hpca.se, 
                    labels = hpca.se$label.fine,
                    method = "clusters",
                    clusters = clusters, 
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")

plotScoreHeatmap(pred.sce)

celltype = data.frame(ClusterID=rownames(pred.sce), celltype=pred.sce$labels
                      , stringsAsFactors = F)

#手动注释
# Specify genes
library(patchwork)
#
VlnPlot(object=sce, features = c("PECAM1"), pt.size = 0, assay = "RNA") + NoLegend()
#内皮细胞 #VWF,CLDN5
VlnPlot(object=sce, features = c("CLDN5"), pt.size = 0, assay = "RNA") + NoLegend()
#上皮细胞 "EPCAM"
VlnPlot(object=sce, features = c("EPCAM"), pt.size = 0, assay = "RNA") + NoLegend()
#周细胞 "MCAM" "RGS5"
VlnPlot(object=sce, features = c("MCAM"), pt.size = 0,assay = "RNA") + NoLegend()
#髓系细胞"CD68"
VlnPlot(object=sce, features = c("CD14"), pt.size = 0,assay = "RNA") + NoLegend()
#T细胞 "CD3D","CD3E","CD3G"
VlnPlot(object=sce, features = c("CD3G"), pt.size = 0,assay = "RNA") + NoLegend()
#B细胞 "CD79A"CD19
VlnPlot(object=sce, features = c("CD79A"), pt.size = 0,assay = "RNA") + NoLegend()
#肥大细胞 "MS4A2"
VlnPlot(object=sce, features = c("MS4A2"), pt.size = 0,assay = "RNA") + NoLegend()
#NK细胞 "NKG7"
VlnPlot(object=sce, features = c("NKG7"), pt.size = 0,assay = "RNA") + NoLegend()
#PDc"LI3RA"
VlnPlot(object=sce, features = c("LI3RA"), pt.size = 0,assay = "RNA") + NoLegend()
#成纤维细胞 "COL1A1" "PDGFRA"PDPN
VlnPlot(object=sce, features = c("COL1A1"), pt.size = 0,assay = "RNA") + NoLegend()
#8,18
VlnPlot(object=sce, features = c("FOXP1"), pt.size = 0, assay = "SCT", ncol = 1) + NoLegend()
#左上角
VlnPlot(object=sce, features = c("KRT8","KRT5"), pt.size = 0,assay = "RNA", ncol = 1) + NoLegend()
celltype<-read.delim('celltype.txt',sep='\t',header = T,check.names = F)
celltype = as.data.frame(celltype)
write.table(celltype,'celltype.txt',quote = F,sep = '\t',row.names = F)
#修改亚群的名称
sce <- RenameIdents(object = sce, 
                    "0" = celltype[1,2],
                    "1" = celltype[2,2],
                    "2" = celltype[3,2],
                    "3" = celltype[4,2],
                    "4" = celltype[5,2],
                    "5" = celltype[6,2],
                    "6" = celltype[7,2],
                    "7" = celltype[8,2],
                    "8" = celltype[9,2],
                    "9" = celltype[10,2],
                    "10" = celltype[11,2],
                    "11" = celltype[12,2],
                    "12" = celltype[13,2],
                    "13" = celltype[14,2],
                    "14" = celltype[15,2])

length(table(sce@active.ident))
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#CCCCFF","#000000","#7B68EE","#9400D3","#A0522D","#800080","#D2B48C","#D2691E",
            "#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A",
            "#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513",
            "#DEB887")
length(table(sce@active.ident))
mycolor = allcolour[1:length(table(sce@active.ident))]

umap_celltype<-DimPlot(sce,cols=mycolor,
                       reduction="umap",
                       label = "T", 
                       pt.size = 0.2,
                       label.size = 5)

umap_celltype
ggplot2::ggsave('umap_celltype_new.pdf',plot = umap_celltype,he=7,wi=10)
save(sce,file = 'sce3_new.RData')

#elc手动注释-------------------------------------------------------
# Specify genes
library(patchwork)

#上皮细胞 "EPCAM"
VlnPlot(object=sce, features = c("EPCAM"), pt.size = 0, assay = "RNA") + NoLegend()

#B细胞 "CD79A"CD19
VlnPlot(object=sce, features = c("CD19"), pt.size = 0,assay = "RNA") + NoLegend()

#macrophage#CD163 CD68 CD14
VlnPlot(object=sce, features = c("CD14"), pt.size = 0, assay = "RNA") + NoLegend()

#T细胞 "CD3D","CD3E","CD3G"
VlnPlot(object=sce, features = c("CD3G"), pt.size = 0,assay = "RNA") + NoLegend()

#成纤维细胞 "COL1A1" "PDGFRA"PDPN
VlnPlot(object=sce, features = c("PDGFRB"), pt.size = 0,assay = "RNA") + NoLegend()

#plasma cells CD79A JSRP1
VlnPlot(object=sce, features = c("JSRP1"), pt.size = 0, assay = "RNA") + NoLegend()

#肥大细胞 TPSAB1 CPA3
VlnPlot(object=sce, features = c("CPA3"), pt.size = 0, assay = "RNA") + NoLegend()

#内皮细胞 #PECAM1
VlnPlot(object=sce, features = c("PECAM1"), pt.size = 0, assay = "RNA") + NoLegend()

#NK细胞 "MCAM" "RGS5"
VlnPlot(object=sce, features = c("RGS5"), pt.size = 0,assay = "RNA") + NoLegend()





# rename identities
a = read.table('Celltype_brca1.txt', header=TRUE)
new.cluster.ids <- as.character(a[, 2])
names(new.cluster.ids) <- levels(tnbc_all_cca)
tnbc_all_cca <- RenameIdents(tnbc_all_cca, new.cluster.ids)
#Idents(tnbc_all_cca) <- tnbc_all_cca$seurat_clusters
library(ggsci)
DimPlot(tnbc_all_cca, reduction = "tsne", label = F, pt.size = 0.5, group.by = "ident") + NoLegend()+ 
  theme_classic() + scale_color_aaas() + scale_fill_aaas(alpha=0.7)

tnbc_all_cca$group = str_split(tnbc_all_cca$orig.ident, '[_]', simplify = TRUE)[ , 1]
tnbc_all_cca$batch = str_split(tnbc_all_cca$orig.ident, '[_]', simplify = TRUE)[ , 2]

DimPlot(tnbc_all_cca, reduction = "tsne", label = TRUE, pt.size = 0.7, group.by = "group") + NoLegend()#+ 

#7、亚群的差异分析---------------------------
load('sce3_new.RData')
cell_num<-as.data.frame(table(sce@active.ident,sce$Group1))
colnames(cell_num)=c('clusters','type','cell_num')
cell_num=tidyr::spread(cell_num,type,cell_num)
cell_num
#T/N
cell_num1=data.frame()
for (i in 1:nrow(cell_num)){
  print(i)
  pval=fisher.test(matrix(c(cell_num[i,3],cell_num[i,2],
                            sum(cell_num[,3])-cell_num[i,3],
                            sum(cell_num[,2])-cell_num[i,2]),
                          nrow = 2,ncol = 2)
                   ,alternative = "two.sided")$p.value
  cell_num1[i,'T_celltype']=cell_num[i,3]
  cell_num1[i,'T_no_celltype']=sum(cell_num[,3])-cell_num[i,3]
  cell_num1[i,'N_celltype']=cell_num[i,2]
  cell_num1[i,'N_no_celltype']=sum(cell_num[,2])-cell_num[i,2]
  cell_num1[i,'p.val']=pval
  cell_num1[i,'fc']=(cell_num[i,3]/(sum(cell_num[,3])-cell_num[i,3]))/(cell_num[i,2]/(sum(cell_num[,2])-cell_num[i,2]))
}
cell_num1$cell_name=cell_num$clusters
cell_num1

#差异分析结果的保存
write.table(cell_num1,'cell_type_statistical_new.txt',quote = F,row.names = F,sep='\t')
#筛选差异的亚群
fc=4
sig_cell_type=as.character(cell_num1[which((cell_num1$fc>fc | cell_num1$fc<(1/fc)) & cell_num1$p.val<0.05),'cell_name'])
sig_cell_type
write.table(sig_cell_type,'sig_cell_type_new.txt',quote = F,row.names = F,sep='\t',col.names = F)

#ReactomeGSA 富集分析########
library(ReactomeGSA)
library(Seurat)
library(readr)
DefaultAssay(sce) <- "RNA"

gsva_result <- analyse_sc_clusters(sce, verbose = TRUE)
pathway_expression <- pathways(gsva_result)
colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))
#write.csv(pathway_expression,file = 'reactomegsa_sce.csv',quote = F,row.names = F)

# find the maximum differently expressed pathway
pathway_expression <-  read_csv("reactomegsa_sce_new.csv")
max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
  values <- as.numeric(row[2:length(row)])
  return(data.frame(name = row[1], min = min(values), max = max(values)))
}))

max_difference$diff <- max_difference$max - max_difference$min
max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]


plot_num = 10
plot_gsva <- pathway_expression[rownames(max_difference[1:plot_num,]),]


pdf("max_difference_top20.pdf")
pheatmap::pheatmap(t(plot_gsva[,-1]),
                   scale ="row",
                   angle_col = 90,cellwidth = 15,cellheight = 15,
                   labels_col = plot_gsva[,1],
                   # cellwidth = 10, cellheight = 8,
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

dev.off()




library(tidyverse)
barplot_pathway <- plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[1]) 
barplot_pathway$data %>% mutate(absmy = ifelse(expr>=0, "Z","Fy")) -> barplot_pathway1

pp1<-barplot_pathway1 %>% ggplot(aes(cluster_id,   expr ,fill=absmy))+ 
  geom_bar(stat='identity') + theme_bw()+ xlab('')+ylab('ES')+
  theme(axis.text.x = element_text(angle =90,hjust = .9,size = 10,vjust = 0.9))+ 
  ggtitle(max_difference$name[1]) + theme(legend.position="none")+coord_flip()

pp1
ggsave(filename = 'pp1.pdf',plot = pp1,he=7,wi=7)

barplot_pathway <- plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[2]) 
barplot_pathway$data %>% mutate(absmy = ifelse(expr>=0, "Z","Fy")) -> barplot_pathway1

pp2<-barplot_pathway1 %>% ggplot(aes(cluster_id,   expr ,fill=absmy    ))+ 
  geom_bar(stat='identity') + theme_bw()+ xlab('')+ylab('ES')+
  theme(axis.text.x = element_text(angle =90,hjust = .9,size = 10,vjust = 0.9))+ 
  ggtitle(max_difference$name[2]) + theme(legend.position="none")+coord_flip()
pp2
ggsave(filename = 'pp2.pdf',plot = pp2,he=7,wi=7)



barplot_pathway <- plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[3]) 
barplot_pathway$data %>% mutate(absmy = ifelse(expr>=0, "Z","Fy")) -> barplot_pathway1

pp3<-barplot_pathway1 %>% ggplot(aes(cluster_id,   expr ,fill=absmy    ))+ 
  geom_bar(stat='identity') + theme_bw()+ xlab('')+ylab('ES')+
  theme(axis.text.x = element_text(angle =90,hjust = .9,size = 10,vjust = 0.9))+ 
  ggtitle(max_difference$name[3]) + theme(legend.position="none")+coord_flip()
pp3
ggsave(filename = 'pp3.pdf',plot = pp3,he=7,wi=7)

barplot_pathway <- plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[4]) 
barplot_pathway$data %>% mutate(absmy = ifelse(expr>=0, "Z","Fy")) -> barplot_pathway1

pp4<-barplot_pathway1 %>% ggplot(aes(cluster_id,   expr ,fill=absmy    ))+ 
  geom_bar(stat='identity') + theme_bw()+ xlab('')+ylab('ES')+
  theme(axis.text.x = element_text(angle =90,hjust = .9,size = 10,vjust = 0.9))+ 
  ggtitle(max_difference$name[4]) + theme(legend.position="none")+coord_flip()
pp4
ggsave(filename = 'pp4.pdf',plot = pp4,he=7,wi=7)

barplot_pathway <- plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[5]) 
barplot_pathway$data %>% mutate(absmy = ifelse(expr>=0, "Z","Fy")) -> barplot_pathway1

pp5<-barplot_pathway1 %>% ggplot(aes(cluster_id,   expr ,fill=absmy    ))+ 
  geom_bar(stat='identity') + theme_bw()+ xlab('')+ylab('ES')+
  theme(axis.text.x = element_text(angle =90,hjust = .9,size = 10,vjust = 0.9))+ 
  ggtitle(max_difference$name[5]) + theme(legend.position="none")+coord_flip()
pp5
ggsave(filename = 'pp5.pdf',plot = pp5,he=7,wi=7)

barplot_pathway <- plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[6]) 
barplot_pathway$data %>% mutate(absmy = ifelse(expr>=0, "Z","Fy")) -> barplot_pathway1

pp6<-barplot_pathway1 %>% ggplot(aes(cluster_id,   expr ,fill=absmy    ))+ 
  geom_bar(stat='identity') + theme_bw()+ xlab('')+ylab('ES')+
  theme(axis.text.x = element_text(angle =90,hjust = .9,size = 10,vjust = 0.9))+ 
  ggtitle(max_difference$name[6]) + theme(legend.position="none")+coord_flip()
pp6
ggsave(filename = 'pp6.pdf',plot = pp6,he=7,wi=7)

barplot_pathway <- plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[7]) 
barplot_pathway$data %>% mutate(absmy = ifelse(expr>=0, "Z","Fy")) -> barplot_pathway1

pp7<-barplot_pathway1 %>% ggplot(aes(cluster_id,   expr ,fill=absmy    ))+ 
  geom_bar(stat='identity') + theme_bw()+ xlab('')+ylab('ES')+
  theme(axis.text.x = element_text(angle =90,hjust = .9,size = 10,vjust = 0.9))+ 
  ggtitle(max_difference$name[7]) + theme(legend.position="none")+coord_flip()
pp7
ggsave(filename = 'pp7.pdf',plot = pp7,he=7,wi=7)

barplot_pathway <- plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[8]) 
barplot_pathway$data %>% mutate(absmy = ifelse(expr>=0, "Z","Fy")) -> barplot_pathway1

pp8<-barplot_pathway1 %>% ggplot(aes(cluster_id,   expr ,fill=absmy    ))+ 
  geom_bar(stat='identity') + theme_bw()+ xlab('')+ylab('ES')+
  theme(axis.text.x = element_text(angle =90,hjust = .9,size = 10,vjust = 0.9))+ 
  ggtitle(max_difference$name[8]) + theme(legend.position="none")+coord_flip()
pp8
ggsave(filename = 'pp8.pdf',plot = pp8,he=7,wi=7)

barplot_pathway <- plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[9]) 
barplot_pathway$data %>% mutate(absmy = ifelse(expr>=0, "Z","Fy")) -> barplot_pathway1

pp9<-barplot_pathway1 %>% ggplot(aes(cluster_id,   expr ,fill=absmy    ))+ 
  geom_bar(stat='identity') + theme_bw()+ xlab('')+ylab('ES')+
  theme(axis.text.x = element_text(angle =90,hjust = .9,size = 10,vjust = 0.9))+ 
  ggtitle(max_difference$name[9]) + theme(legend.position="none")+coord_flip()
pp9
ggsave(filename = 'pp9.pdf',plot = pp9,he=7,wi=7)

barplot_pathway <- plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[10]) 
barplot_pathway$data %>% mutate(absmy = ifelse(expr>=0, "Z","Fy")) -> barplot_pathway1

pp10<-barplot_pathway1 %>% ggplot(aes(cluster_id,   expr ,fill=absmy    ))+ 
  geom_bar(stat='identity') + theme_bw()+ xlab('')+ylab('ES')+
  theme(axis.text.x = element_text(angle =90,hjust = .9,size = 10,vjust = 0.9))+ 
  ggtitle(max_difference$name[10]) + theme(legend.position="none")+coord_flip()
pp10
ggsave(filename = 'pp10.pdf',plot = pp10,he=7,wi=7)


#8、细胞轨迹-------
setwd("F:/BC/GEO/3/细胞轨迹/")
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(monocle)
load('sce2.RData')
sig_cell_type<-read.delim('sig_cell_type_new.txt',sep='\t',header = F)
sig_cell_type=as.character(sig_cell_type$V1)
cell_anno<-read.delim('celltype.txt',sep='\t',header = T)
sce.markers<-read.delim('scRNA_marker_gene.txt',sep='\t',header = T)

sce$cell_type=sce$seurat_clusters
for (i in 1:nrow(cell_anno)){
  sce$cell_type=gsub(paste0('^',cell_anno$ClusterID[i],'$'),
                     as.character(cell_anno$celltype[i]),sce$cell_type)
}

#提取差异的亚群的sce文件
sce_sig<-subset(sce,idents = cell_anno[cell_anno$celltype %in% sig_cell_type,1])
save(sce_sig,file = 'sce_sig_new.RData')

rm(sce)

exp.rawdata <- as(as.matrix(GetAssayData(sce_sig,slot = "counts")), 'sparseMatrix')
#构建featuredata，一般featuredata需要两个col，一个是gene_id,一个是gene_short_name,row对应counts的rownames
feature_ann<-data.frame(gene_short_name=rownames(sce_sig))
rownames(feature_ann) <- rownames(exp.rawdata)
scRNA_fd <-new("AnnotatedDataFrame", data = feature_ann)
sample_ann <- sce_sig@meta.data
#rownames(sample_ann)<-colnames(Mono_matrix)
rm(sce_sig)
scRNA_pd<-new("AnnotatedDataFrame", data =sample_ann)
#build new cell data set

scRNA.cds<-newCellDataSet(exp.rawdata,
                          phenoData =scRNA_pd,
                          featureData =scRNA_fd,
                          expressionFamily=negbinomial.size())
rm(exp.rawdata)
#查看phenodata、featuredata
head(pData(scRNA.cds))
head(fData(scRNA.cds))

#计算size factors 和 dispersions，用于后期分析；
scRNA.cds <- estimateSizeFactors(scRNA.cds)
scRNA.cds <- estimateDispersions(scRNA.cds)

# Filtering low-quality cells
scRNA.cds <- detectGenes(scRNA.cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(scRNA.cds),
                                    num_cells_expressed >= 0.01*ncol(scRNA.cds))) 

#筛选基因,这里可以根据自己的需要筛选特定的基因
Top500 <- sce.markers %>% group_by(cluster) %>% 
  slice_max(n =500, order_by = avg_log2FC) %>% 
  slice_max(n =500, order_by =pct.diff)  

ordering_genes <- intersect(unique(Top500$gene),rownames(scRNA.cds))
scRNA.cds <- setOrderingFilter(scRNA.cds, ordering_genes)
plot_ordering_genes(scRNA.cds)

#用DDRtree 进行降维分析
scRNA.cds <- reduceDimension(scRNA.cds,
                             norm_method ='log', 
                             max_components = 2,
                             method = 'DDRTree')

#计算psudotime值
scRNA.cds <- orderCells(scRNA.cds)
head(pData(scRNA.cds))
head(scRNA.cds$cell_type)

save(scRNA.cds,file = 'scRNA.cds.RData')

mono_cell=plot_cell_trajectory(scRNA.cds,cell_size = 1, color_by = "cell_type")+
  theme(legend.position = "right",
        text=element_text(size=10),
        legend.title=element_blank(),
        panel.background = element_rect(fill = "white", colour = "black",size = 0.2), 
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.background = (element_rect(colour= "white",fill = "white")))+
  facet_wrap("~cell_type", ncol=2)+
  guides(shape = guide_legend(override.aes = list(size = 3),nrow = 10),
         color = guide_legend(override.aes = list(size = 3),nrow = 10))
mono_cell
ggsave(filename = 'mono_cell.pdf',mono_cell,he=9,wi=9)

mono_state<-plot_cell_trajectory(scRNA.cds,cell_size = 1,
                                 color_by = "State")+
  theme(legend.position = "right",
        text=element_text(size=10),
        panel.background = element_rect(fill = "white", 
                                        colour = "black",size = 0.2), 
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.background = (element_rect(colour= "white",fill = "white")))
mono_state
ggsave(filename = 'mono_state.pdf',mono_state,he=5,wi=5)

mono_time=plot_cell_trajectory(scRNA.cds,cell_size = 1, color_by = "Pseudotime")+
  scale_color_gradient(low = 'red',high = 'blue')+
  theme(legend.position = "right",
        text=element_text(size=10),
        panel.background = element_rect(fill = "white", colour = "black",size = 0.2), 
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.background = (element_rect(colour= "white",fill = "white")))
mono_time
ggsave(filename = 'mono_time.pdf',mono_time,he=5,wi=5)

#四种亚群，每一个做轨迹
sig_cell_type[1]
cell_anno[cell_anno$celltype %in% sig_cell_type[1],1]
sig_cell_type[2]
cell_anno[cell_anno$celltype %in% sig_cell_type[2],1]
sig_cell_type[3]
cell_anno[cell_anno$celltype %in% sig_cell_type[3],1]
sig_cell_type[4]
cell_anno[cell_anno$celltype %in% sig_cell_type[4],1]

mono_celltype=plot_cell_trajectory(scRNA.cds,cell_size = 1, color_by = "seurat_clusters")+
  theme(legend.position = "right",
        text=element_text(size=10),
        legend.title=element_blank(),
        panel.background = element_rect(fill = "white", colour = "black",size = 0.2), 
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.background = (element_rect(colour= "white",fill = "white")))+
  facet_wrap("~cell_type", ncol=2)+
  guides(shape = guide_legend(override.aes = list(size = 3),ncol = 2),
         color = guide_legend(override.aes = list(size = 3),ncol = 2))
mono_celltype
ggsave(filename = 'mono_celltype.pdf',plot = mono_celltype,width = 9,height = 5)

mono_cluster=plot_cell_trajectory(scRNA.cds,cell_size = 1, color_by = "seurat_clusters")+
  theme(legend.position = "right",
        text=element_text(size=10),
        legend.title=element_blank(),
        panel.background = element_rect(fill = "white", colour = "black",size = 0.2), 
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.background = (element_rect(colour= "white",fill = "white")))+
  facet_wrap("~seurat_clusters", ncol=3)+
  guides(shape = guide_legend(override.aes = list(size = 3),ncol = 2),
         color = guide_legend(override.aes = list(size = 3),ncol = 2))
mono_cluster
ggsave(filename = 'mono_cluster.pdf',plot = mono_cluster,width = 15,height = 15)

#BEAM统计分析
BEAM_res1 <- BEAM(scRNA.cds[ordering_genes,], branch_point = 1, cores = 1)
BEAM_res1 <- BEAM_res1[order(BEAM_res1$qval),]
BEAM_res1 <- BEAM_res1[,c("gene_short_name", "pval", "qval")]
#选100前个基因可视化
BEAM_genes1<-BEAM_res1[order(-BEAM_res1$qval),][c(1:100),'gene_short_name']

length(BEAM_genes1)
BEAM_p1 <- plot_genes_branched_heatmap(scRNA.cds[BEAM_genes1,], 
                                       branch_point = 1, 
                                       num_clusters = 3, show_rownames = T, 
                                       return_heatmap = T)

BEAM_p1$ph_res

anno_row<-BEAM_p1$annotation_row
BEAM_genes_dat<-t(as.matrix(scRNA.cds@assayData$exprs[BEAM_genes1,]))
meta_clust<-data.frame(cell=colnames(scRNA.cds),
                       cell_type=scRNA.cds$cell_type,
                       seraut_cluster=scRNA.cds$seurat_clusters,
                       Pseudotime=scRNA.cds$Pseudotime)
colnames(BEAM_genes_dat)=gsub('-','__',colnames(BEAM_genes_dat))
BEAM_genes_dat<-merge(meta_clust,
                      data.frame(cell=rownames(BEAM_genes_dat),BEAM_genes_dat),
                      by='cell')
BEAM_genes_dat[1:4,1:4]
rownames(BEAM_genes_dat)=BEAM_genes_dat$cell
BEAM_genes_dat=BEAM_genes_dat[,-1]
BEAM_genes_dat[1:4,1:4]
#排序
BEAM_genes_dat=BEAM_genes_dat[order(BEAM_genes_dat$cell_type,
                                    BEAM_genes_dat$seraut_cluster,
                                    BEAM_genes_dat$Pseudotime),]
anno_row=data.frame(gene=rownames(anno_row),Cluster=anno_row$Cluster)
anno_row=anno_row[order(anno_row$Cluster),]
dim(anno_row)
write.table(anno_row,'state_100_gene.txt',quote = F,sep='\t',row.names = F)
colnames(BEAM_genes_dat)=gsub('__','-',colnames(BEAM_genes_dat))

anno_col=BEAM_genes_dat[,c(1,2,3)]
gene_clust=t(BEAM_genes_dat[,-c(1,2,3)])
bk<-c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
anno_row1=data.frame(Cluster=anno_row$Cluster)
rownames(anno_row1)=anno_row$gene

pheatmap::pheatmap(gene_clust[rownames(anno_row1),],scale = 'row',
                   show_colnames = F,annotation_row = anno_row1,
                   annotation_col = anno_col,
                   show_rownames = F,cluster_rows = F,cluster_cols = F,
                   color=c(colorRampPalette(colors=c("blue","white"))(length(bk)/2),
                           colorRampPalette(colors=c("white","red"))(length(bk)/2)),
                   legend_breaks=seq(-2,2,1),breaks=bk,
                   filename = 'state.gene.pdf')

pheatmap::pheatmap(gene_clust[rownames(anno_row1),],scale = 'row',
                   show_colnames = F,annotation_row = anno_row1,
                   annotation_col = anno_col,
                   show_rownames = F,cluster_rows = F,cluster_cols = F,
                   color=c(colorRampPalette(colors=c("blue","white"))(length(bk)/2),
                           colorRampPalette(colors=c("white","red"))(length(bk)/2)),
                   legend_breaks=seq(-2,2,1),breaks=bk,
                   filename = 'state.gene.png')


head(anno_row)

write.table(anno_row,'state_100_gene.txt',quote = F,row.names = F,sep='\t')

#9、细胞通讯-------------------
load('sce2.RData')
library(Seurat)
library(CellChat)
library(patchwork)

data.input=normalizeData(as.matrix(GetAssayData(sce,slot = "counts")),
                         scale.factor = 10000, do.log = TRUE)
cell_type=sce@meta.data
#读取亚群注释文件，并添加至细胞注释文件中
cell_anno<-read.delim('celltype.txt',sep='\t',header = T,)
colnames(cell_anno)=c('seurat_clusters','cell_type')
cell_type=merge(data.frame(cell=rownames(cell_type),cell_type),
                cell_anno,by='seurat_clusters')
head(cell_type)
rownames(cell_type)=cell_type$cell

#创建cellchat对象
cellchat <- createCellChat(object = data.input, meta = cell_type, group.by = "cell_type")
#添加meta，细胞的注释信息
cellchat <- addMeta(cellchat, meta = cell_type)
#设置细胞的默认的分组
cellchat <- setIdent(cellchat, ident.use = "cell_type") 
levels(cellchat@idents) 
#计算每一个分组的细胞个数
groupSize <- as.numeric(table(cellchat@idents)) 
#导入人的配受体数据库
CellChatDB=CellChatDB.human
##小鼠的
#CellChatDB <- CellChatDB.mouse


#使用配受体数据库中全部的信息进行后续分析
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
#在矩阵的所有的基因中提取signaling gene ,结果保存在data.signaling
cellchat <- subsetData(cellchat) 
future::plan("multiprocess", workers = 4)
#寻找每一个亚群中高表达的配受体
cellchat <- identifyOverExpressedGenes(cellchat)
#将上一步运行的结果储存在cellchat@LR$LRsig
cellchat <- identifyOverExpressedInteractions(cellchat)
#projectData将配受体对的表达值投射到PPI上,来对cellchat@data.signaling 中的表达值进行校正并保存在cellchat@data.project
cellchat <- projectData(cellchat, PPI.human)




#根据表达值推测细胞互做的概率
cellchat <- computeCommunProb(cellchat)
#过滤，某些细胞群中只有少数细胞，则过滤掉细胞间的通信
cellchat <- filterCommunication(cellchat,min.cells=10)
df.net <- subsetCommunication(cellchat)
#推断信号通路水平的细胞通讯网络
cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat,slot.name='netP')


#统计细胞和细胞之间通信的数量（有多少个配体-受体对）和强度（概率）
cellchat <- aggregateNet(cellchat)
#计算每种细胞各有多少个
groupSize <- as.numeric(table(cellchat@idents))
pdf('net_number_strength_new_22.pdf',he=9,wi=15,onefile = F)
par(mfrow = c(1,2)) 
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions") 
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Interaction weights/strength")

dev.off()

# count也可以改权重
mat <- cellchat@net$count
#mat <- cellchat@net$weight
cell_num=length(as.character(unique(cell_type$cell_type)))
cell_col=ceiling(sqrt(cell_num))
cell_row=ceiling(cell_num/cell_col)
pdf('net_weight_individual_new_weight_22.pdf',he=15,wi=15)
par(mfrow = c(cell_row,cell_col), xpd=TRUE) 
for (i in 1:nrow(mat)) { 
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))   
  mat2[i, ] <- mat[i, ]   
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, arrow.width = 0.2,                     
                   arrow.size = 0.1, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}
dev.off()
write.table(df.net,file = 'cellchat_result_new.txt',quote = F,row.names = F,sep='\t')

#单个信号通路或配体-受体介导的细胞互作可视化###

#1、层次图
#查看都有哪些信号通路
cellchat@netP$pathways  
# 选择其中一个信号通路，比如说 MHC-II
pathways.show <- c("VISFATIN")
#有哪些细胞亚群
levels(cellchat@idents)
#选择哪些细胞亚群进行绘图
vertex.receiver = c(3,8) 
p1=netVisual_aggregate(cellchat, signaling = pathways.show, 
                       vertex.receiver = vertex.receiver,
                       layout = 'hierarchy')
p1
pdf('MHC-I_hierarchy.pdf',he=7,wi=15)
print(p1)
dev.off()

#2、网络图
p2=netVisual_aggregate(cellchat, signaling = pathways.show,
                       vertex.receiver=vertex.receiver,
                       layout = "circle")

p2
pdf('MHC-I_circle.pdf',he=7,wi=12)
print(p2)
dev.off()

#3、和弦图（Chord diagram）
pdf('MHC-I_chord.pdf',he=15,wi=15)
netVisual_aggregate(cellchat, signaling = pathways.show,
                    layout = "chord")
dev.off()

#4、热图
pdf('MHC-I_heatmap.pdf',he=7,wi=7)
par(mfrow=c(1,1)) 
netVisual_heatmap(cellchat,
                  signaling = pathways.show, 
                  color.heatmap = "Reds")
dev.off()

#计算配体受体对选定信号通路的贡献值（在这里就是查看哪条配体-受体对MHC-II贡献最大）
p3<-netAnalysis_contribution(cellchat,
                             signaling = pathways.show) 
p3
pairLR.path <- extractEnrichedLR(cellchat,
                                 signaling = pathways.show
                                 , geneLR.return = FALSE) 
pairLR.path
pdf('VISFATIN_ontribution.pdf')
print(p3)
dev.off()
#
#提取对这个通路贡献最大的配体受体对来展示（也可以选择其他的配体受体对）
LR.show <- pairLR.path[1,]
pdf('VISFATIN_hierarchy.pdf',he=7,wi=7)
netVisual_individual(cellchat, 
                     signaling = pathways.show,
                     pairLR.use = LR.show, layout = "circle")
dev.off()
#多个配受体介导的细胞互作关系可视化#####




#有哪些细胞亚群
levels(cellchat@idents)
#配受体介导的亚群-亚群作用的气泡图
netVisual_bubble(cellchat, 
                 sources.use = c(1,7,8),  
                 targets.use = c(3,4,5), 
                 remove.isolate = FALSE)
#查看都有哪些信号通路
cellchat@netP$pathways

#####指定信号通路的
netVisual_bubble(cellchat, 
                 sources.use = c(3,5,7,8,9), 
                 targets.use = c(1,2,4,6),                  
                 signaling = c("SEMA4","GAS"), 
                 remove.isolate = FALSE)
