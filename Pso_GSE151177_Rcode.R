#对单细胞数据进行质控
rm(list=ls())
setwd("~/Psoriasis/")
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
dir_name=c('P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11','P12','P13')
datalist=list()
for (i in 1:length(dir_name)){
  dir.10x = paste0("GSE151177/",dir_name[i])
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
              nFeature_RNA < 5000 & 
              quantile(percent.mt, 0.98) > percent.mt & percent.mt < 25 &
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
meta1[grep("P1",meta1$Sample),]$Group1="Pso"
meta1[grep("P2",meta1$Sample),]$Group1="Pso"
meta1[grep("P3",meta1$Sample),]$Group1="Pso"
meta1[grep("P4",meta1$Sample),]$Group1="Pso"
meta1[grep("P5",meta1$Sample),]$Group1="Pso"
meta1[grep("P6",meta1$Sample),]$Group1="Pso"
meta1[grep("P7",meta1$Sample),]$Group1="Pso"
meta1[grep("P8",meta1$Sample),]$Group1="Pso"
meta1[grep("P9",meta1$Sample),]$Group1="Pso"
meta1[grep("P10",meta1$Sample),]$Group1="Pso"
meta1[grep("P11",meta1$Sample),]$Group1="Pso"
meta1[grep("P12",meta1$Sample),]$Group1="Pso"
meta1[grep("P13",meta1$Sample),]$Group1="Pso"

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


#选择锚点和分辨率
Dims <- 40
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
                  #reduction="umap",
                  reduction="tsne",
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
#NK：GNLY  KLRB1
VlnPlot(object=sce, features = c("GNLY","KLRB1"), pt.size = 0, assay = "RNA") + NoLegend()
#CD161_T_cell：CD3D TRAC
VlnPlot(object=sce, features = c("CD3D","TRBC1","TRAC","GZMK","CD8A","GNLY"), pt.size = 0, assay = "RNA") + NoLegend()
#Treg：TRBC1 CD3D TRAC
VlnPlot(object=sce, features = c("TRBC1"), pt.size = 0, assay = "RNA") + NoLegend()
#CD8_T_cell：GZMK CD3D TRAC CD8A CD8B
VlnPlot(object=sce, features = c("CD8B"), pt.size = 0,assay = "RNA") + NoLegend()
#Mature_DC：HLA-DQA1 HLA-DQB1 HLA-DRB1 LAMP3 LY75 CIITA CD40 HLA-DRA
VlnPlot(object=sce, features = c("HLA-DQA1","LAMP3","LYZ","LY75"), pt.size = 0,assay = "RNA") + NoLegend()
#Semimature_DC:LYZ CD14
VlnPlot(object=sce, features = c("CD14"), pt.size = 0,assay = "RNA") + NoLegend()
#Macrophage: SPRR2G LCE3D CD163 CDSN
VlnPlot(object=sce, features = c("CDSN","CD163","LCE3D","SPRR2G","CD14"), pt.size = 0,assay = "RNA") + NoLegend()
#Melanocyte: DCT TYRP1 MLANA
VlnPlot(object=sce, features = c("DCT"), pt.size = 0,assay = "RNA") + NoLegend()
#KC_S.Granulosum:FABP5 KRT10
VlnPlot(object=sce, features = c("KRT10","FABP5"), pt.size = 0,assay = "RNA") + NoLegend()
#KC_S.Spinosum:KRT1
VlnPlot(object=sce, features = c("KRT1"), pt.size = 0,assay = "RNA") + NoLegend()
#KC_S.Basale:KRT14
VlnPlot(object=sce, features = c("KRT14"), pt.size = 0,assay = "RNA") + NoLegend()
#KC_S.Corneum:SPRR2G LCE3D CDSN
VlnPlot(object=sce, features = c("SPRR2G","LCE3D","KRT10","FABP5"), pt.size = 0, assay = "RNA", ncol = 1) + NoLegend()
#左上角
VlnPlot(object=sce, features = c("KRT8","KRT5"), pt.size = 0,assay = "RNA", ncol = 1) + NoLegend()
celltype<-read.delim('Pso_celltype.txt',sep='\t',header = T,check.names = F)
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
                    "11" = celltype[12,2])

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
ggplot2::ggsave('umap_celltype.pdf',plot = umap_celltype,he=7,wi=10)

tsne_celltype<-DimPlot(sce,cols=mycolor,
                       reduction="tsne",
                       label = "T", 
                       pt.size = 0.2,
                       label.size = 5)

tsne_celltype
ggplot2::ggsave('tsne_celltype.pdf',plot = tsne_celltype,he=7,wi=10)

save(sce,file = 'sce3.RData')
# 单细胞中基因的分布---------------
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(RColorBrewer)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(harmony)
library(ggraph)

load("sce3.RData")
#TYMP, PTTG1, IRF7, CD274, LTF, S100A8, S100A9
FeaturePlot(sce, reduction = "tsne", features = c("S100A9"), slot = "scale.data", 
            coord.fixed = T, order = T, cols = c("lightgrey", "brown2"), min.cutoff = "q10", max.cutoff = "q90")



DotPlot(subset(sce), 
        features = c("TYMP","PTTG1","IRF7","CD274","LTF","S100A8","S100A9"), assay='RNA') +  coord_flip() + scale_color_viridis()




