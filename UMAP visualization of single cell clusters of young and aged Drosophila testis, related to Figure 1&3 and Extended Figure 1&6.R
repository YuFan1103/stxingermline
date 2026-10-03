library(Seurat)
library(dplyr)
library(cowplot)
library(Matrix)
library(ggplot2)
library(patchwork)
library(future)
library(harmony)
library(SCpubr)

getwd()

plan("multisession", workers = 12)
count <- Read10X(data.dir = "/Users/admin/Documents/T1_T60_filtered_feature_bc_matrix")
cells <- colnames(count)

#
T1_data1 <- count[,grepl('1$', cells)]  
T60_data1 <- count[,grepl('2$', cells)]  

T1_1 <- CreateSeuratObject(counts = T1_data1, project = "T1",min.cells = 2, min.features = 200)
T1_1@meta.data$group_ID <- "T1_1"

T60_1 <- CreateSeuratObject(counts = T60_data1, project = "T60",min.cells = 2, min.features = 200)
T60_1@meta.data$group_ID <- "T60_1"

T60_1 <- subset(T60_1,cells = sample(colnames(T60_1),8000))

#T1.combined.raw <- merge(T1_1, y = c(T60_1, T120_1),  project = "T1")
#Quality control
T1_1[["percent.mt"]] <- PercentageFeatureSet(T1_1,pattern = "^mt")
T1_1[["percent.ribo"]] <- PercentageFeatureSet(T1_1,pattern = "^Rp[SL]")
T60_1[["percent.mt"]] <- PercentageFeatureSet(T60_1,pattern = "^mt")
T60_1[["percent.ribo"]] <- PercentageFeatureSet(T60_1,pattern = "^Rp[SL]")

VlnPlot(T1_1,features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 5)
VlnPlot(T60_1,features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 5)

T1_1 <- subset(T1_1, nFeature_RNA>200 & percent.mt <15 ) 
T60_1 <- subset(T60_1, nFeature_RNA>800 & percent.mt <5 ) 
T1.combined.filted <- merge(T1_1, y = T60_1,  project = "T1_T60")

VlnPlot(T1.combined.raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 5)
T1.combined.filted <- subset(T1.combined.raw, nFeature_RNA>500 & percent.mt <5 ) 
VlnPlot(T1.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 5)

options(future.globals.maxSize = 1000 * 1024^2*20)

T1.combined <- NormalizeData(T1.combined.filted, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(verbose=TRUE)

#Run Harmony to correct batch effects
T1.combined <- RunHarmony (T1.combined,"orig.ident")

#elbowplot decide dims

T1.combined <- FindNeighbors(T1.combined,  dims = 1:20) %>% FindClusters(resolution = 0.2)

T1.combined <- RunUMAP(T1.combined,  dims = 1:20, reduction ="harmony")

DimPlot(T1.combined, reduction = "umap",group.by = "orig.ident",cols = c("#F6A69F","#196EC9"), pt.size = 0.1)
do_DimPlot(T1.combined,border.color = "white",border.size = 0.1,pt.size = 0.7)
do_FeaturePlot(T1.combined,features = "stx")

#Subclusters of cyst cells and cyst stem cells
sce=T1.combined
Idents(sce)
levels(sce)
head(sce@meta.data)
genes_to_check = c('stx')
DotPlot(sce, group.by = 'seurat_clusters',
        features = unique(genes_to_check)) + RotatedAxis()
p1=DimPlot(sce, reduction = 'tsne', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()
p2=DotPlot(sce, group.by = 'seurat_clusters',
           features = unique(genes_to_check)) + RotatedAxis()

p1+p2
do_DimPlot(T1.combined,reduction = "umap")

stx_sce1 = sce[,sce@meta.data$seurat_clusters %in% c(8)]
stx_sce2 = sce[, Idents(sce) %in% c( "8")]

sce=stx_sce1
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1e4) 
sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
sce <- FindNeighbors(sce, dims = 1:20)
sce <- FindClusters(sce, resolution = 0.3) #resolution of 0.3 can distinguish CySCs #0 is CySCs

head(Idents(sce), 5)

sce <- RunUMAP(sce, dims = 1:20)
DimPlot(sce, reduction = "umap")
CySC_aging_markers <- FindMarkers(sce, ident.1 = "T60",ident.2 = "T1", group.by = "orig.ident", subset.ident = 0)
de_genes <- CySC_aging_markers
do_VolcanoPlot(sce, de_genes = de_genes, n_genes = 5, FC_cutoff = 0.2, pt.size = 0.8, colors.use = "darkgreen")
write.table(CySC_aging_markers, 'CySC_aging_markers.xls', row.names=T, col.names=NA, sep="\t")
do_ViolinPlot(sce,features = "stx", group.by = "orig.ident")

FeaturePlot(sce, features = c("tj","zfh1","eya","stx"))

#plots and tables （not all plots were used in the paper）
do_BoxPlot(sce, feature = "stx", group.by = "orig.ident")
FeaturePlot(sce, features = c("zfh1","eya"))
do_ViolinPlot(sce, feature = "stx", group.by = "orig.ident")
VlnPlot(sce,feature = "stx", group.by = "orig.ident",cols =c("#9579DB","#4C9A7D"),pt.size = 1)
VlnPlot(sce,feature = c("tj","zfh1","aub","zpg"), group.by = "seurat_clusters", cols =c("#9579DB","#4C9A7D"),pt.size = 1)
do_ViolinPlot(sce,feature = c("tj","zfh1","aub","zpg"), group.by = "seurat_clusters", cols =c("#9579DB","#4C9A7D"),pt.size = 1)
do_NebulosaPlot(T1.combined,features = c("stg","zpg","aub","His2Av","bam","kmg"),viridis_color_map = "E")
do_NebulosaPlot(T1.combined,features = c("stg","zpg","aub","His2Av","bam","kmg","Rbp4","aly","can","sa","CycB","fzo","twe","Dic61B"))
Harmony_markers <- FindAllMarkers(T1.combined,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10markers<-Harmony_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
do_ExpressionHeatmap(T1.combined, feature=top10markers$gene, viridis_color_map = "A")
DoHeatmap(T1.combined, feature=top10markers$gene)
write.table(Harmony_markers, 'Harmony_markers.xls', row.names=T, col.names=NA, sep="\t")
DimPlot(T1.combined, group.by = "orig.ident",cols =c("#68A1E5","#F49F97"),pt.size = 0.1)
do_NebulosaPlot(sce,features = c("zfh1","tj","Stat92E","piwi","eya"),pt.size = 0.2,border.color = "grey")
do_ExpressionHeatmap(sce, features = c("stx","zfh1","Stat92E","tj"))


