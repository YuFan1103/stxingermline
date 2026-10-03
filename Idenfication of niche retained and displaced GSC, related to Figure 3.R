#code for Figure 4 
sce=T1.combined
Idents(sce)
levels(sce)
head(sce@meta.data)
genes_to_check = c('stx')
DotPlot(sce, group.by = 'seurat_clusters',
        features = unique(genes_to_check),cols = c("white","#3300ff")) + RotatedAxis()
p1=DimPlot(sce, reduction = 'umap', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()
p2=DotPlot(sce, group.by = 'seurat_clusters',
           features = unique(genes_to_check)) + RotatedAxis()
do_NebulosaPlot(sce, features = "stx")

p1+p2
do_DimPlot(T1.combined,reduction = "umap")

#Subset of early germ cells of cluster 2
stx_sce1 = sce[,sce@meta.data$seurat_clusters %in% c(2)]
stx_sce2 = sce[, Idents(sce) %in% c( "2")]
sce=stx_sce1
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1e4) 
sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
sce <- FindNeighbors(sce, dims = 1:20)
sce <- FindClusters(sce, resolution = 0.1)

sce <- RunUMAP(sce, dims = 1:20)
DimPlot(sce, reduction = "umap")
FeaturePlot(sce, features = c("stg","zpg","aub","bam","kmg","aly"))   

#cluster1 is GSCs and spermatogonia,stg and zpg is highly expressed 
#subset of cluster1
stx_sce1 = sce[,sce@meta.data$seurat_clusters %in% c(1)]
stx_sce2 = sce[, Idents(sce) %in% c( "1")]
sce=stx_sce1
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1e4) 
sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
sce <- FindNeighbors(sce, dims = 1:30)
sce <- FindClusters(sce, resolution = 0.1)

sce <- RunUMAP(sce, dims = 1:30)
DimPlot(sce, reduction = "umap")
FeaturePlot(sce, features = c("stg","zpg","aub","bam","kmg","aly"))

#subset of cluster0 in cluster1 - identify GSC
stx_sce1 = sce[,sce@meta.data$seurat_clusters %in% c(0)]
stx_sce2 = sce[, Idents(sce) %in% c( "0")]
sce=stx_sce1
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1e4) 
sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
sce <- FindNeighbors(sce, dims = 1:30)
sce <- FindClusters(sce, resolution = 0.2)

#ientify stg and zpg positive GSCs
stx_sce1 = sce[,sce@meta.data$seurat_clusters %in% c(0)]
stx_sce2 = sce[, Idents(sce) %in% c( "0")]

sce=stx_sce1
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1e4) 
sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
sce <- FindNeighbors(sce, dims = 1:30)
sce <- FindClusters(sce, resolution = 0.5)
sce <- RunUMAP(sce, dims = 1:30)

cluster_cells <- Cells(sce)[sce$orig.ident == "T60"]
sampled_cells <- sample(cluster_cells, size = 40,replace = FALSE)
cluster2_cells <- Cells(sce)[sce$orig.ident == "T1"]
combined_cells <- c(sampled_cells, cluster2_cells)
subset_seurat <- subset(sce, cells = combined_cells)
DimPlot(subset_seurat,reduction = "umap")

#plots for figure 4 (not all plots were used in paper)

FeaturePlot(subset_seurat, genes = c("stg","piwi","how","zpg","bam","kmg","aly","Rbp4"))
do_FeaturePlot(subset_seurat, features = c("stg","piwi","how","zpg","bam","kmg","aly","Rbp4"))
do_FeaturePlot(subset_seurat, features = c("stg","piwi","how","zpg","bam","kmg","aly","Rbp4"),use_viridis = TRUE, viridis.palette = "D",border.size = 0)
FeaturePlot(subset_seurat, features = c("stg","piwi","how","zpg","bam","kmg","aly","Rbp4"), cols = c("#450657","#f6e729"))
FeaturePlot_scCustom(subset_seurat,features = c("stg","Drep2","aub","ovo","Pxt","cher","Mov10","Rbp4"),colors_use = viridis_light_high)

# markers for retained and displaced GSCs
Adherent_markers <- FindMarkers(sce, ident.1 = 0,ident.2 = 1, logfc.threshold = 0.25,subset.ident = "T60")
de_genes <- Adherent_markers
do_VolcanoPlot(sce, de_genes = de_genes, n_genes = 2, pval_cutoff = 0.05,FC_cutoff = 0.25, pt.size = 4, colors.use = "darkgreen",border.size = 0.1, order_tags_by = "logfc")
write.table(Adherent_markers, 'adherent_markers-2023.xls', row.names=T, col.names=NA, sep="\t")

