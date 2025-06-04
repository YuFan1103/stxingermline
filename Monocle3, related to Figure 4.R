library(monocle3)
library(Seurat)

data <- GetAssayData(rGSC, assay = "RNA", slot = "counts")
cell_metadata <- rGSC@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))

rownames(gene_annotation) <- rownames(data)

cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 100)

plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds,preprocess_method = "PCA")

plot_cells(cds, reduction_method = "UMAP") 

cds <- cluster_cells(cds, reduction_method = "UMAP")

cds <- learn_graph(cds, use_partition = TRUE)

FeaturePlot(cds, features = "aub")

cds <- order_cells(cds)

ciliated_genes <- c("S","Myc","Atg8a")

plot_cells(cds, genes=ciliated_genes,    label_cell_groups=FALSE,    show_trajectory_graph=FALSE, cell_size = 1)


-----------------------------------------------------------------------------------
#plots
plot_cells(
    cds = cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE,
    group_cells_by = 'cluster',cell_size = 1,label_leaves=FALSE,) + scale_color_viridis(option = "D")


plot_cells(cds,
           group_label_size = 6,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 3,reduction_method = "UMAP",cell_size = 1)




  
