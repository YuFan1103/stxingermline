library(dplyr)
library(ComplexHeatmap)
library(data.table)
library(SCopeLoomR)
library(SCENIC)
library(AUCell)
library(scales)

scenicLoomPath <- "/lustre/user/taowlab/bijl/fly/pyscenic/tt/wt_SCENIC.loom"

loom <- open_loom(scenicLoomPath)
    # Read information from loom file:
    regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
        regulons <- regulonsToGeneLists(regulons_incidMat)
    regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
    regulonAucThresholds <- get_regulon_thresholds(loom, only.selected = TRUE)

nn <- as.numeric(names(regulonAucThresholds))
names(nn) <- regulonAucThresholds
regulonAucThresholds <- nn



# get regulons
rr <- unlist(regulons) %>% as.data.frame()
colnames(rr) <- "target"
rr$TF <- rownames(rr) %>% gsub("\\(\\+\\).*", "", .)
rownames(rr) <- NULL

co <- fread("/lustre/user/taowlab/bijl/fly/pyscenic/tt/WT.tsv")
co <- as.data.frame(co)
# get intersection of regulons and results of grn
ss <- merge(rr, co, by = c("TF", "target"), all.x = TRUE)


tfs <- ss %>% filter(target == "stx") %>% with(TF)

ff <- ss %>% filter(target != "stx")
indexx <- sample(1:nrow(ff), 1000)
final <- ff[indexx,] 
final <- rbind(final, filter(ss, target == "stx"))


write.table(final, "tt/stx_ss.txt", quote = FALSE, sep = "\t")


# heatmap

# This function will be included in the next version of AUCell
binarizeAUC <- function(auc, thresholds)
{
  thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
  regulonsCells <- setNames(lapply(names(thresholds), 
                                   function(x) {
                                     trh <- thresholds[x]
                                     names(which(getAUC(auc)[x,]>trh))
                                   }),names(thresholds))
  
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"  
  
  return(binaryRegulonActivity)
}


binaryRegulonActivity <- binarizeAUC(regulonAUC, regulonAucThresholds)
dim(binaryRegulonActivity)



nCells <- 1000
set.seed(131)
cellsSelected <- sample(colnames(binaryRegulonActivity), nCells) 
binAct_subset <- binaryRegulonActivity[, which(colnames(binaryRegulonActivity) %in% cellsSelected)]
dim(binAct_subset)


wt <- readRDS("/lustre/user/taowlab/bijl/fly/pyscenic/WT.rds")
gg <- wt@meta.data[colnames(binAct_subset),"seurat_clusters"]

cc <- hue_pal()(length(unique(gg)))
names(cc) <- levels(gg)

ann <- HeatmapAnnotation(Group = gg,
                  #simple_anno_size = unit(2, 'mm'), 
                  col = list(Group = cc),
                  show_annotation_name = FALSE)

tfs <- c("E(bx)(+)", "Trf2(+)", "Zif(+)", "maf-S(+)", "vis(+)", "Cnot4(+)", "Cnx99A(+)",
        "NFAT(+)", "Mnt(+)", "Abl(+)", "E2f1(+)")
# get the matched index
binAct_subset <- binAct_subset[c(tfs, rownames(binAct_subset)[!rownames(binAct_subset) %in% tfs]),]

lab = rowAnnotation(ano = anno_mark(at = 1:length(tfs),
                    labels = tfs,
                    labels_gp = gpar(fontsize = 10, fontface="bold")))

pdf("../autism/tmp/heatmap1.pdf", width = 14, height = 14)
ComplexHeatmap::Heatmap(binAct_subset, name="Binarized activity", 
                      col = c("white", "black"),
                      cluster_rows = TRUE, cluster_columns=TRUE,
                      show_column_names = FALSE,
                      column_title = NULL,
                      row_names_gp=grid::gpar(fontsize=6), # row font size
                      show_row_names = TRUE,
                      top_annotation = ann,
                      #column_split = length(unique(gg)),
                      #right_annotation = lab,
                      show_row_dend = FALSE,
                      show_column_dend = FALSE
                      )

dev.off()













