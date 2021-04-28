library(SingleCellExperiment)
library(Seurat)
library(slingshot)
library(SeuratWrappers)
library(viridis)
library(RColorBrewer)

Epi.integrated <- readRDS('/node210data/project/mouse_NK_cell/processed_data/Epi.sct_integrated.rds')
Epi.integrated

Epi.integrated_sce <- as.SingleCellExperiment(Epi.integrated)
DimPlot(Epi.integrated, pt.size = 0.5, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

Epi.integrated.sling <- slingshot(Epi.integrated_sce, reducedDim= 'UMAP', clusterLabels = Epi.integrated$seurat_clusters, 
                                  approx_points=100)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
colors2 <- brewer.pal(9,"Set1")[Epi.integrated$seurat_clusters]
plotcol <- colors[cut(Epi.integrated.sling$slingPseudotime_1, breaks=100)]

plot(reducedDims(Epi.integrated.sling)$UMAP, col = colors2, pch=16, asp = 0.5)
lines(SlingshotDataSet(Epi.integrated.sling), lwd = 1, type = 'lineage', col = 'black')
