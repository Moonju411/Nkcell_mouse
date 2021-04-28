library(Seurat)
library(SeuratWrappers)
library(cowplot)
library(dplyr)
library(ggplot2)

sample2_H_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/2-H-Epi/outs/filtered_feature_bc_matrix")
sample2_H_Epi <- CreateSeuratObject(counts = sample2_H_Epi.data, project = "2-H_Epi")
sample2_H_Epi

sample3_H_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/3-H-Epi/outs/filtered_feature_bc_matrix")
sample3_H_Epi <- CreateSeuratObject(counts = sample3_H_Epi.data, project = "3-H_Epi")
sample3_H_Epi

sample4_H_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/4-H-Epi/outs/filtered_feature_bc_matrix")
sample4_H_Epi <- CreateSeuratObject(counts = sample4_H_Epi.data, project = "4-H_Epi")
sample4_H_Epi
####################################################################
sample2_L_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/2-L-Epi/outs/filtered_feature_bc_matrix")
sample2_L_Epi <- CreateSeuratObject(counts = sample2_L_Epi.data, project = "2-L_Epi")
sample2_L_Epi

sample3_L_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/3-L-Epi/outs/filtered_feature_bc_matrix")
sample3_L_Epi <- CreateSeuratObject(counts = sample3_L_Epi.data, project = "3-L_Epi")
sample3_L_Epi

sample4_L_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/4-L-Epi/outs/filtered_feature_bc_matrix")
sample4_L_Epi <- CreateSeuratObject(counts = sample4_L_Epi.data, project = "4-L_Epi")
sample4_L_Epi

Epi.list <- list(sample2_H_Epi, sample3_H_Epi, sample4_H_Epi, sample2_L_Epi, sample3_L_Epi, sample4_L_Epi)

#do Per sample sctransform
for (i in 1:length(Epi.list)) {
  Epi.list[[i]] <- SCTransform(Epi.list[[i]], verbose = FALSE)
}

Epi.features <- SelectIntegrationFeatures(object.list = Epi.list, nfeatures = 3000)

Epi.list <- PrepSCTIntegration(object.list = Epi.list, anchor.features = Epi.features, verbose = TRUE)

Epi.anchors <- FindIntegrationAnchors(object.list = Epi.list, normalization.method = "SCT", anchor.features = Epi.features, verbose = TRUE)
Epi.integrated <- IntegrateData(anchorset = Epi.anchors, normalization.method = "SCT", verbose = TRUE)

Epi.integrated <- RunPCA(Epi.integrated, verbose = FALSE)
Epi.integrated <- RunUMAP(Epi.integrated, dims = 1:30)
Epi.integrated <- FindNeighbors(Epi.integrated, dims = 1:30, verbose = FALSE)
Epi.integrated <- FindClusters(Epi.integrated, verbose = FALSE)

Epi.integrated[['group_id']] <- lapply(list(Epi.integrated@meta.data[,1]), function(x) substr(x,3,8))
FeaturePlot(Epi.integrated, "Lpl", split.by = "group_id")
DimPlot(Epi.integrated, label = TRUE)

DefaultAssay(Epi.integrated) <- "RNA"
Epi.integrated <- NormalizeData(Epi.integrated, verbose = FALSE)
