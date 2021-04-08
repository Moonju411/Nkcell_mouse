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

H_Epi.combined <- merge(sample2_H_Epi, y = c(sample3_H_Epi, sample4_H_Epi), add.cell.ids = c("sample2_H_Epi", "sample3_H_Epi", "sample4_H_Epi"))
H_Epi.combined

H_Epi.combined[["percent.mt"]] <- PercentageFeatureSet(H_Epi.combined, pattern = "^mt-")
H_Epi.combined <- subset(H_Epi.combined, subset = nCount_RNA < 30000 & nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 30)
H_Epi.combined <- NormalizeData(H_Epi.combined, verbose = FALSE)
H_Epi.combined <- FindVariableFeatures(H_Epi.combined)
H_Epi.combined <- RunFastMNN(object.list = SplitObject(H_Epi.combined, split.by = "orig.ident"))
H_Epi.combined <- RunUMAP(H_Epi.combined, reduction = "mnn", dims = 1:15)
H_Epi.combined <- FindNeighbors(H_Epi.combined, reduction = "mnn", dims = 1:15)
H_Epi.combined <- FindClusters(H_Epi.combined, resolution = 0.5)

DimPlot(H_Epi.combined, reduction = "umap",split.by = "orig.ident")
H_Epi.combined[['seurat_clusters_with']] <- lapply(list(H_Epi.combined@meta.data[,6]), function(x) paste(x,"_H"))

################################################
sample2_L_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/2-L-Epi/outs/filtered_feature_bc_matrix")
sample2_L_Epi <- CreateSeuratObject(counts = sample2_L_Epi.data, project = "2-L_Epi")
sample2_L_Epi

sample3_L_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/3-L-Epi/outs/filtered_feature_bc_matrix")
sample3_L_Epi <- CreateSeuratObject(counts = sample3_L_Epi.data, project = "3-L_Epi")
sample3_L_Epi

sample4_L_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/4-L-Epi/outs/filtered_feature_bc_matrix")
sample4_L_Epi <- CreateSeuratObject(counts = sample4_L_Epi.data, project = "4-L_Epi")
sample4_L_Epi

L_Epi.combined <- merge(sample2_L_Epi, y = c(sample3_L_Epi, sample4_L_Epi), add.cell.ids = c("sample2_L_Epi", "sample3_L_Epi", "sample4_L_Epi"))
L_Epi.combined

L_Epi.combined[["percent.mt"]] <- PercentageFeatureSet(L_Epi.combined, pattern = "^mt-")
L_Epi.combined <- subset(L_Epi.combined, subset = nCount_RNA < 30000 & nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 30)
L_Epi.combined <- NormalizeData(L_Epi.combined, verbose = FALSE)
L_Epi.combined <- FindVariableFeatures(L_Epi.combined)
L_Epi.combined <- RunFastMNN(object.list = SplitObject(L_Epi.combined, split.by = "orig.ident"))
L_Epi.combined <- RunUMAP(L_Epi.combined, reduction = "mnn", dims = 1:15)
L_Epi.combined <- FindNeighbors(L_Epi.combined, reduction = "mnn", dims = 1:15)
L_Epi.combined <- FindClusters(L_Epi.combined, resolution = 0.5)

DimPlot(L_Epi.combined, reduction = "umap",split.by = "orig.ident")
L_Epi.combined[['seurat_clusters_with']] <- lapply(list(L_Epi.combined@meta.data[,6]), function(x) paste(x,"_L"))

Epi_batchremoved_merge <- merge(H_Epi.combined, L_Epi.combined)
Epi_batchremoved_merge <- NormalizeData(Epi_batchremoved_merge, verbose = FALSE)
Epi_batchremoved_merge <- FindVariableFeatures(Epi_batchremoved_merge)

all.genes <- rownames(Epi_batchremoved_merge)
Epi_batchremoved_merge <- ScaleData(Epi_batchremoved_merge, features = all.genes)
Epi_batchremoved_merge <- RunPCA(Epi_batchremoved_merge, features = VariableFeatures(object = Epi_batchremoved_merge))
Epi_batchremoved_merge <- RunUMAP(Epi_batchremoved_merge, reduction = "pca", dims = 1:15)

Epi_batchremoved_merge <- FindNeighbors(Epi_batchremoved_merge, reduction = "pca", dims = 1:15)
Epi_batchremoved_merge <- FindClusters(Epi_batchremoved_merge, resolution = 0.5)

Epi_batchremoved_merge[['group_id']] <- lapply(list(Epi_batchremoved_merge@meta.data[,1]), function(x) substr(x,3,8))

p1 <- DimPlot(Epi_batchremoved_merge, group.by = "seurat_clusters", split.by = "group_id")
ggplot2::ggsave(filename = "Epi_batchremoved_merge.png", plot = p1, dpi = 600)

DimPlot(Epi_batchremoved_merge, )
