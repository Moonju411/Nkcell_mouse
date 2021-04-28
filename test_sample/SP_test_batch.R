library(Seurat)
library(SeuratWrappers)
library(cowplot)
library(dplyr)
library(ggplot2)

sample2_H_SP.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/2-H-SP/outs/filtered_feature_bc_matrix")
sample2_H_SP <- CreateSeuratObject(counts = sample2_H_SP.data, project = "2-H_SP")

sample3_H_SP.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/3-H-SP/outs/filtered_feature_bc_matrix")
sample3_H_SP <- CreateSeuratObject(counts = sample3_H_SP.data, project = "3-H_SP")

sample4_H_SP.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/4-H-SP/outs/filtered_feature_bc_matrix")
sample4_H_SP <- CreateSeuratObject(counts = sample4_H_SP.data, project = "4-H_SP")

H_SP.combined <- merge(sample2_H_SP, y = c(sample3_H_SP, sample4_H_SP), add.cell.ids = c("sample2_H_SP", "sample3_H_SP", "sample4_H_SP"))
H_SP.combined

H_SP.combined[["percent.mt"]] <- PercentageFeatureSet(H_SP.combined, pattern = "^mt-")
H_SP.combined <- subset(H_SP.combined, subset = nCount_RNA < 30000 & nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 30)
H_SP.combined <- NormalizeData(H_SP.combined, verbose = FALSE)
H_SP.combined <- FindVariableFeatures(H_SP.combined)
H_SP.combined <- RunFastMNN(object.list = SplitObject(H_SP.combined, split.by = "orig.ident"))
H_SP.combined <- RunUMAP(H_SP.combined, reduction = "mnn", dims = 1:15)
H_SP.combined <- FindNeighbors(H_SP.combined, reduction = "mnn", dims = 1:15)
H_SP.combined <- FindClusters(H_SP.combined, resolution = 0.5)

DimPlot(H_SP.combined, reduction = "umap",split.by = "orig.ident")
H_SP.combined[['seurat_clusters_with']] <- lapply(list(H_SP.combined@meta.data[,6]), function(x) paste(x,"_H"))
################################################
sample2_L_SP.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/2-L-SP/outs/filtered_feature_bc_matrix")
sample2_L_SP <- CreateSeuratObject(counts = sample2_L_SP.data, project = "2-L_SP")

sample3_L_SP.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/3-L-SP/outs/filtered_feature_bc_matrix")
sample3_L_SP <- CreateSeuratObject(counts = sample3_L_SP.data, project = "3-L_SP")

sample4_L_SP.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/4-L-SP/outs/filtered_feature_bc_matrix")
sample4_L_SP <- CreateSeuratObject(counts = sample4_L_SP.data, project = "4-L_SP")

L_SP.combined <- merge(sample2_L_SP, y = c(sample3_L_SP, sample4_L_SP), add.cell.ids = c("sample2_L_SP", "sample3_L_SP", "sample4_L_SP"))
L_SP.combined

L_SP.combined[["percent.mt"]] <- PercentageFeatureSet(L_SP.combined, pattern = "^mt-")
L_SP.combined <- subset(L_SP.combined, subset = nCount_RNA < 30000 & nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 30)
L_SP.combined <- NormalizeData(L_SP.combined, verbose = FALSE)
L_SP.combined <- FindVariableFeatures(L_SP.combined)
L_SP.combined <- RunFastMNN(object.list = SplitObject(L_SP.combined, split.by = "orig.ident"))
L_SP.combined <- RunUMAP(L_SP.combined, reduction = "mnn", dims = 1:15)
L_SP.combined <- FindNeighbors(L_SP.combined, reduction = "mnn", dims = 1:15)
L_SP.combined <- FindClusters(L_SP.combined, resolution = 0.5)

DimPlot(L_SP.combined, reduction = "umap",split.by = "orig.ident")
L_SP.combined[['seurat_clusters_with']] <- lapply(list(L_SP.combined@meta.data[,6]), function(x) paste(x,"_L"))

SP_batchremoved_merge <- merge(H_SP.combined, L_SP.combined)
SP_batchremoved_merge <- NormalizeData(SP_batchremoved_merge, verbose = FALSE)
SP_batchremoved_merge <- FindVariableFeatures(SP_batchremoved_merge)

all.genes <- rownames(SP_batchremoved_merge)
SP_batchremoved_merge <- ScaleData(SP_batchremoved_merge, features = all.genes)
SP_batchremoved_merge <- RunPCA(SP_batchremoved_merge, features = VariableFeatures(object = SP_batchremoved_merge))
SP_batchremoved_merge <- RunUMAP(SP_batchremoved_merge, reduction = "pca", dims = 1:15)

SP_batchremoved_merge[['group_id']] <- lapply(list(SP_batchremoved_merge@meta.data[,1]), function(x) substr(x,3,8))

p1 <- DimPlot(SP_batchremoved_merge, group.by = "seurat_clusters_with", split.by = "group_id")

ggplot2::ggsave(filename = "SP_batchremoved_merge.png", plot = p1, dpi = 600)

DimPlot(SP_batchremoved_merge, )
