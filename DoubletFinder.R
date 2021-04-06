library(Seurat)

sample2_H_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/2-H-Epi/outs/filtered_feature_bc_matrix")
sample2_H_Epi <- CreateSeuratObject(counts = sample2_H_Epi.data, project = "2-H_Epi")

sample3_H_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/3-H-Epi/outs/filtered_feature_bc_matrix")
sample3_H_Epi <- CreateSeuratObject(counts = sample3_H_Epi.data, project = "3-H_Epi")

sample4_H_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/4-H-Epi/outs/filtered_feature_bc_matrix")
sample4_H_Epi <- CreateSeuratObject(counts = sample4_H_Epi.data, project = "4-H_Epi")

H_Epi.combined <- merge(sample2_H_Epi, y = c(sample3_H_Epi, sample4_H_Epi), add.cell.ids = c("sample2_H_Epi", "sample3_H_Epi", "sample4_H_Epi"))
H_Epi.combined

sample2_L_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/2-L-Epi/outs/filtered_feature_bc_matrix")
sample2_L_Epi <- CreateSeuratObject(counts = sample2_L_Epi.data, project = "2-L_Epi")

sample3_L_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/3-L-Epi/outs/filtered_feature_bc_matrix")
sample3_L_Epi <- CreateSeuratObject(counts = sample3_L_Epi.data, project = "3-L_Epi")

sample4_L_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/4-L-Epi/outs/filtered_feature_bc_matrix")
sample4_L_Epi <- CreateSeuratObject(counts = sample4_L_Epi.data, project = "4-L_Epi")

L_Epi.combined <- merge(sample2_L_Epi, y = c(sample3_L_Epi, sample4_L_Epi), add.cell.ids = c("sample2_L_Epi", "sample3_L_Epi", "sample4_L_Epi"))
L_Epi.combined

Epi.combined <- merge(H_Epi.combined, c(L_Epi.combined))
Epi.combined[["percent.mt"]] <- PercentageFeatureSet(Epi.combined, pattern = "^mt-")
Epi.combined_qc <- subset(Epi.combined, subset = nCount_RNA < 30000 & nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 30)

Epi.combined_qc <- NormalizeData(Epi.combined_qc, normalization.method = "LogNormalize", scale.factor = 10000)
Epi.combined_qc <- FindVariableFeatures(Epi.combined_qc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Epi.combined_qc)
Epi.combined_qc <- ScaleData(Epi.combined_qc, features = all.genes)

Epi.combined_qc <- RunPCA(Epi.combined_qc, features = VariableFeatures(object = Epi.combined_qc))
DimPlot(Epi.combined_qc, reduction = "pca")

ElbowPlot(H_Epi.combined_qc)

Epi.combined_qc <- RunUMAP(Epi.combined_qc, dims = 1:15)
Epi.combined_qc <- FindNeighbors(Epi.combined_qc, dims = 1:15)
Epi.combined_qc <- FindClusters(Epi.combined_qc, resolution = 0.5)

#Doublet finder
suppressMessages(require(DoubletFinder))

# Can run parameter optimization with paramSweep, but skip for now.
# sweep.res <- paramSweep_v3(data.filt) sweep.stats <- summarizeSweep(sweep.res,
# GT = FALSE) bcmvn <- find.pK(sweep.stats) barplot(bcmvn$BCmetric, names.arg =
# bcmvn$pK, las=2)

# define the expected number of doublet cellscells.
nExp <- round(ncol(Epi.combined_qc) * 0.07)  # expect 7% doublets following 10x vigentes
data.filt <- doubletFinder_v3(Epi.combined_qc, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:15)

# name of the DF prediction can change, so extract the correct column name.
DF.name = colnames(data.filt@meta.data)[grepl("DF.classification", colnames(data.filt@meta.data))]
cowplot::plot_grid(ncol = 2, DimPlot(data.filt, group.by = "orig.ident") + NoAxes(), 
                   DimPlot(data.filt, group.by = DF.name) + NoAxes())

