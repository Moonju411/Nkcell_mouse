library(Seurat)

################################################
sample2_H_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/2-H-Epi/outs/filtered_feature_bc_matrix")
sample2_H_Epi <- CreateSeuratObject(counts = sample2_H_Epi.data, project = "2-H_Epi")

sample3_H_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/3-H-Epi/outs/filtered_feature_bc_matrix")
sample3_H_Epi <- CreateSeuratObject(counts = sample3_H_Epi.data, project = "3-H_Epi")

sample4_H_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/4-H-Epi/outs/filtered_feature_bc_matrix")
sample4_H_Epi <- CreateSeuratObject(counts = sample4_H_Epi.data, project = "4-H_Epi")

H_Epi.combined <- merge(sample2_H_Epi, y = c(sample3_H_Epi, sample4_H_Epi), add.cell.ids = c("sample2_H_Epi", "sample3_H_Epi", "sample4_H_Epi"))
H_Epi.combined

################################################
sample2_L_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/2-L-Epi/outs/filtered_feature_bc_matrix")
sample2_L_Epi <- CreateSeuratObject(counts = sample2_L_Epi.data, project = "2-L_Epi")

sample3_L_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/3-L-Epi/outs/filtered_feature_bc_matrix")
sample3_L_Epi <- CreateSeuratObject(counts = sample3_L_Epi.data, project = "3-L_Epi")

sample4_L_Epi.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/4-L-Epi/outs/filtered_feature_bc_matrix")
sample4_L_Epi <- CreateSeuratObject(counts = sample4_L_Epi.data, project = "4-L_Epi")

L_Epi.combined <- merge(sample2_L_Epi, y = c(sample3_L_Epi, sample4_L_Epi), add.cell.ids = c("sample2_L_Epi", "sample3_L_Epi", "sample4_L_Epi"))
L_Epi.combined

################################################
sample2_H_SQ.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/2-H-SQ/outs/filtered_feature_bc_matrix")
sample2_H_SQ <- CreateSeuratObject(counts = sample2_H_SQ.data, project = "2-H_SQ")

sample3_H_SQ.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/3-H-SQ/outs/filtered_feature_bc_matrix")
sample3_H_SQ <- CreateSeuratObject(counts = sample3_H_SQ.data, project = "3-H_SQ")

sample4_H_SQ.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/4-H-SQ/outs/filtered_feature_bc_matrix")
sample4_H_SQ <- CreateSeuratObject(counts = sample4_H_SQ.data, project = "4-H_SQ")

H_SQ.combined <- merge(sample2_H_SQ, y = c(sample3_H_SQ, sample4_H_SQ), add.cell.ids = c("sample2_H_SQ", "sample3_H_SQ", "sample4_H_SQ"))
H_SQ.combined

################################################
sample2_L_SQ.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/2-L-SQ/outs/filtered_feature_bc_matrix")
sample2_L_SQ <- CreateSeuratObject(counts = sample2_L_SQ.data, project = "2-L_SQ")

sample3_L_SQ.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/3-L-SQ/outs/filtered_feature_bc_matrix")
sample3_L_SQ <- CreateSeuratObject(counts = sample3_L_SQ.data, project = "3-L_SQ")

sample4_L_SQ.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/4-L-SQ/outs/filtered_feature_bc_matrix")
sample4_L_SQ <- CreateSeuratObject(counts = sample4_L_SQ.data, project = "4-L_SQ")

L_SQ.combined <- merge(sample2_L_SQ, y = c(sample3_L_SQ, sample4_L_SQ), add.cell.ids = c("sample2_L_SQ", "sample3_L_SQ", "sample4_L_SQ"))
L_SQ.combined

################################################
sample2_H_SP.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/2-H-SP/outs/filtered_feature_bc_matrix")
sample2_H_SP <- CreateSeuratObject(counts = sample2_H_SP.data, project = "2-H_SP")

sample3_H_SP.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/3-H-SP/outs/filtered_feature_bc_matrix")
sample3_H_SP <- CreateSeuratObject(counts = sample3_H_SP.data, project = "3-H_SP")

sample4_H_SP.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/4-H-SP/outs/filtered_feature_bc_matrix")
sample4_H_SP <- CreateSeuratObject(counts = sample4_H_SP.data, project = "4-H_SP")

H_SP.combined <- merge(sample2_H_SP, y = c(sample3_H_SP, sample4_H_SP), add.cell.ids = c("sample2_H_SP", "sample3_H_SP", "sample4_H_SP"))
H_SP.combined

################################################
sample2_L_SP.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/2-L-SP/outs/filtered_feature_bc_matrix")
sample2_L_SP <- CreateSeuratObject(counts = sample2_L_SP.data, project = "2-L_SP")

sample3_L_SP.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/3-L-SP/outs/filtered_feature_bc_matrix")
sample3_L_SP <- CreateSeuratObject(counts = sample3_L_SP.data, project = "3-L_SP")

sample4_L_SP.data <- Read10X(data.dir = "/node210data/project/mouse_NK_cell/cellranger_out/4-L-SP/outs/filtered_feature_bc_matrix")
sample4_L_SP <- CreateSeuratObject(counts = sample4_L_SP.data, project = "4-L_SP")

L_SP.combined <- merge(sample2_L_SP, y = c(sample3_L_SP, sample4_L_SP), add.cell.ids = c("sample2_L_SP", "sample3_L_SP", "sample4_L_SP"))
L_SP.combined

Epi.combined <- merge(H_Epi.combined, c(L_Epi.combined))
SP.combined <- merge(H_SP.combined, c(L_SP.combined))
SQ.combined <- merge(H_SQ.combined, c(L_SQ.combined))

#Standard workflow QC
Epi.combined[["percent.mt"]] <- PercentageFeatureSet(Epi.combined, pattern = "^mt-")
VlnPlot(Epi.combined, features = c("nFeature_RNA"))
VlnPlot(Epi.combined, features = c("nCount_RNA"))
VlnPlot(Epi.combined, features = c("percent.mt"))

Epi.combined[["percent.lpl"]] <- PercentageFeatureSet(Epi.combined, pattern = "^Lpl")
VlnPlot(Epi.combined, features = c("percent.lpl"))

Epi.combined_qc <- subset(Epi.combined, subset = nCount_RNA < 30000 & nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 30)
VlnPlot(Epi.combined_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Epi.combined_qc, features = c("percent.lpl"))

#Save seurat obj to rds
saveRDS(Epi.combined_qc, file = "/node210data/project/mouse_NK_cell/cellranger_out/ALL_combined_qc.rds")
