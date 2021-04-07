library(Seurat)
library(cowplot)

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

############################################
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

# Set up control object
H_Epi.combined[["percent.mt"]] <- PercentageFeatureSet(H_Epi.combined, pattern = "^mt-")
H_Epi_qc <- subset(H_Epi.combined, subset = nCount_RNA < 30000 & nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 30)
H_Epi_qc <- NormalizeData(H_Epi_qc, verbose = FALSE)
H_Epi_qc <- FindVariableFeatures(H_Epi_qc, selection.method = "vst", nfeatures = 2000)

# Set up control object
L_Epi.combined[["percent.mt"]] <- PercentageFeatureSet(L_Epi.combined, pattern = "^mt-")
L_Epi_qc <- subset(L_Epi.combined, subset = nCount_RNA < 30000 & nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 30)
L_Epi_qc <- NormalizeData(L_Epi_qc, verbose = FALSE)
L_Epi_qc <- FindVariableFeatures(L_Epi_qc, selection.method = "vst", nfeatures = 2000)

# Perform integration
Epi.anchors <- FindIntegrationAnchors(object.list = list(H_Epi_qc, L_Epi_qc), dims = 1:20)
Epi.integrated <- IntegrateData(anchorset = Epi.anchors, dims = 1:20)

DefaultAssay(Epi.integrated) <- "integrated"
Epi.integrated[['group_id']] <- lapply(list(Epi.integrated@meta.data[,1]), function(x) substr(x,3,8))

# Run the standard workflow for visualization and clustering
Epi.integrated <- ScaleData(Epi.integrated, verbose = FALSE)
Epi.integrated <- RunPCA(Epi.integrated, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
Epi.integrated <- RunUMAP(Epi.integrated, reduction = "pca", dims = 1:15)
Epi.integrated <- FindNeighbors(Epi.integrated, reduction = "pca", dims = 1:15)
Epi.integrated <- FindClusters(Epi.integrated, resolution = 0.5)

# Visualization
p1 <- DimPlot(Epi.integrated, reduction = "umap", group.by = "group_id")
p2 <- DimPlot(Epi.integrated, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

# Visualization splited by group_id
DimPlot(Epi.integrated, reduction = "umap", split.by = "group_id")
                                       
# For performing differential expression after integration, we switch back to the original data
DefaultAssay(Epi.integrated) <- "RNA"
                                       
# FeaturePlot
FeaturePlot()

# Dotplot
markers.to.plot <- c("Gzmb","Prf1","Gzma","Cpt1b",'Fabp5',"Lpl","Cd36",'Fabp4')
DotPlot(Epi.integrated, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "group_id") + 
  RotatedAxis()

# CorrelationPlot
cluster0_H <- subset(Epi.integrated, (seurat_clusters == '0') & (group_id == 'H_Epi'))
cluster0_L <- subset(Epi.integrated, (seurat_clusters == '0') & (group_id == 'L_Epi'))
Idents(cluster0_H) <- "H_Epi"
Idents(cluster0_L) <- "L_Epi"
df_cluster0_H <- as.data.frame(log1p(AverageExpression(cluster0_H, verbose = FALSE)$RNA))
df_cluster0_L <- as.data.frame(log1p(AverageExpression(cluster0_L, verbose = FALSE)$RNA))
df_cluster0 <- cbind(df_cluster0_H, df_cluster0_L)
df_cluster0$gene <- rownames(df_cluster0)
df_cluster0$average_fc <- df_cluster0$H_Epi/df_cluster0$L_Epi

genes.to.label <- c("Gzmb","Prf1","Gzma","Cpt1b",'Fabp5',"Lpl","Cd36",'Fabp4')
p1 <- ggplot(df_cluster0, aes(H_Epi, L_Epi)) + geom_point()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p1
                                      
                                               
                                       
