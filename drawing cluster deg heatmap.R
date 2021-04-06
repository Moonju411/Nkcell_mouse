#install.packages("heatmaply")
library(heatmaply)
library(Seurat)
library(Matrix)

draw_cluster_deg_heatmap <- function(input_seuratobj, thr_p_val_adj=0.01, thr_avg_logFC=1) {
  # cluter per marker_gene_number heatmap
  intput_cluster_number = length(levels(input_seuratobj@active.ident))    #클러스터 갯수(if 0~8이면 9)
  mtx_all <<- matrix(, nrow = intput_cluster_number, ncol = intput_cluster_number)    #FindMarkers결과값을 저장할 empty matrix
  
  for(i in (intput_cluster_number-1):0){    #중첩된 loop문 FindMarkers에 숫자 넣어줌
    for(j in 0:(i-1)){
      if ((i-1) >= 0){
        bim <- FindMarkers(input_seuratobj, ident.1 = i,  ident.2 =j, only.pos = F, test.use = "bimod")
        bim_thresh <- subset(bim, (bim$p_val_adj < thr_p_val_adj)&(abs(bim$avg_logFC) >= thr_avg_logFC))
        num_deg <- dim(bim_thresh)[1]
        
        mtx_all[i+1,j+1] <<- dim(bim_thresh)[1]    #FindMarkers결과값을 mtx_all에 순서에 맞춰서 넣어줌
      }
    }
  }
  mtx_all[is.na(mtx_all)] <<- 0
  
  mtx_sym <<- forceSymmetric(mtx_all,uplo = 'L')
  df_mtx <<- data.frame(mtx_sym)
  
  rownames(df_mtx) <<- c(0:(intput_cluster_number-1))
  colnames(df_mtx) <<- c(0:(intput_cluster_number-1))
  
  library(heatmaply)
  heatmaply(
    df_mtx, 
    symm=TRUE,
    dendrogram = FALSE,
    title= 'Statistically significant DEGs between clusters',
    xlab= "Cluster_id",
    ylab= "Cluster_id",
    cellnote = mtx_all,
    cellnote_textposition='middle center')
}


#For example,
draw_cluster_deg_heatmap(H_Epi.combined_qc, 0.01, 1)
