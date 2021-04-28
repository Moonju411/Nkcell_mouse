#install.packages("heatmaply")
library(heatmaply)
library(Seurat)
library(Matrix)

#seurat이 계산한 deg df도 저장하게 코드 수정함.
draw_cluster_deg_heatmap_ver2 <- function(input_seuratobj, thr_p_val_adj=0.01, thr_avg_logFC=1) {
  # cluter per marker_gene_number heatmap
  output_deg_per_clusters <<- list()
  intput_cluster_number = length(levels(input_seuratobj@active.ident))    #클러스터 갯수(if 0~8이면 9)
  mtx_all <<- matrix(, nrow = intput_cluster_number, ncol = intput_cluster_number)    #FindMarkers결과값을 저장할 empty matrix
  
  for(i in (intput_cluster_number-1):0){    #중첩된 loop문 FindMarkers에 숫자 넣어줌
    for(j in 0:(i-1)){
      if ((i-1) >= 0){
        bim <- FindMarkers(input_seuratobj, ident.1 = i,  ident.2 =j, only.pos = F, test.use = "bimod")
        bim_thresh <- subset(bim, (bim$p_val_adj < thr_p_val_adj)&(abs(bim$avg_logFC) >= thr_avg_logFC))
        output_deg_per_clusters[[length(output_deg_per_clusters)+1]] <<- bim_thresh
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

#위에 함수 돌리면 output_deg_per_clusters에 df들이 순서대로 들어가는데, 그거 index를 찾는 함수임
Find_deg_df_key <- function(input_seuratobj, ident1, ident2){
  idx = 1
  intput_cluster_number = length(levels(input_seuratobj@active.ident))    #클러스터 갯수(if 0~8이면 9)
  mtx_key <<- matrix(, nrow = intput_cluster_number, ncol = intput_cluster_number) 
  for(i in (intput_cluster_number-1):0){    
    for(j in 0:(i-1)){
      if ((i-1) >= 0){
        mtx_key[i+1,j+1] <<- idx
        idx <- idx+1
      }
    }
  }
  mtx_key <<- forceSymmetric(mtx_key,uplo = 'L')
  return(mtx_key[ident1+1, ident2+1])
}

#For example,
draw_cluster_deg_heatmap(H_Epi.combined_qc, 0.01, 1)
#만약에 0번과 10번 cluster간의 deg를 보고 싶으면,
output_deg_per_clusters[[Find_deg_df_key(Epi.integrated.F,0,10)]] # 순서 상관없음.
