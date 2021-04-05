# cluter per marker_gene_number heatmap
intput_cluster_number = length(levels(H_Epi.combined_qc@active.ident))    #클러스터 갯수(if 0~8이면 9)
mtx_all <- matrix(, nrow = intput_cluster_number, ncol = intput_cluster_number)    #FindMarkers결과값을 저장할 empty matrix

for(i in (intput_cluster_number-1):0){    #중첩된 loop문 FindMarkers에 숫자 넣어줌
  for(j in 0:(i-1)){
    if ((i-1) >= 0){
      bim <- FindMarkers(H_Epi.combined_qc, ident.1 = i,  ident.2 =j, only.pos = F, test.use = "bimod")
      bim_thresh <- subset(bim, (bim$p_val_adj < 0.01)&(abs(bim$avg_logFC) >= 1))
      num_deg <- dim(bim_thresh)[1]
      
      mtx_all[i+1,j+1] <- dim(bim_thresh)[1]    #FindMarkers결과값을 mtx_all에 순서에 맞춰서 넣어줌
    }
  }
}
mtx_all

mtx_all[is.na(mtx_all)] <- 0

df_mtx <- data.frame(mtx_all)
rownames(df_mtx) <- c(0:(intput_cluster_number-1))
colnames(df_mtx) <- c(0:(intput_cluster_number-1))

#install.packages("heatmaply")
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
