library(ggplot2)
library(dplyr)
draw_condition_correlation_plot_ver3 <- function(input_seuratobj, cluster= NULL, top_n= 10 ) {
  DefaultAssay(input_seuratobj) <- "RNA"
  subset_cluster <- subset(input_seuratobj, idents = cluster) # 보고 싶은 클러스터로 subset
  
  Idents(subset_cluster) <- subset_cluster@meta.data$orig.ident
  
  diff_condition_marker <<- FindMarkers(subset_cluster, ident.1 = levels(Idents(subset_cluster))[1], ident.2 = levels(Idents(subset_cluster))[2])
  H_enriched_top10 <- diff_condition_marker %>% top_n(n = 10, wt = avg_logFC) # 표시할 유전자 고르기
  L_enriched_top10 <- diff_condition_marker %>% top_n(n = 10, wt = -avg_logFC) # 표시할 유전자 고르기
  
  subset_cluster <- as.data.frame(log1p(AverageExpression(subset_cluster, verbose = FALSE)$RNA))
  subset_cluster$gene <- rownames(subset_cluster)
  
  df_mark_H = subset(subset_cluster, subset_cluster[['gene']] %in% rownames(H_enriched_top10))
  df_mark_L = subset(subset_cluster, subset_cluster[['gene']] %in% rownames(L_enriched_top10))
  
  #여기서부터 플롯 그리기
  p1 <- ggplot(subset_cluster, aes_string(names(subset_cluster)[1], names(subset_cluster)[2])) + geom_point(color="gray")
  p1 <- LabelPoints(plot = p1, points = rownames(H_enriched_top10), repel = TRUE)
  p1 <- LabelPoints(plot = p1, points = rownames(L_enriched_top10), repel = TRUE)
  
  p1 <- p1 + geom_point(data= df_mark_H, aes_string(names(df_mark_H)[1], names(df_mark_H)[2]), color= "darkblue")
  p1 <- p1 + geom_point(data= df_mark_L, aes_string(names(df_mark_L)[1], names(df_mark_L)[2]), color= "darkred")
  p1 <- p1 + xlab(paste('log1p(Avg.expression(' , names(df_mark_H)[1] , '))')) + ylab(paste('log1p(Avg.expression(' , names(df_mark_H)[2] , '))')) + theme(legend.position="none")
  
  if (length(cluster) > 1){
    p_out <<- p1 + ggtitle(paste("cluster", paste( unlist(cluster), collapse=' and ')))
  } else {
    p_out <<- p1 + ggtitle(paste('cluster', cluster))
  }
  p_out
}

plist <- list()  #플롯들을 저장할 리스트
markerlist <- list() #마커들을 저장할 리스트
for(i in 0:(length(levels(Idents(BM_integrated)))-1)){
  draw_condition_correlation_plot_ver3(BM_integrated, cluster = i)
  plist[[i+1]] <- p_out
  markerlist[[i+1]] <- diff_condition_marker
}
pdf("test.pdf") #플롯들 pdf로 한번에 저장
plist
dev.off() 

cluster_19 <- subset(BM_integrated, idents = "19")
Idents(cluster_19) <- cluster_19@meta.data$orig.ident
cluster_19_H <- subset(cluster_19, idents = "HBM")
cluster_19_L <- subset(cluster_19, idents = "LBM")
cluster_19_H.avg <- AverageExpression(cluster_19_H)
cluster_19_L.avg <- AverageExpression(cluster_19_L)
cluster_19_merged <- merge(cluster_19_H.avg$RNA, cluster_19_L.avg$RNA)
cluster_19_merged
