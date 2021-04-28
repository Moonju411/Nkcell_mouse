library(ggplot2)
library(reshape2)
library(Seurat)

#Drawing stacked-barplot per cluster on each group ( H-Epi vs L-Epi )
draw_stacked_bar_per_cluster <- function(input_seuratobj ){
  #first, read seurat metadata and generate dataframe for generate plot
  freq_table <- prop.table(x = table(input_seuratobj@active.ident, input_seuratobj@meta.data[, "group_id"]), margin = 2)
  freq_table <- melt(freq_table)
  colnames(freq_table)  <- c('cluster_id','group_id','freq')
  freq_table$cluster_id <- as.character(freq_table$cluster_id)
  freq_table$freq <- round(freq_table$freq*100, digits = 1)
  
  draw_stacked_bar_per_cluster_freq_table <<- freq_table
  
  ggplot(freq_table, aes(fill = cluster_id, y = freq, x = group_id, label = freq)) + #label제거할 때는 label=freq을 제거하시면 됩니다. 
    geom_bar(stat = "identity" ) + 
    geom_text(size = 3, position = position_stack(vjust = 0.5))
}
#For example,
draw_stacked_bar_per_cluster(Epi.integrated.F) # freq table is generated as draw_stacked_bar_per_cluster_freq_table


##Drawing stacked-barplot per group on each cluster ( H-Epi vs L-Epi )
draw_stacked_bar_per_group <- function(input_seuratobj ){
  freq_table <- prop.table(x = table(input_seuratobj@active.ident, input_seuratobj@meta.data[, "group_id"]), margin = 1)
  freq_table <- melt(freq_table)
  colnames(freq_table)  <- c('cluster_id','group_id','freq')
  freq_table$cluster_id <- as.character(freq_table$cluster_id)
  freq_table$freq       <- round(freq_table$freq*100, digits = 1)
  
  draw_stacked_bar_per_group_freq_table <<- freq_table
  
  ggplot(freq_table, aes(fill = group_id, y = freq, x =reorder(cluster_id, sort(as.numeric(cluster_id))), label = freq, order = as.numeric(cluster_id))) + #label제거할 때는 label=freq을 제거하시면 됩니다. 
    geom_bar(stat = "identity" )+
    geom_text(size = 3, position = position_stack(vjust = 0.5))+
    ggtitle("stacked_bar_per_group")+
    xlab("cluster id")
}
#For example,
draw_stacked_bar_per_group(Epi.integrated.F) # freq table is generated as draw_stacked_bar_per_group_freq_table
