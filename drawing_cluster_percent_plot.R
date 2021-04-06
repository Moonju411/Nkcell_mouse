#Drawing cluster per barplot by group
library(ggplot2)
library(reshape2)

#Drawing cluster per barplot
library(ggplot2)
library(reshape2)

freq_table <- prop.table(x = table(Epi.integrated@active.ident, Epi.integrated@meta.data[, "group_id"]),
                         margin = 2)

freq_table2 <- melt(freq_table)
colnames(freq_table2)  <- c('cluster_id','group_id','freq')
freq_table2$cluster_id <- as.character(freq_table2$cluster_id)
freq_table2$freq       <- round(freq_table2$freq*100, digits = 1)

ggplot(freq_table2, aes(fill = cluster_id, y = freq, x = group_id, label = freq, order = as.numeric(cluster_id))) + 
  geom_bar(stat = "identity" )+
  geom_text(size = 3, position = position_stack(vjust = 0.5))


ggplot(freq_table2, aes(fill = cluster_id, y = freq, x = group_id, label = freq, order = as.numeric(freq)))+
  geom_bar(stat = "identity")
