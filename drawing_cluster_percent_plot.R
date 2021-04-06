#Drawing cluster per barplot by group
library(ggplot2)
library(reshape2)

freq_table <- prop.table(x = table(Epi.integrated@active.ident, Epi.integrated@meta.data[, "group_id"]),
                         margin = 2)

freq_table2 <- melt(freq_table)
freq_table2$Var1 <- as.character(freq_table2$Var1)
colnames(freq_table2) <- c('cluster_id','group_id','freq')

ggplot(freq_table2, aes(fill = cluster_id, y = freq, x = group_id)) + geom_bar(stat = "identity")
