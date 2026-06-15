library(tidyverse)

table <- read.table("short-variant_tr.csv", header=T, sep=",")
table2 <- table[, colnames(table) != "subclonal"]
write.table(table2, "short-variant_tr.cut.csv", sep=",", row.names = F, quote = F, col.names = T)
table3 <- table[, "subclonal"] %>% as.data.frame() %>% dplyr::rename(subclonal=".")
write.table(table3, "short-variant_tr.sub.csv", sep=",", row.names = F, quote = F, col.names = T)

