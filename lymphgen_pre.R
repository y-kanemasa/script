library(tidyverse)
library(org.Hs.eg.db)

# 遺伝子リスト読み込み
genes <- read.table("genes319.txt", stringsAsFactors = FALSE)$V1

# シンボル修正
genes[genes == "FAM46C"]   <- "TENT5C"
genes[genes == "HIST1H1E"] <- "H1-4"

# ENTREZ ID へ変換
mapped <- AnnotationDbi::select(org.Hs.eg.db,
                                keys = genes,
                                keytype = "SYMBOL",
                                columns = c("ENTREZID"))

# 保存
write.table(mapped[2], "genes319_entrez.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)