library(tidyverse)
library(org.Hs.eg.db)

# 遺伝子リスト読み込み
shortVariants <- read.table("shortVariants3.snpeff.sort.annotate2.hg19_3.csv", stringsAsFactors = FALSE, header=T, sep=",")
shortVariants <- shortVariants %>% filter(!is.na(.[[2]]) & grepl("^MT", as.character(.[[2]])))
genes <- shortVariants$geneSymbol

# シンボル修正
genes[genes == "FAM46C"]   <- "TENT5C"
genes[genes == "HIST1H1E"] <- "H1-4"

# ENTREZ ID へ変換
mapped <- AnnotationDbi::select(org.Hs.eg.db,
                                keys = genes,
                                keytype = "SYMBOL",
                                columns = c("ENTREZID"))

shortVariants$Type <- NA  # 新しい列を初期化

# 1. TRUNC
shortVariants$Type[grepl("frameshift_variant|stop_gained|splice", shortVariants$calculatedEffects)] <- "TRUNC"

# 2. MUTATION
shortVariants$Type[grepl("missense_variant", shortVariants$calculatedEffects)] <- "MUTATION"

# 3. MYD88 L265P 特例
shortVariants$Type[shortVariants$geneSymbol == "MYD88" & shortVariants$aminoAcidsChange == "p.L252P"] <- "L265P"

shortVariants2 <- cbind(shortVariants, mapped)
colnames(shortVariants2)[colnames(shortVariants2) == "ENTREZID"] <- "ENTREZ.ID"
colnames(shortVariants2)[colnames(shortVariants2) == "POS"] <- "Location"

# 今いるフォルダ名を取得
current_dir <- basename(getwd())

# "_" で分割して最後の要素を取得
sample_name <- tail(strsplit(current_dir, "-")[[1]], 1)

# データフレームに列を追加（全行に同じ値）
shortVariants2$Sample <- sample_name

shortVariants3 <- shortVariants2[, c("Sample", "ENTREZ.ID", "Type", "Location")]

shortVariants3 <- shortVariants3 %>% unique()

# 保存
write.table(shortVariants3, "Mutation_flat.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)