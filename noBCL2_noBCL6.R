library(tidyverse)

# 今いるフォルダ名を取得
current_dir <- basename(getwd())

# "_" で分割して最後の要素を取得
sample_name <- tail(strsplit(current_dir, "-")[[1]], 1)

sample <- data.frame(
  Sample.ID     = sample_name[1],
  Copy.Number   = 0L,
  BCL2.transloc = 0L,
  BCL6.transloc = 0L,
  check.names   = FALSE
)


# 保存
write.table(sample, "samp.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)