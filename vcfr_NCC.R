library(tidyverse)
library(vcfR)

vcf <- read.vcfR("Book9.snpeff.sort.annotate.vcf")
vcf_df <- vcfR2tidy(vcf)
write.table(vcf_df$fix, "Book9.snpeff.sort.annotate.tsv", sep = "\t", row.names = F, na = ".", quote=F)