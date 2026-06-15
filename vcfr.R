library(tidyverse)
library(vcfR)

vcf <- read.vcfR("short-variant_tr.change2.combined6.snpeff.sort.annotate.vcf")
vcf_df <- vcfR2tidy(vcf)
write.table(vcf_df$fix, "short-variant_tr.change2.combined6.snpeff.sort.annotate.tsv", sep = "\t", row.names = F, na = ".", quote=F)