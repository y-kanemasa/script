library(tidyverse)
library(vcfR)

vcf <- read.vcfR("annotation.snpeff.sort.annotate.vcf")
vcf_df <- vcfR2tidy(vcf)
write.table(vcf_df$fix, "annotation.snpeff.sort.annotate.tsv", sep = "\t", row.names = F, na = ".", quote=F)