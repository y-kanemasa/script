library(tidyverse)
library(vcfR)

vcf <- read.vcfR("shortVariants3.snpeff.sort.annotate.vcf")
vcf_df <- vcfR2tidy(vcf)
write.table(vcf_df$fix, "shortVariants3.snpeff.sort.annotate.tsv", sep = "\t", row.names = F, na = ".", quote=F)