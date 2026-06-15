library(tidyverse)
library(biomaRt)
options(timeout = 50000)


ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version=108)
gene_list <- read.delim("gene_list.txt", header=F, sep='\t') %>% dplyr::distinct() %>% as.matrix() %>% as.character()
results = getBM(attributes=c('hgnc_symbol', 'refseq_mrna', 'strand', 'ensembl_transcript_id', 'transcript_mane_select'), filters = 'refseq_mrna', values = gene_list, mart = ensembl)
results <- results %>% group_by(refseq_mrna) %>% filter(n() == 1 | # 重複していない行はそのまま残す
    (n() > 1 & transcript_mane_select != "") # 重複している場合、transcript_mane_selectが空でない行を残す
  ) %>%
  ungroup()

write.table(results, "annotation_NM.tsv", row.names=F, quote=F, sep="\t")
