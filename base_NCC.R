library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)


#del_short の処理
try(del_short <- read.table("Book4.del2.csv", header=F, sep=","))

if(exists("del_short")){
  
  del_short_chr <- strsplit(del_short[,1], ":")
  del_short_chr2 <-  t(as.data.frame(del_short_chr))
  row.names(del_short_chr2) <- NULL

  del_short_base <- NULL
  chr_number <- NULL
  base_number <- NULL
  base <- NULL
  row_length_del_short <- as.numeric(count(del_short))
  
  for(i in 1:row_length_del_short){
  
    chr_number <- del_short_chr2[i,1]
    base_number <- as.numeric(del_short_chr2[i,2])
    base <- as.data.frame(eval(parse(text=paste0("Hsapiens$",chr_number)))[base_number])[1,1] %>% as.character()
    del_short_base <- rbind(del_short_base, base)
  }

  row.names(del_short_base) <- NULL
  del_short2 <- cbind(del_short, del_short_base)
  write.table(del_short2, "Book4.del2.join.csv",sep=",", row.names = F, quote = F, col.names = F)
  
} else {
  
  del_short2 <- NULL
  write.table(del_short2, "Book4.del2.join.csv",sep=",", row.names = F, quote = F, col.names = F)
  
}



#delins の処理
try(delins <- read.table("Book4.delins2.csv", header=F, sep=","))

if(exists("delins")){
  
  delins_chr <- strsplit(delins[,1], ":")
  delins_chr2 <-  t(as.data.frame(delins_chr))
  row.names(delins_chr2) <- NULL
  
  delins_del <- gsub("[a-zA-Z.]", "", delins[,4])
  delins_del_list <- strsplit(delins_del, "_") %>% as.data.frame()
  delins_del1 <- delins_del_list[1,1] %>% as.numeric()
  delins_del2 <- delins_del_list[2,1] %>% as.numeric()
  delins_del3 <- delins_del2-delins_del1

  delins_base <- NULL
  chr_number <- NULL
  base_number1 <- NULL
  base_number2 <- NULL
  base <- NULL
  row_length_delins <- as.numeric(count(delins))
  
  for(i in 1:row_length_delins){
  
    chr_number <- delins_chr2[i,1]
    base_number1 <- as.numeric(delins_chr2[i,2])
    base_number2 <- as.numeric(delins_chr2[i,2]) + delins_del3
    base <- eval(parse(text=paste0("Hsapiens$",chr_number)))[base_number1:base_number2] %>% as.character()
    delins_base <- rbind(delins_base, base)
  }

  row.names(delins_base) <- NULL
  delins2 <- cbind(delins, delins_base)
  write.table(delins2, "Book4.delins2.join.csv",sep=",", row.names = F, quote = F, col.names = F)
  
} else {
  
  delins2 <- NULL
  write.table(delins2, "Book4.delins2.join.csv",sep=",", row.names = F, quote = F, col.names = F)
  
}





#insの処理
try(ins <- read.table("Book4.ins2.csv", header=F, sep=","))

if(exists("ins")){

  ins_chr <- strsplit(ins[,1], ":")
  ins_chr2 <-  t(as.data.frame(ins_chr))
  row.names(ins_chr2) <- NULL

  ins_base <- NULL
  chr_number <- NULL
  base_number <- NULL
  base <- NULL
  row_length_ins <- as.numeric(count(ins))
  
  for(i in 1:row_length_ins){
  
    chr_number <- ins_chr2[i,1]
    base_number <- as.numeric(ins_chr2[i,2])
    base <- as.data.frame(eval(parse(text=paste0("Hsapiens$",chr_number)))[base_number])[1,1] %>% as.character()
    ins_base <- rbind(ins_base, base)
  }

  row.names(ins_base) <- NULL
  ins2 <- cbind(ins, ins_base)
  write.table(ins2, "Book4.ins2.join.csv",sep=",", row.names = F, quote = F, col.names = F)
  
} else {
  
  ins2 <- NULL
  write.table(ins2, "Book4.ins2.join.csv",sep=",", row.names = F, quote = F, col.names = F)
  
}
  