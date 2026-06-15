library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

CC      <- rep(NA, 16) # NAだけのベクトル
CC[16]   <- "character" # 1つ目だけ"character"に書き換える

#del_short の処理
try(del_short <- read.table("short-variant_tr.change2.combined2.del2.short.csv", header=F, sep=",", colClasses=CC))

if(exists("del_short")){
  
  del_short$V16 <- sub("TRUE", "T", del_short$V16)
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
  write.table(del_short2, "short-variant_tr.change2.combined2.del2.short.join.csv",sep=",", row.names = F, quote = F, col.names = F)
  
} else {
  
  del_short2 <- NULL
  write.table(del_short2, "short-variant_tr.change2.combined2.del2.short.join.csv",sep=",", row.names = F, quote = F, col.names = F)
  
}



#del_large の処理
try(del_large <- read.table("short-variant_tr.change2.combined2.del2.large2.csv", header=F, sep=",", colClasses=CC))

if(exists("del_large")){
  
  del_large$V16 <- sub("TRUE", "T", del_large$V16)
  del_large_chr <- strsplit(del_large[,1], ":")
  del_large_chr2 <-  t(as.data.frame(del_large_chr))
  row.names(del_large_chr2) <- NULL

  del_large_base <- NULL
  chr_number <- NULL
  base_number1 <- NULL
  base_number2 <- NULL
  base <- NULL
  base1 <- NULL
  base2 <- NULL
  row_length_del_large <- as.numeric(count(del_large))
  
  for(i in 1:row_length_del_large){
  
    chr_number <- del_large_chr2[i,1]
    base_number1 <- as.numeric(del_large_chr2[i,2])
    base_number2 <- as.numeric(del_large_chr2[i,2]) + as.numeric(del_large[i,17])
    base1 <- eval(parse(text=paste0("Hsapiens$",chr_number)))[base_number1:base_number2] %>% as.character()
    base2 <- as.data.frame(eval(parse(text=paste0("Hsapiens$",chr_number)))[base_number1])[1,1] %>% as.character()
    base <- cbind(base1, base2)
    del_large_base <- rbind(del_large_base, base)
  }

  row.names(del_large_base) <- NULL
  del_large2 <- cbind(del_large[,c(1:15)], del_large_base)
  write.table(del_large2, "short-variant_tr.change2.combined2.del2.large.join.csv",sep=",", row.names = F, quote = F, col.names = F)
  
} else {
  
  del_large2 <- NULL
  write.table(del_large2, "short-variant_tr.change2.combined2.del2.large.join.csv",sep=",", row.names = F, quote = F, col.names = F)
  
}


#splice site の処理
try(del_splice <- read.table("short-variant_tr.change2.combined2.splice5.csv", header=F, sep=",", colClasses=CC))

if(exists("del_splice")){
  
  del_splice$V16 <- sub("TRUE", "T", del_splice$V16)
  del_splice_chr <- strsplit(del_splice[,1], ":")
  del_splice_chr2 <-  t(as.data.frame(del_splice_chr))
  row.names(del_splice_chr2) <- NULL

  del_splice_base <- NULL
  chr_number <- NULL
  base_number1 <- NULL
  base_number2 <- NULL
  base <- NULL
  base1 <- NULL
  base2 <- NULL
  row_length_del_splice <- as.numeric(count(del_splice))
  
  for(i in 1:row_length_del_splice){
  
    chr_number <- del_splice_chr2[i,1]
    base_number1 <- as.numeric(del_splice_chr2[i,2])
    base_number2 <- as.numeric(del_splice_chr2[i,2]) + as.numeric(del_splice[i,17])
    base <- eval(parse(text=paste0("Hsapiens$",chr_number)))[base_number1:base_number2] %>% as.character()
    del_splice_base <- rbind(del_splice_base, base)
  }

  row.names(del_splice_base) <- NULL
  del_splice2 <- cbind(del_splice[,c(1:16)], del_splice_base)
  write.table(del_splice2, "short-variant_tr.change2.combined2.splice5.join.csv",sep=",", row.names = F, quote = F, col.names = F)
  
} else {
  
  del_splice2 <- NULL
  write.table(del_splice2, "short-variant_tr.change2.combined2.splice5.join.csv",sep=",", row.names = F, quote = F, col.names = F)
  
}




#insの処理
try(ins <- read.table("short-variant_tr.change2.combined2.ins2.csv", header=F, sep=",", colClasses=CC))

if(exists("ins") && !all(is.na(ins$V16) | ins$V16 == "")){

  ins$V16 <- sub("TRUE", "T", ins$V16)
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
  write.table(ins2, "short-variant_tr.change2.combined2.ins2.join.csv",sep=",", row.names = F, quote = F, col.names = F)
  
} else if (exists("ins") && all(is.na(ins$V16) | ins$V16 == "")) {
  
  generate_random_string <- function(length) {
    letters <- c("A", "T", "G", "C")
    return(paste0(sample(letters, size = length, replace = TRUE), collapse = ""))
  }

  extracted_numbers <- as.numeric(gsub("^.*ins([0-9]+).*", "\\1", ins[,5]))
  random_strings <- sapply(extracted_numbers, generate_random_string)
  ins$V16 <- random_strings

  ins$V16 <- sub("TRUE", "T", ins$V16)
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
  write.table(ins2, "short-variant_tr.change2.combined2.ins2.join.csv",sep=",", row.names = F, quote = F, col.names = F)

} else {

  ins2 <- NULL
  write.table(ins2, "short-variant_tr.change2.combined2.ins2.join.csv",sep=",", row.names = F, quote = F, col.names = F)
  
}
  