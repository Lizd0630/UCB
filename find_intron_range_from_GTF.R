
#### Find all of the introns range of gencode.v15.annotation.gtf
library(dplyr)

setwd("/mnt/data5/BGI/UCB/tangchao/IR/")

data <- read.table("/mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.87.gtf",sep="\t",header=F)
data <- data[data$V3 == "exon", ]

test <- as.data.frame(do.call(rbind, strsplit(as.character(data$V9), split = "[; ]")))


valueFind <- function(x, type = "gene_id"){
  for(i in 1:length(x)){
    if(sum(x == type) == 0){
      return(NA)
    }else{
      return(x[which(x == type)[1]+1])
    }
  }
}


gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
           transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
           exon_numble = as.vector(apply(test,1,function(x){valueFind(x,type = "exon_number")})),
           gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
           gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
           exon_id = as.vector(apply(test,1,function(x){valueFind(x,type = "havana_transcript")})), stringsAsFactors=F)


test2 <- data.frame(chr = as.character(data$V1), type = as.character(data$V3),
                    start = as.character(data$V4), end = as.character(data$V5),
                    strand = as.character(data$V7)) 

data <- cbind(test2, gtf_infor)
data <- data.frame(data, stringsAsFactors = F)



library(plyr)
data_split <- dlply(data, .variables = "transcript_id")

library(parallel)
intron_range <- mclapply(data_split, function(x) {
  if(nrow(x) == 1){
    return(NULL)
  } else {
    temp = data.frame(apply(x, 2, as.character), stringsAsFactors=F)

    test <- sort(c(temp$end, temp$start))[-c(1,length(c(temp$end, temp$start)))]

    result <- data.frame(chr = unique(temp$chr), start = test[seq(1,length(test),2)],
                         end = test[seq(2,length(test),2)], strand = unique(temp$strand),
                         gene_id = unique(temp$gene_id), transcript_id = unique(temp$transcript_id),
                         intron_numble = 1:length(seq(1,length(test),2)))

    return(result)
  }
}, mc.cores = 4)

## for loop
intron_range <- list()
for(i in 1:length(data_split)){
  print(paste(i, "of", length(data_split)))
  temp = data.frame(apply(data_split[[i]], 2, as.character), stringsAsFactors=F)

  if(nrow(temp) == 1){
    intron_range[i] <- NULL
  } else {
    test <- sort(c(temp$end, temp$start))[-c(1,length(c(temp$end, temp$start)))]

    intron_range[[i]] <- data.frame(chr = unique(temp$chr), start = test[seq(1,length(test),2)],
                                    end = test[seq(2,length(test),2)], strand = unique(temp$strand),
                                    gene_id = unique(temp$gene_id), transcript_id = unique(temp$transcript_id),
                                    intron_numble = 1:length(seq(1,length(test),2)))
  }
  rm(list =c("temp","test"))
}



## strand - need to reverse
for(i in 1:length(intron_range)){
  print(paste(i,"of",length(intron_range)))
  if(!is.null(unique(intron_range[[i]]$strand))){
    if(unique(intron_range[[i]]$strand) == "-"){
      intron_range[[i]]$intron_numble <- rev(intron_range[[i]]$intron_numble)
    }
  }
}

library(tidyr)

intron_range_bind <- bind_rows(intron_range)

#write.table(intron_range_bind, "./intron_range_bind.txt", sep = "\t", row.names = F, quote = F)

intron_range_bind_allInfo <- merge(x = intron_range_bind, y = data[,c(7,9,10)], by.x = "transcript_id", by.y = "transcript_id", all = F)
intron_range_bind_allInfo <- intron_range_bind_allInfo[!duplicated(intron_range_bind_allInfo),]
intron_range_bind_allInfo <- intron_range_bind_allInfo[intron_range_bind_allInfo$start < intron_range_bind_allInfo$end, ]# Remov the fusion gene
write.table(intron_range_bind_allInfo, "/mnt/data5/BGI/UCB/tangchao/IR/intron_range_bind_allInfo_nofusion.txt", sep = "\t", row.names = F, quote = F)

intron_range_bind_allInfo[!duplicated(intron_range_bind_allInfo[,c("chr","start","end")]),c("chr","start","end")] -> all_introns

all_introns <- all_introns[order(all_introns$chr,all_introns$start,all_introns$end),]

write.table(all_introns, "/mnt/data5/BGI/UCB/tangchao/IR/all_introns.txt", sep = "\t", row.names = F, quote = F, col.names=F)























