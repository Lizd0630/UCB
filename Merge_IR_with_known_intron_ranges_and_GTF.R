## Merge IR with known intron ranges and GTF
library(GenomicRanges)
library(dplyr)
library(plyr)
library(pryr)

#### Retained introns merge found introns range ============================================================================================================


gtf <- read.table("/mnt/data5/BGI/UCB/tangchao/IR/intron_range_bind_allInfo_nofusion.txt",head=T,sep="\t", stringsAsFactors = F)

gtf <- gtf[order(gtf$chr, gtf$start, gtf$end),]
sum(gtf$start > gtf$end)
# [1] 180
gtf <- gtf[gtf$start < gtf$end, ]

intron_ca <- dlply(gtf, "chr")

intron_c <- gtf[,2:4]

mode(intron_c)

intron_cna <- split(IRanges(intron_c[,2], intron_c[,3]), intron_c[,1])


load(file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/IR_irl_tu.RData")

all_ir_id <- row.names(IR_irl_tu)

IR <- data.frame(chr = do.call(rbind,strsplit(all_ir_id, split = "[:-]"))[,1], 
  start = do.call(rbind,strsplit(all_ir_id, split = "[:-]"))[,2],
  end = do.call(rbind,strsplit(all_ir_id, split = "[:-]"))[,3], event_id = all_ir_id, stringsAsFactors = F)

IR <- IR[IR$chr %in% c(1:22,"X","Y"),]
IR <- IR[order(IR$chr,IR$start,IR$end),]## It's important!

unique(IR$chr)
chr <- unique(IR$chr)

IR_type <- dlply(IR, .variables = "chr")

for(j in 1:length(chr)){
  print(chr[j])
  eval(parse(text = paste("intron_cna$'",chr[j],"'->intron_cn",sep="")))
  temp <- IR[IR$chr == chr[j],]
  if(nrow(temp) > 1){
    test <- as.matrix(temp[,1:3])
  }else{
    test <- t(as.matrix(temp[,1:3]))
  }
  
  mode(test) <- "numeric"
  for(i in 1:nrow(test)){
    ASrange <- IRanges(test[i,2]-1, test[i,3]+1)
    ## because the positions of SJ are intron sets, so we  +- 1 to match exon sets
    ASlen <- diff(c(test[i,2], test[i,3]))
    ov_within<-countOverlaps(intron_cn, ASrange, type="within")
    ov_start<-countOverlaps(intron_cn, ASrange, type="start")
    ov_end<-countOverlaps(intron_cn, ASrange, type="end")
    ov_equal<-countOverlaps(intron_cn, ASrange, type="equal")
    ov_any<-countOverlaps(intron_cn,ASrange, type="any")
    if (sum(ov_equal)>0) {
      AStype <- "intron_retention"
      
    } else {
      AStype <- "Others" 
    } 

    IR_type[[j]][i,"ASrange"] <- paste("chr", chr[j], ":", test[i,2],"-",test[i,3],sep="")
    IR_type[[j]][i,"AStype"] <- AStype
    IR_type[[j]][i,"OverlapCount"] <- paste(sum(ov_within),"|",sum(ov_start),"|",sum(ov_end),"|",sum(ov_equal),"|",sum(ov_any),sep="")
    IR_type[[j]][i,"ASlen"] <- ASlen
    print(paste(i," chr", chr[j], ":", test[i,2],"-",test[i,3]," ", AStype," - ",ASlen,sep=""))
  }
}
IR_type <- do.call(rbind, IR_type)
head(IR_type)

table(IR_type$AStype)
# intron_retention           Others
#             4120             2586


IR_exactly <- IR_type[IR_type$AStype == "intron_retention", ]
IR_novel <- IR_type[IR_type$AStype == "Others", ]


#### orignal intron's gene finding (merge with known introns) ============================================================================================================================


## all chr


chr <- unique(IR_exactly$chr)

IR_gtf_chr <- list()

system.time(for(j in 1:length(chr)){
  print(paste("chr",chr[j]))
  eval(parse(text = paste("intron_cna$'",chr[j],"'->intron_chr",sep="")))
  temp <- IR_exactly[IR_exactly$chr == chr[j],]
  
  eval(parse(text = paste("intron_ca$'",chr[j],"'->intron_ca_chr",sep="")))
  
  if(nrow(temp) > 1){
    test <- as.matrix(temp[,1:3])
  }else{
    test <- t(as.matrix(temp[,1:3]))
  }
  colnames(test) <- paste("V", 1:3, sep = "")
  mode(test) <- "numeric"
  
  ir_gtf_chr <- list()
  
  for(i in 1:nrow(test)){
    #print(i)
    ASrange <- IRanges(test[i,2]-1, test[i,3]+1)
    ## because the positions of SJ are intron sets, so we  +- 1 to match exon sets
    
    hit_any <- queryHits(findOverlaps(intron_chr, ASrange, type="equal"))
    
    if(length(hit_any) > 0){
      within_all <- data.frame(OverType = "equal", intron_ca_chr[hit_any,])
    }
    
    ir_gtf_chr[[i]] <- suppressWarnings(cbind(temp[i,], within_all))
    
  }
  
  IR_gtf_chr[[j]] <- suppressMessages(do.call(rbind,ir_gtf_chr))
})
#     user   system  elapsed 
# 4079.416   40.844 4119.096 

system.time(IR_GTF <- suppressMessages(do.call(rbind,IR_gtf_chr)))

dim(IR_GTF)
# [1] 16027    18
length(unique(IR_GTF$gene_id))
# [1] 2320

save(IR_GTF, file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/IR_GTF_merge/IR_known_introns_match.RData")



#### Overlap with GTF =======================================================================================================================================



gtf <- read.table("/mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.87.gtf",sep="\t",header=F)

gene <- subset(gtf,gtf$V3=="gene")
gene <- gene[!duplicated(gene[,c(1,4,5)]),]
gene <- gene[gene$V1 %in% unique(IR_novel$chr), ]
gene <- gene[order(gene$V1,gene$V4,gene$V5),]

gene_ca <- dlply(gene, "V1")
gene_c <- gene[,c(1,4,5)]
mode(gene_c)
row.names(gene_c) <- 1:nrow(gene_c)
gene_cn <- split(IRanges(gene_c[,2], gene_c[,3]), gene_c[,1])


IR_novel$start <- as.integer(IR_novel$start)
IR_novel$end <- as.integer(IR_novel$end)
IR_novel <- IR_novel[order(IR_novel$chr, IR_novel$start, IR_novel$end),]
row.names(IR_novel) <- 1:nrow(IR_novel)
IR_ca <- dlply(IR_novel, "chr")

IR_novel[,1:3] -> IR_c
IR_c$start <- as.integer(IR_c$start)
IR_c$end <- as.integer(IR_c$end)
IR_c <- IR_c[order(IR_c$chr, IR_c$start, IR_c$end),]
mode(IR_c)
row.names(IR_c) <- 1:nrow(IR_c)
IR_cn <- split(IRanges(IR_c[,2], IR_c[,3]), IR_c[,1])





#### orignal gene finding ####


## all chr


chr <- unique(IR_novel$chr)

IR_gtf_chr <- list()

system.time(for(j in 1:length(chr)){
  print(paste("chr",chr[j]))
  eval(parse(text = paste("IR_cn$'",chr[j],"'->IR_chr",sep="")))

  temp <- gene[gene$V1 == chr[j],]
  
  eval(parse(text = paste("IR_ca$'",chr[j],"'->IR_ca_chr",sep="")))
  
  if(nrow(temp) > 1){
    test <- as.matrix(temp[,4:5])
  }else{
    test <- t(as.matrix(temp[,4:5]))
  }
  colnames(test) <- paste("V", 1:2, sep = "")
  mode(test) <- "numeric"
  
  ir_gtf_chr <- list()
  
  for(i in 1:nrow(test)){
    #print(i)
    Gene_range <- IRanges(test[i,1], test[i,2])
    ## because the positions of SJ are intron sets, so we  +- 1 to match exon sets
    
    hit_within <- queryHits(findOverlaps(IR_chr, Gene_range, type="within")) ## We must change the position of "ASrange" and "gene_chr".
    
    if(length(hit_within) > 0){
      within_all <- data.frame(OverType = "within", IR_ca_chr[hit_within,])
      ir_gtf_chr[[i]] <- suppressWarnings(cbind(temp[i,], within_all))
    }    
  }
  
  IR_gtf_chr[[j]] <- suppressMessages(do.call(rbind, ir_gtf_chr))
})
#     user   system  elapsed 
# 4079.416   40.844 4119.096 

system.time(IR_GTF_novel <- do.call(rbind,IR_gtf_chr))

dim(IR_GTF_novel)
# [1] 2902   18
dim(IR_novel)
# [1] 2586    8

length(table(IR_GTF_novel$event_id))
# [1] 2426

## 2586 - 2426 = 160 intergenic IR

sum(table(IR_GTF_novel$event_id) == 1)
# [1] 1991
## 1991 of 2426 IR just match with only one gene.
## 2426 - 1991 = 435 IR match more than one gene.


save(IR_GTF_novel, file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/IR_GTF_merge/IR_GTF_match.RData")


table(table(IR_GTF_novel$event_id))
#    1    2    3    4
# 1991  399   31    5
IR_GTF_novel_match_multigene <- IR_GTF_novel[IR_GTF_novel$event_id %in% names(table(IR_GTF_novel$event_id)[table(IR_GTF_novel$event_id)>1]),]
dim(IR_GTF_novel_match_multigene)
# [1] 911  18


















