## pipeline of DSU for single cell

#### 0. basic settings ---------------------------------------------------------------------------------------------------------------------------
require(parallel)
library(dplyr)
library(plyr)
library(ggplot2)
library(cowplot)
library(GenomicRanges)
library(reshape2)
library(vioplot)

depar <- par()

## path and name of SJ merged RData

## path of figure/file output
pffo <- file.path("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/")

## path of RData output
pfro <- file.path("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/RData/")

## path of gtf
pfgtf <- file.path("/mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.87.gtf")


#### 1. prepare data ---------------------------------------------------------------------------------------------------------------------------



gtf <- read.table(pfgtf,head=F,sep="\t")
exon <- subset(gtf,gtf$V3=="exon")
exon <- exon[!duplicated(exon[,c(1,4,5)]),]
exon_ca <- dlply(exon, "V1")
exon_c <- exon[,c(1,4,5)]
mode(exon_c)
exon_cn <- split(IRanges(exon_c[,2], exon_c[,3]), exon_c[,1])

pfsj <- file.path("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te.RData")
load(pfsj)

te_table <- as.data.frame(te)

row.names(te_table) <- te_table$name
te_table <- te_table[,-1]

colnames(te_table) <- substr(colnames(te_table), 1, 10)

require(data.table)  # v1.6.6
require(gdata) 
f_dowle3 = function(DT) {
  # either of the following for loops
  
  # by name :
  for (j in names(DT))
    set(DT,which(is.na(DT[[j]])),j,0)
  
  # or by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    set(DT,which(is.na(DT[[j]])),j,0)
}

f_dowle3(te_table)

save(te_table, file = "/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te(no_NA).RData")

SJ_tu <- te_table


## We only care about chr1-22 and X Y
chr_tu <- do.call(rbind, strsplit(rownames(SJ_tu), split = ":"))[,1] %in% c(1:22,"X","Y")
sum(chr_tu)
# [1] 2053238

SJ_tu <- SJ_tu[chr_tu,]
dim(SJ_tu)
# [1] 2053238    3574


dim(SJ_tu)
   

SE_result <- list()
A3SS_result <- list()
A5SS_result <- list()
MXE_result <- list()


for(bl1 in 481:ncol(SJ_tu)){
print(paste(bl1,"of",ncol(SJ_tu)))
na.omit(subset.data.frame(SJ_tu, select=bl1)) -> count_use

subset.data.frame(count_use, subset = count_use[,1]>1) -> count_use


#### 1.find splicing junction ----------------------------------------------------------------------------------------------------------
print(paste("Step1: find splicing junction"))

junc=count_use

junc.names=do.call(rbind,strsplit(rownames(junc),split="[:-]"))

colnames(junc.names) <- c("chr","start","end")
rownames(junc.names) <- rownames(junc)
junc.names=data.frame(junc.names,stringsAsFactors=F)
junc.names$start=as.integer(junc.names$start)
junc.names$end=as.integer(junc.names$end)

junc.names=junc.names[order(junc.names$chr,junc.names$start,junc.names$end),]
junc.names$names=rownames(junc.names)


junc.as=do.call(c,mclapply(unique(junc.names$chr),function(chr) {
  junc.chr=junc.names[junc.names$chr==chr,]
  same.start=dlply(junc.chr,c("start"),function(x) x)
  same.end=dlply(junc.chr,c("end"),function(x) {if(nrow(x)>1) {return(x)} else {return(NULL)}})
  #print(chr)
  return(c(same.start,same.end[!sapply(same.end,is.null)]))
},mc.cores=10))

juncs <- junc.as

junc.as <- junc.as[sapply(junc.as,nrow)>1]
names(junc.as) <- sapply(junc.as,function(x) paste(x$names,collapse="_"))

## junction just has two alternative starts/ends:
junc.as.2=junc.as[sapply(junc.as,nrow)==2]

#### 2.calculate PSI of all junction 2 ----------------------------------------------------------------------------------------------------------
print(paste("Step2: calculate PSI of all junction 2"))


psi_sj2=mclapply(junc.as.2,function(sjs){
  sjs.names=sjs$names[order(as.numeric(sjs$end)-as.numeric(sjs$start))]# find the splicing_in isoform to calculate the psi
  tab=as.matrix(t(junc[sjs.names,]))
  distc=max(c(diff(as.numeric(sjs$start)),diff(as.numeric(sjs$end))))
  rs<-apply(tab, 1, function(a) sum(a))
  tsrs<-tab/rs
  return(tsrs[,1])
},mc.cores=8)

psi.co<-t(do.call(cbind,psi_sj2))
colnames(psi.co) <- colnames(junc)



#### 3.find SE ----------------------------------------------------------------------------------------------------------
print(paste("Step3: find SE"))


junc.as2.df=do.call(rbind,junc.as.2)
junc.as.2.ot <- junc.as2.df[order(junc.as2.df$chr,junc.as2.df$start,junc.as2.df$end),]

temp <- do.call(list, mclapply(1:length(unique(junc.as.2.ot$chr)), function(i){
  
  ## Calculated separately by chr
  junc.chr=junc.as.2.ot[junc.as.2.ot$chr==unique(junc.as.2.ot$chr)[i],]
  if(length(names(table(junc.chr$start)[table(junc.chr$start) == 3])) == 0 ){
    print(paste("The chr",unique(junc.as.2.ot$chr)[i],"have no SE", sep = " "))
  }else{
    ## left:
    ## because the longest SJ has been counted twice in SE, so the same starts/ends will appear three times
    se_start_left <- names(table(junc.chr$start)[table(junc.chr$start) == 3])
    
    se_left_all <- junc.chr[junc.chr$start %in% se_start_left, ]
    
    ## because we have ordered the data table by chr, start and end, 
    ## so the order of same start and end in SE will appear in a certain order: 1,2,1,2
    se_start_left_yn <- as.vector(rep(NA,length(unique(se_left_all$start))))
    for(j in 1:length(unique(se_left_all$start))){
      se_left_all_j <- se_left_all[se_left_all$start == unique(se_left_all$start)[j],]
      se_start_left_yn[j] <- do.call(rbind, strsplit(rownames(se_left_all_j), split = "\\."))[1,2] == "1" &
        do.call(rbind, strsplit(rownames(se_left_all_j), split = "\\."))[2,2] == "2" &
        do.call(rbind, strsplit(rownames(se_left_all_j), split = "\\."))[3,2] == "1"
    }
    se_left_all <- se_left_all[rep(se_start_left_yn,rep(3,length(se_start_left_yn))),]
    
    se_left <- se_left_all[seq(2,nrow(se_left_all),3),]
    se_left <- se_left[do.call(rbind,strsplit(rownames(se_left), split = "\\."))[,2] == 2, ]
    rownames(se_left) <- do.call(rbind,strsplit(rownames(se_left), split = "\\."))[,1]
    ## right:
    se_end_right <- names(table(junc.chr$end)[table(junc.chr$end) == 3])
    
    se_right_all <- junc.chr[junc.chr$end %in% se_end_right, ]
    
    se_end_right_yn <- as.vector(rep(NA,length(unique(se_right_all$end))))
    for(j in 1:length(unique(se_right_all$end))){
      se_right_all_j <- se_right_all[se_right_all$end == unique(se_right_all$end)[j],]
      se_end_right_yn[j] <- do.call(rbind, strsplit(rownames(se_right_all_j), split = "\\."))[1,2] == "2" &
        do.call(rbind, strsplit(rownames(se_right_all_j), split = "\\."))[2,2] == "1" &
        do.call(rbind, strsplit(rownames(se_right_all_j), split = "\\."))[3,2] == "2"
    }
    se_right_all <- se_right_all[rep(se_end_right_yn,rep(3,length(se_end_right_yn))),]
    
    se_right <- se_right_all[seq(2,nrow(se_right_all),3),]
    
    se_right <- se_right[do.call(rbind,strsplit(rownames(se_right), split = "\\."))[,2] == 1, ]
    
    rownames(se_right) <- do.call(rbind,strsplit(rownames(se_right), split = "\\."))[,1]
    ## merge:
    se_left_tu <- se_left[se_left$names %in% se_right$names,]
    se_right_tu <- se_right[se_right$names %in% se_left$names,]
    
    se_raw_name <- paste(row.names(se_left_tu),row.names(se_right_tu),sep = "_")
    
    se_names <- do.call(rbind, strsplit(se_raw_name,split = "[-:_\\.]"))[,c(1,2,3,11,12)]
    if(length(se_names) == 5){
      se_name <- paste(se_names[1], paste(se_names[2],se_names[3],se_names[4],se_names[5], sep = "-"), sep = ":")
    }else{
      se_name <- paste(se_names[,1], paste(se_names[,2],se_names[,3],se_names[,4],se_names[,5], sep = "-"), sep = ":")
    }
    se <- data.frame(loci = se_name, left = rownames(se_left_tu), right = rownames(se_right_tu), row.names = se_left_tu$names)
    
  
    return(se)
    }
  }, mc.cores = 10))

SE <- na.omit(suppressWarnings(do.call(rbind, temp)))

loci <- do.call(rbind, strsplit(as.vector(SE$loci), split = "[-:]"))

SE$chr <- loci[,1]
SE$dist_left <- as.numeric(loci[,5]) - as.numeric(loci[,3])
SE$dist_right <- as.numeric(loci[,4]) - as.numeric(loci[,2])
SE$dist_insert <- as.numeric(loci[,4]) - as.numeric(loci[,3])

SE <- SE[SE$dist_insert > 0,]


#### 4.SE type ----------------------------------------------------------------------------------------------------------
print(paste("Step4: SE type"))

chr <- unique(SE$chr)
SE_type <- dlply(SE, .variables = "chr")

for(j in 1:length(chr)){
  #print(chr[j])
  eval(parse(text = paste("exon_cn$'",chr[j],"'->exon_chr",sep="")))
  temp <- SE[SE$chr == chr[j],]
  if(nrow(temp) > 1){
    test <- as.matrix(do.call(rbind, strsplit(as.vector(temp$loci), split = "[-:]"))[,2:5])
  }else{
    test <- t(as.matrix(do.call(rbind, strsplit(as.vector(temp$loci), split = "[-:]"))[,2:5]))
  }
  
  mode(test) <- "numeric"
  for(i in 1:nrow(test)){
    ASrange <- IRanges(test[i,2]+1, test[i,3]-1)
    ## because the positions of SJ are intron sets, so we  +- 1 to match exon sets
    ASlen <- diff(c(test[i,2], test[i,3]))
    ov_within<-countOverlaps(exon_chr, ASrange, type="within")
    ov_start<-countOverlaps(exon_chr, ASrange, type="start")
    ov_end<-countOverlaps(exon_chr, ASrange, type="end")
    ov_equal<-countOverlaps(exon_chr, ASrange, type="equal")
    ov_any<-countOverlaps(exon_chr,ASrange, type="any")
    if (sum(ov_equal)>0) {
      AStype <- "exon-skipping_exactly"
      
    } else if (sum(ov_equal)==0 & sum(ov_within)>0) {
      AStype <- "exon-skipping_multipul-exon"
      
    } else if (sum(ov_equal)==0 & sum(ov_within)==0 & sum(ov_start)>0 | sum(ov_equal)==0 & sum(ov_within)==0 & sum(ov_end)>0){
      AStype <- "exon-skipping_half-exon"
      
    } else if (sum(ov_any)==0){
      AStype <- "intergenic-splicing_site"
      
    } else {
      AStype <- "Other"
    }
    SE_type[[j]][i,"ASrange"] <- paste("chr", chr[j], ":", test[i,2],"-",test[i,3],sep="")
    SE_type[[j]][i,"AStype"] <- AStype
    SE_type[[j]][i,"OverlapCount"] <- paste(sum(ov_within),"|",sum(ov_start),"|",sum(ov_end),"|",sum(ov_equal),"|",sum(ov_any),sep="")
    SE_type[[j]][i,"ASlen"] <- ASlen
    #print(paste(i," chr", chr[j], ":", test[i,2],"-",test[i,3]," ", AStype," - ",ASlen,sep=""))
  }
}
SE_type <- do.call(rbind, SE_type)


SE_type[,"PSI.R"] <- psi.co[as.character(SE_type$right),1]
SE_type[,"PSI.L"] <- psi.co[as.character(SE_type$left),1]
SE_type[,"PSI"] <- rowMeans(cbind(psi.co[as.character(SE_type$left),1], psi.co[as.character(SE_type$right),1]), na.rm = T)


#### 5. find A3SS -------------------------------------------------------------------------------------------------------------------------------
print(paste("Step5: find A3SS"))


junc.as.2.ot_uniq <- junc.as.2.ot[!duplicated(junc.as.2.ot),]
junc.as.2.ot_uniq$event <- do.call(rbind, strsplit(rownames(junc.as.2.ot_uniq), split = "\\."))[,1]

a3_raw <- junc.as.2.ot_uniq[junc.as.2.ot_uniq$start %in% names(table(junc.as.2.ot_uniq$start)[table(junc.as.2.ot_uniq$start)==2]), ]

a3_raw <- a3_raw[a3_raw$event %in% SE$left == FALSE, ]

a3_raw <- dlply(a3_raw, c("chr","start"))
#sum(unlist(lapply(a3_raw,function(x) length(unique(x$chr)) != 1)))
# if != 0, there are some problem!!!

names(a3_raw) <- sapply(a3_raw,function(x) paste(x$names,collapse="_"))

a3 <- data.frame(row.names = names(a3_raw), name = names(a3_raw),
                 distan = as.character(sapply(a3_raw,function(x) diff(as.numeric(x$end)))),
                 chr = as.character(sapply(a3_raw, function(x) unique(x$chr))),
                 start = sapply(a3_raw, function(x) mean(as.numeric(x$start))),
                 end1 = sapply(a3_raw, function(x) min(as.numeric(x$end))),
                 end2 = sapply(a3_raw, function(x) max(as.numeric(x$end))),stringsAsFactors=F)


A3SS <- a3[rownames(a3) %in% rownames(psi.co),]
#identical(row.names(psi.co)[which(row.names(psi.co) %in% row.names(A3SS))],rownames(A3SS))

psi <- psi.co[as.character(row.names(A3SS)),]

A3SS_psi <- cbind(A3SS, psi)

A3SS_psi <- A3SS_psi[order(A3SS_psi$chr,A3SS_psi$start,A3SS_psi$end1,A3SS_psi$end2),]


#### 6.validate A3SS using GTF --------------------------------------------------------
print(paste("Step6: validate A3SS using GTF"))




a3_od <- A3SS_psi


#unique(a3_od$chr)
chr=unique(a3_od$chr)
A3SS_type <- dlply(a3_od[,1:7], .variables = "chr")

for(j in 1:length(chr)){
  #print(chr[j])
  eval(parse(text = paste("exon_cn$'",chr[j],"'->exon_chr",sep="")))
  temp <- A3SS[A3SS$chr == chr[j],]
  if(nrow(temp) > 1){
    test <- as.matrix(temp[,3:6])
  }else{
    test <- as.matrix(temp[,3:6])
  }
  
  mode(test) <- "numeric"
  for(i in 1:nrow(test)){
    ASrange <- IRanges(test[i,3]+1, test[i,4]+1)
    ASlen <- diff(c(test[i,3], test[i,4]))
    ov_within<-countOverlaps(exon_chr, ASrange, type="within")
    ov_start<-countOverlaps(exon_chr, ASrange, type="start")
    ov_end<-countOverlaps(exon_chr, ASrange, type="end")
    ov_equal<-countOverlaps(exon_chr, ASrange, type="equal")
    ov_any<-countOverlaps(exon_chr,ASrange, type="any")
    
    if(sum(ov_start)>0 & sum(ov_within) == 0){
      AStype<-"alt_5/3-splicing_site"
    }else if(sum(ov_start)>0 & sum(ov_within) == 1){
      AStype<-"exon-skipping_exactly"
    }else if(sum(ov_start)>0 & sum(ov_within) > 1){
      AStype<-"exon-skipping_multiple-exon"
    }else if(sum(ov_start) == 0 & sum(ov_within) == 0 & sum(ov_any) == 1){
      AStype<-"exon-skipping_skipped_one_unannotated_exon"
    }else if(sum(ov_start) == 0 & sum(ov_within) > 0 & sum(ov_any) == sum(ov_within) + 1){
      AStype<-"exon-skipping_multiple-exon_and_first_exon_unannotated"
    }else if(sum(ov_start) == 0 & sum(ov_within) > 0 & sum(ov_any) == sum(ov_within)){
      AStype<-"exon-skipping_multiple-exon_and_first_and_last_exon_unannotated"
    }else if (sum(ov_any) == 0){
      AStype<-"intergenic-alt_5/3-splicing_site"
    }else{
      AStype<-"Others"
    }
    
    A3SS_type[[j]][i,"ASrange"] <- paste("chr", chr[j], ":", test[i,3],"-",test[i,4],sep="")
    A3SS_type[[j]][i,"AStype"] <- AStype
    A3SS_type[[j]][i,"OverlapCount"] <- paste(sum(ov_within),"|",sum(ov_start),"|",sum(ov_end),"|",sum(ov_equal),"|",sum(ov_any),sep="")
    A3SS_type[[j]][i,"ASlen"] <- ASlen
    #print(paste(i," chr", chr[j], ":", test[i,3],"-",test[i,4]," ", AStype," - ",ASlen,sep=""))
  }
}

A3SS_type <- do.call(rbind, A3SS_type)

#table(A3SS_type$AStype)



#### 7. find A5SS -------------------------------------------------------------------------------------------------------------------------------
print(paste("Step7: find A5SS "))





junc.as.2.ot_uniq <- junc.as.2.ot[!duplicated(junc.as.2.ot, fromLast = TRUE),]

junc.as.2.ot_uniq$event <- do.call(rbind, strsplit(rownames(junc.as.2.ot_uniq), split = "\\."))[,1]

a5_raw <- junc.as.2.ot_uniq[junc.as.2.ot_uniq$end %in% names(table(junc.as.2.ot_uniq$end)[table(junc.as.2.ot_uniq$end)==2]), ]

a5_raw <- a5_raw[order(a5_raw$chr, a5_raw$end, a5_raw$start),]

a5_raw <- a5_raw[a5_raw$event %in% SE$right == FALSE, ]

a5_raw <- dlply(a5_raw, c("chr","end"))
#sum(unlist(lapply(a5_raw,function(x) length(unique(x$chr)) != 1)))
# if != 0, there are some problem!!!


names(a5_raw) <- sapply(a5_raw,function(x) paste(x$names,collapse="_"))

a5 <- data.frame(name = names(a5_raw), row.names = names(a5_raw),
                 distan = as.integer(sapply(a5_raw,function(x) diff(as.numeric(x$start)))), 
                 chr = as.character(sapply(a5_raw, function(x) unique(x$chr))),
                 start1 = sapply(a5_raw, function(x) min(as.numeric(x$start))),
                 start2 = sapply(a5_raw, function(x) max(as.numeric(x$start))),
                 end = sapply(a5_raw, function(x) mean(as.numeric(x$end))), stringsAsFactors = F)


A5SS <- a5[rownames(a5) %in% rownames(psi.co),]

#identical(rownames(psi.co[row.names(A5SS),]),rownames(A5SS))
#identical(row.names(psi.co)[which(row.names(psi.co) %in% row.names(A5SS))],rownames(A5SS))

psi <- psi.co[as.character(row.names(A5SS)),]
A5SS_psi <- cbind(A5SS, psi)

A5SS_psi <- A5SS_psi[order(A5SS_psi$chr,A5SS_psi$end,A5SS_psi$start1,A5SS_psi$start2),]

A5SS <- A5SS[order(A5SS$chr,A5SS$end,A5SS$start1,A5SS$start2),]



#### 8.validate A5SS using GTF --------------------------------------------------------
print(paste("Step8: validate A5SS using GTF "))




a5_od <- A5SS_psi
chr=unique(a5_od$chr)
A5SS_type <- dlply(a5_od[,1:7], .variables = "chr")

for(j in 1:length(chr)){
  #print(chr[j])
  eval(parse(text = paste("exon_cn$'",chr[j],"'->exon_chr",sep="")))
  temp <- A5SS[A5SS$chr == chr[j],]
  if(nrow(temp) > 1){
    test <- as.matrix(temp[,3:6])
  }else{
    test <- as.matrix(temp[,3:6])
  }
  
  mode(test) <- "numeric"
  for(i in 1:nrow(test)){
    ASrange <- IRanges(test[i,2]-1, test[i,3]-1)
    ASlen <- diff(c(test[i,2], test[i,3]))
    ov_within<-countOverlaps(exon_chr, ASrange, type="within")
    ov_start<-countOverlaps(exon_chr, ASrange, type="start")
    ov_end<-countOverlaps(exon_chr, ASrange, type="end")
    ov_equal<-countOverlaps(exon_chr, ASrange, type="equal")
    ov_any<-countOverlaps(exon_chr,ASrange, type="any")
    
    if(sum(ov_end)>0 & sum(ov_within) == 0){
      AStype<-"alt_5/3-splicing_site"
    }else if(sum(ov_end)>0 & sum(ov_within) == 1){
      AStype<-"exon-skipping"
    }else if(sum(ov_end)>0 & sum(ov_within) > 1){
      AStype<-"exon-skipping_multiple-exon"
    }else if(sum(ov_end) == 0 & sum(ov_within) == 0 & sum(ov_any) == 1){
      AStype<-"exon-skipping_skipped_one_unannotated_exon"
    }else if(sum(ov_end) == 0 & sum(ov_within) > 0 & sum(ov_any) == sum(ov_within) + 1){
      AStype<-"exon-skipping_multiple-exon_and_last_exon_unannotated"
    }else if(sum(ov_end) == 0 & sum(ov_within) > 0 & sum(ov_any) == sum(ov_within)){
      AStype<-"exon-skipping_multiple-exon_and_first_and_last_exon_unannotated"
    }else if (sum(ov_any) == 0){
      AStype<-"intergenic-alt_5/3-splicing_site"
    }else{
    	AStype<-"Others"
    }
    
    A5SS_type[[j]][i,"ASrange"] <- paste("chr", chr[j], ":", test[i,3],"-",test[i,4],sep="")
    A5SS_type[[j]][i,"AStype"] <- AStype
    A5SS_type[[j]][i,"OverlapCount"] <- paste(sum(ov_within),"|",sum(ov_start),"|",sum(ov_end),"|",sum(ov_equal),"|",sum(ov_any),sep="")
    A5SS_type[[j]][i,"ASlen"] <- ASlen
    #print(paste(i," chr", chr[j], ":", test[i,3],"-",test[i,4]," ", AStype," - ",ASlen,sep=""))
  }
}


A5SS_type <- do.call(rbind, A5SS_type)




#### 9.find MXE ---------------------------------------------------------------------------------
print(paste("Step8: find MXE"))




start_pos <- names(table(junc.as.2.ot$start)[table(junc.as.2.ot$start)==2])

end_pos <- names(table(junc.as.2.ot$end)[table(junc.as.2.ot$end)==2])


## start1 < end2 < start2
## end1 < start1 < end2
## start < end1 < start1 < end2 < start2 < end


a3_for_mxe <- data.frame(a3[order(a3$chr, a3$start, a3$end1, a3$end2),])
a5_for_mxe <- data.frame(a5[order(a5$chr, a5$start1, a5$start2, a5$end),])

library(dplyr)

a5_for_mxe_rep <- bind_rows(replicate(nrow(a3_for_mxe), a5_for_mxe, simplify = FALSE))

library(tidyr)

a3_for_mxe_rep <- bind_rows(replicate(nrow(a5_for_mxe), a3_for_mxe, simplify = FALSE))

a3_for_mxe_rep <- data.frame(a3_for_mxe_rep[order(a3_for_mxe_rep$chr, a3_for_mxe_rep$start, a3_for_mxe_rep$end1, a3_for_mxe_rep$end2),])


mxe_raw <- cbind(a3_for_mxe_rep,a5_for_mxe_rep)

MXE <- mxe_raw[mxe_raw[,"start1"] < mxe_raw[,"end2"] & mxe_raw[,"start1"] > mxe_raw[,"end1"] &
                 mxe_raw[,"end2"] < mxe_raw[,"start2"] & mxe_raw[,"end2"] > mxe_raw[,"start1"], ]


## we just neeed the A3 and A5 from the same chr
MXE <- MXE[apply(MXE[,colnames(MXE) == "chr"], 1, 
                         function(x) sum(duplicated(as.character(x))) > 0),]

if(nrow(MXE) != 0){


#### 9.validate MXE using GTF ---------------------------------------------------------------------------------
print(paste("Step9: validate MXE using GTF"))



unique(MXE$chr)
chr=unique(MXE$chr)

MXE_type <- dlply(MXE, .variables = "chr")

for(j in 1:length(chr)){
  #print(chr[j])
  eval(parse(text = paste("exon_cn$'",chr[j],"'->exon_chr",sep="")))
  temp <- MXE[MXE$chr == chr[j],]
  
  test <- as.matrix(temp[,c(3:6,10:12)])
  
  mode(test) <- "numeric"
  
  for(i in 1:nrow(test)){
    ASrange <- IRanges(test[i,5], test[i,4])
    ASlen <- diff(c(test[i,5], test[i,4]))
    ov_within<-countOverlaps(exon_chr, ASrange, type="within")
    ov_start<-countOverlaps(exon_chr, ASrange, type="start")
    ov_end<-countOverlaps(exon_chr, ASrange, type="end")
    ov_equal<-countOverlaps(exon_chr, ASrange, type="equal")
    ov_any<-countOverlaps(exon_chr,ASrange, type="any")
    
    if(sum(ov_any) == 0){
      AStype<-"Mutually-exclusive-exon"
    }else {
      AStype<-"Mixed-splicing"
    }
    
    MXE_type[[j]][i,"ASrange"] <- paste("chr", chr[j], ":", test[i,5],"-",test[i,4],sep="")
    MXE_type[[j]][i,"AStype"] <- AStype
    MXE_type[[j]][i,"OverlapCount"] <- paste(sum(ov_within),"|",sum(ov_start),"|",sum(ov_end),"|",sum(ov_equal),"|",sum(ov_any),sep="")
    MXE_type[[j]][i,"ASlen"] <- ASlen
    #print(paste(i," chr", chr[j], ":", test[i,5],"-",test[i,4]," ", AStype," - ",ASlen,sep=""))
  }
}


MXE_type <- do.call(rbind, MXE_type)

}else{MXE_type <- NULL}

SE_result[[bl1]] <- SE_type
A3SS_result[[bl1]] <- A3SS_type
A5SS_result[[bl1]] <- A5SS_type
MXE_result[[bl1]] <- MXE_type

}


save(SE_result, A3SS_result, A5SS_result, MXE_result, file = paste(pfro,"DSU_result_cell_by_cell.RData", sep=""))



#### merge SE result =======================================================================

SE_result -> SE_result_1

for (i in 1:length(SE_result)) {
	SE_result[[i]]$id <- colnames(SJ_tu)[i]
}

test <- do.call(rbind, SE_result)

test[!duplicated(test[,1:11]),][,1:11] -> test2

test2$loci <- as.character(test2$loci)

test[test$loci %in% test2$loci[1], c("PSI","id")]


test3 <- data.frame(matrix(NA, nrow = nrow(test2), ncol = ncol(SJ_tu)))
colnames(test3) <- colnames(SJ_tu)

cbind(test2, test3) -> test3


test[test$loci %in% test2$loci[1], c("PSI","id")]

for(i in 1:sum(test$loci %in% test2$loci[1])){
	test[test$loci %in% test2$loci[1], c("PSI","id")] -> test4
	test3[1,test4[i,"id"]] <- test4[i,"PSI"]
}



for(j in 1:nrow(test3)){
  
  test[test$loci %in% test3$loci[j], c("PSI","id")] -> test4
  
  for(i in 1:nrow(test4)){
    
    test3[j,test4[i,"id"]] <- test4[i,"PSI"]
  }
  
}


test3 -> all_cells_SE_type_and_PSI

test3[,12:ncol(test3)] -> all_cells_SE_PSI
row.names(all_cells_SE_PSI) <- as.character(test3$loci)

rowSums(!is.na(all_cells_SE_PSI)) -> all_cells_SE_PSI_rowsum

length(all_cells_SE_PSI_rowsum)
# [1] 15845
sum(all_cells_SE_PSI_rowsum>1)
# [1] 7096



test3[,1:11] -> all_cells_SE_type


save(all_cells_SE_type_and_PSI, file = "/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/SE/RData/all_cells_SE_type_and_PSI.RData")


#### 6.Find the skipping exon and original gene ----------------------------------------------------------------------------------------------------------


## all chr

all_cells_SE_type -> SE_type

chr <- unique(SE_type$chr)

SE_gtf_chr <- list()

system.time(for(j in 1:length(chr)){
  print(paste("chr",chr[j]))
  eval(parse(text = paste("exon_cn$'",chr[j],"'->exon_chr",sep="")))
  temp <- SE_type[SE_type$chr == chr[j],]
  
  eval(parse(text = paste("exon_ca$'",chr[j],"'->exon_ca_chr",sep="")))
  
  if(nrow(temp) > 1){
    test <- as.matrix(do.call(rbind, strsplit(as.vector(temp$loci), split = "[-:]"))[,2:5])
  }else{
    test <- t(as.matrix(do.call(rbind, strsplit(as.vector(temp$loci), split = "[-:]"))[,2:5]))
  }
  colnames(test) <- paste("V", 1:4, sep = "")
  mode(test) <- "numeric"
  
  se_gtf_chr <- list()
  
  for(i in 1:nrow(test)){
    #print(i)
    ASrange <- IRanges(test[i,2]+1, test[i,3]-1)
    ## because the positions of SJ are intron sets, so we  +- 1 to match exon sets
    
    hit_within <- queryHits(findOverlaps(exon_chr, ASrange, type="within"))
    hit_start <- queryHits(findOverlaps(exon_chr, ASrange, type="start"))
    hit_end <- queryHits(findOverlaps(exon_chr, ASrange, type="end"))
    hit_equal <- queryHits(findOverlaps(exon_chr, ASrange, type="equal"))
    
    if(length(hit_within) > 0){
      within_all <- data.frame(OverType = "within", exon_ca_chr[hit_within,])
    }else{
      within_all <- data.frame(OverType = "within", t(data.frame(rep(NA,9), row.names = paste("V", 1:9, sep = ""))))
    }
    
    if(length(hit_start) >0 ){
      start_all <- data.frame(OverType = "start", exon_ca_chr[hit_start,])
    }else{
      start_all <- data.frame(OverType = "start", t(data.frame(rep(NA,9), row.names = paste("V", 1:9, sep = ""))))
    }
    
    if(length(hit_end) > 0){
      end_all <- data.frame(OverType = "end", exon_ca_chr[hit_end,])
    }else{
      end_all <- data.frame(OverType = "end", t(data.frame(rep(NA,9), row.names = paste("V", 1:9, sep = ""))))
    }
    
    if(length(hit_equal) > 0){
      equal_all <- data.frame(OverType = "equal", exon_ca_chr[hit_equal,])
    }else{
      equal_all <- data.frame(OverType = "equal", t(data.frame(rep(NA,9), row.names = paste("V", 1:9, sep = ""))))
    }
    
    se_gtf_chr[[i]] <- suppressWarnings(cbind(temp[i,], rbind(within_all,start_all,end_all,equal_all)))
    
  }
  
  SE_gtf_chr[[j]] <- suppressMessages(bind_rows(se_gtf_chr))
  gc()
})
#     user   system  elapsed 
# 4079.416   40.844 4119.096 

system.time(SE_GTF <- bind_rows(SE_gtf_chr))

dim(SE_type)
# [1] 15845   14
dim(SE_GTF)
# [1] 111400    24



save(SE_GTF, file = "/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/SE/RData/all_cells_SE_GTF.RData")





#### merge A3SS result =======================================================================


for (i in 1:length(A3SS_result)) {
	A3SS_result[[i]]$id <- colnames(SJ_tu)[i]
}

test <- do.call(rbind, A3SS_result)

test[!duplicated(test[,c(1:6,8:11)]),][,c(1:6,8:11)] -> test2
table(test2$AStype)
#
#                                          alt_5/3-splicing_site
#                                                           4995
#                                          exon-skipping_exactly
#                                                           5539
#                                    exon-skipping_multiple-exon
#                                                          11524
#exon-skipping_multiple-exon_and_first_and_last_exon_unannotated
#                                                           5249
#         exon-skipping_multiple-exon_and_first_exon_unannotated
#                                                           1082
#                     exon-skipping_skipped_one_unannotated_exon
#                                                           1937
#                               intergenic-alt_5/3-splicing_site
#                                                           1669
#                                                         Others
#                                                          4284

test2$AStype <- as.character(test2$AStype)
test2 <- test2[test2$AStype %in% c("alt_5/3-splicing_site","exon-skipping_exactly","intergenic-alt_5/3-splicing_site"),]

test2$name <- as.character(test2$name)

test[test$name %in% test2$name[1], c("psi","id")]


test3 <- data.frame(matrix(NA, nrow = nrow(test2), ncol = ncol(SJ_tu)))
colnames(test3) <- colnames(SJ_tu)

cbind(test2, test3) -> test3


test[test$loci %in% test2$loci[1], c("PSI","id")]

for(i in 1:sum(test$loci %in% test2$loci[1])){
	test[test$loci %in% test2$loci[1], c("PSI","id")] -> test4
	test3[1,test4[i,"id"]] <- test4[i,"PSI"]
}



for(j in 1:nrow(test3)){
  print(paste(j, "of", nrow(test3)))
  test[test$name %in% test3$name[j], c("psi","id")] -> test4
  
  for(i in 1:nrow(test4)){
    
    test3[j,test4[i,"id"]] <- test4[i,"psi"]
    
  }
  
}


test3 -> all_cells_A3SS_type_and_PSI

test3[,11:ncol(test3)] -> all_cells_A3SS_PSI
row.names(all_cells_A3SS_PSI) <- as.character(test3$name)

rowSums(!is.na(all_cells_A3SS_PSI)) -> all_cells_A3SS_PSI_rowsum

length(all_cells_A3SS_PSI_rowsum)
# [1] 12203
sum(all_cells_A3SS_PSI_rowsum>1)
# [1] 5856
sum(all_cells_A3SS_PSI_rowsum[all_cells_A3SS_type_and_PSI$AStype == "alt_5/3-splicing_site"]>1)
# [1] 2767


test3[,1:10] -> all_cells_A3SS_type


save(all_cells_A3SS_type_and_PSI, file = "/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/A3SS/RData/all_cells_A3SS_type_and_PSI.RData")



#### 9.Find the A3SS exon and original gene ----------------------------------------------------------------------------------------------------------


all_cells_A3SS_type -> A3SS_type


chr <- unique(A3SS_type$chr)

A3SS_gtf_chr <- list()

system.time(for(j in 1:length(chr)){
  print(paste("chr",chr[j]))
  eval(parse(text = paste("exon_cn$'",chr[j],"'->exon_chr",sep="")))
  temp <- A3SS_type[A3SS_type$chr == chr[j],]
  
  eval(parse(text = paste("exon_ca$'",chr[j],"'->exon_ca_chr",sep="")))
  
  test <- as.matrix(temp[,3:6])
  
  mode(test) <- "numeric"
  
  a3_gtf_chr <- list()
  
  for(i in 1:nrow(test)){
    #print(i)
    ASrange <- IRanges(test[i,3]+1, test[i,4]+1)
    ## because the positions of SJ are intron sets, so we  +- 1 to match exon sets
    
    hit_within <- queryHits(findOverlaps(exon_chr, ASrange, type="within"))
    hit_start <- queryHits(findOverlaps(exon_chr, ASrange, type="start"))
    hit_end <- queryHits(findOverlaps(exon_chr, ASrange, type="end"))
    hit_equal <- queryHits(findOverlaps(exon_chr, ASrange, type="equal"))
    
    if(length(hit_within) > 0){
      within_all <- data.frame(OverType = "within", exon_ca_chr[hit_within,])
    }else{
      within_all <- data.frame(OverType = "within", t(data.frame(rep(NA,9), row.names = paste("V", 1:9, sep = ""))))
    }
    
    if(length(hit_start) >0 ){
      start_all <- data.frame(OverType = "start", exon_ca_chr[hit_start,])
    }else{
      start_all <- data.frame(OverType = "start", t(data.frame(rep(NA,9), row.names = paste("V", 1:9, sep = ""))))
    }
    
    if(length(hit_end) > 0){
      end_all <- data.frame(OverType = "end", exon_ca_chr[hit_end,])
    }else{
      end_all <- data.frame(OverType = "end", t(data.frame(rep(NA,9), row.names = paste("V", 1:9, sep = ""))))
    }
    
    if(length(hit_equal) > 0){
      equal_all <- data.frame(OverType = "equal", exon_ca_chr[hit_equal,])
    }else{
      equal_all <- data.frame(OverType = "equal", t(data.frame(rep(NA,9), row.names = paste("V", 1:9, sep = ""))))
    }
    
    a3_gtf_chr[[i]] <- suppressWarnings(cbind(temp[i,], rbind(within_all,start_all,end_all,equal_all)))
    
  }
  
  A3SS_gtf_chr[[j]] <- suppressMessages(bind_rows(a3_gtf_chr))
  gc()
})
#     user   system  elapsed 
# 4079.416   40.844 4119.096 

system.time(A3SS_GTF <- bind_rows(A3SS_gtf_chr))


save(A3SS_GTF, file = "/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/A3SS/RData/all_cells_A3SS_GTF.RData")





#### merge A5SS result =======================================================================


for (i in 1:length(A5SS_result)) {
	A5SS_result[[i]]$id <- colnames(SJ_tu)[i]
}

test <- do.call(rbind, A5SS_result)

test[!duplicated(test[,c(1:6,8:11)]),][,c(1:6,8:11)] -> test2
table(test2$AStype)
#
#                                          alt_5/3-splicing_site
#                                                           4785
#                                          exon-skipping_exactly
#                                                           5481
#                                    exon-skipping_multiple-exon
#                                                          11868
#exon-skipping_multiple-exon_and_first_and_last_exon_unannotated
#                                                           5300
#         exon-skipping_multiple-exon_and_first_exon_unannotated
#                                                           1054
#                     exon-skipping_skipped_one_unannotated_exon
#                                                           1783
#                               intergenic-alt_5/3-splicing_site
#                                                           1707
#                                                         Others
#                                                           4354

test2$AStype <- as.character(test2$AStype)
test2 <- test2[test2$AStype %in% c("alt_5/3-splicing_site","exon-skipping","intergenic-alt_5/3-splicing_site"),]

test2$name <- as.character(test2$name)

test[test$name %in% test2$name[1], c("psi","id")]


test3 <- data.frame(matrix(NA, nrow = nrow(test2), ncol = ncol(SJ_tu)))
colnames(test3) <- colnames(SJ_tu)

cbind(test2, test3) -> test3




for(j in 1:nrow(test3)){
  print(paste(j, "of", nrow(test3)))
  test[test$name %in% test3$name[j], c("psi","id")] -> test4
  
  for(i in 1:nrow(test4)){
    
    test3[j,test4[i,"id"]] <- test4[i,"psi"]
    
  }
  
}


test3 -> all_cells_A5SS_type_and_PSI

test3[,11:ncol(test3)] -> all_cells_A5SS_PSI
row.names(all_cells_A5SS_PSI) <- as.character(test3$name)

rowSums(!is.na(all_cells_A5SS_PSI)) -> all_cells_A5SS_PSI_rowsum

length(all_cells_A5SS_PSI_rowsum)
# [1] 11973
sum(all_cells_A5SS_PSI_rowsum>1)
# [1] 5728
sum(all_cells_A5SS_PSI_rowsum[all_cells_A5SS_type_and_PSI$AStype == "alt_5/3-splicing_site"]>1)
# [1] 2667


test3[,1:10] -> all_cells_A5SS_type


save(all_cells_A5SS_type_and_PSI, file = "/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/A5SS/RData/all_cells_A5SS_type_and_PSI.RData")





#### 9.Find the A3SS exon and original gene ----------------------------------------------------------------------------------------------------------




all_cells_A5SS_type -> A5SS_type


chr <- unique(A5SS_type$chr)

A5SS_gtf_chr <- list()

system.time(for(j in 1:length(chr)){
  print(paste("chr",chr[j]))
  eval(parse(text = paste("exon_cn$'",chr[j],"'->exon_chr",sep="")))
  temp <- A5SS_type[A5SS_type$chr == chr[j],]
  
  eval(parse(text = paste("exon_ca$'",chr[j],"'->exon_ca_chr",sep="")))
  
  test <- as.matrix(temp[,3:6])
  
  mode(test) <- "numeric"
  
  a5_gtf_chr <- list()
  
  for(i in 1:nrow(test)){
    #print(i)
    ASrange <- IRanges(test[i,2]-1, test[i,3]-1)
    ## because the positions of SJ are intron sets, so we  +- 1 to match exon sets
    
    hit_within <- queryHits(findOverlaps(exon_chr, ASrange, type="within"))
    hit_start <- queryHits(findOverlaps(exon_chr, ASrange, type="start"))
    hit_end <- queryHits(findOverlaps(exon_chr, ASrange, type="end"))
    hit_equal <- queryHits(findOverlaps(exon_chr, ASrange, type="equal"))
    
    if(length(hit_within) > 0){
      within_all <- data.frame(OverType = "within", exon_ca_chr[hit_within,])
    }else{
      within_all <- data.frame(OverType = "within", t(data.frame(rep(NA,9), row.names = paste("V", 1:9, sep = ""))))
    }
    
    if(length(hit_start) >0 ){
      start_all <- data.frame(OverType = "start", exon_ca_chr[hit_start,])
    }else{
      start_all <- data.frame(OverType = "start", t(data.frame(rep(NA,9), row.names = paste("V", 1:9, sep = ""))))
    }
    
    if(length(hit_end) > 0){
      end_all <- data.frame(OverType = "end", exon_ca_chr[hit_end,])
    }else{
      end_all <- data.frame(OverType = "end", t(data.frame(rep(NA,9), row.names = paste("V", 1:9, sep = ""))))
    }
    
    if(length(hit_equal) > 0){
      equal_all <- data.frame(OverType = "equal", exon_ca_chr[hit_equal,])
    }else{
      equal_all <- data.frame(OverType = "equal", t(data.frame(rep(NA,9), row.names = paste("V", 1:9, sep = ""))))
    }
    
    a5_gtf_chr[[i]] <- suppressWarnings(cbind(temp[i,], rbind(within_all,start_all,end_all,equal_all)))
    
  }
  
  A5SS_gtf_chr[[j]] <- suppressMessages(bind_rows(a5_gtf_chr))
  gc()
})
#     user   system  elapsed 
# 4079.416   40.844 4119.096 

system.time(A5SS_GTF <- bind_rows(A5SS_gtf_chr))


save(A5SS_GTF, file = "/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/A5SS/RData/all_cells_A5SS_GTF.RData")






#### merge MXE result =======================================================================


for (i in 1:length(MXE_result)) {
	MXE_result[[i]]$id <- colnames(SJ_tu)[i]
}

unlist(lapply(MXE_result,is.data.frame)) -> cells_have_MXE
sum(cells_have_MXE)
#[1] 3113
length(cells_have_MXE)
#[1] 3574

test <- do.call(rbind, MXE_result[cells_have_MXE])

test[!duplicated(test[,1:16]),][,1:16] -> test2
table(test2$AStype)
#          Mixed-splicing Mutually-exclusive-exon
#                    1172                     376

test2 <- test2[test2$AStype == "Mutually-exclusive-exon",]

test2$name <- as.character(test2$name)


test3 <- data.frame(matrix(0, nrow = nrow(test2), ncol = ncol(SJ_tu)))
colnames(test3) <- colnames(SJ_tu)

cbind(test2, test3) -> test3




for(j in 1:nrow(test3)){
  print(paste(j, "of", nrow(test3)))
  unique(test[test$name %in% test3$name[j], "id"]) -> test4
  
  test3[j,test4] <- 1
  
}

length(as.numeric(rowSums(test3[,17:ncol(test3)])))
# [1] 376
sum(as.numeric(rowSums(test3[,17:ncol(test3)]))>1)
# [1] 154
sum(as.numeric(rowSums(test3[,17:ncol(test3)]))>=10)
# [1] 38

test3 -> all_cells_MXE_type_and_existing


save(all_cells_MXE_type_and_existing, file = "/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/MXE/RData/all_cells_MXE_type_and_existing.RData")





























