## calculate PSI of AS3
require(parallel)
library(dplyr)
library(plyr)
library(pryr)
require(parallel)
library(ggplot2)
library(cowplot)
library(GenomicRanges)
library(reshape2)
library(vioplot)

depar <- par()

## path and name of SJ merged RData
pfsj <- file.path("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ(all_more_than_10)_merged_touse.RData")

## path of figure/file output
pffo <- file.path("/mnt/data5/BGI/UCB/tangchao/DSU/")

## path of RData output
pfro <- file.path("/mnt/data5/BGI/UCB/tangchao/DSU/RData/")

## path of gtf
pfgtf <- file.path("/mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.87.gtf")


#### 1.load data --------------------------------------------------------------------------------------------------------------------------------

load(pfsj)

SJ_tu -> count_use

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

f_dowle3(count_use)




#### 2. Identify alternative splicing events ----------------------------------------------------------------------------------------------------


junc=count_use
#junc.names=do.call(rbind,(strsplit(sub(x=rownames(junc),
#                                       pattern="^([[:alnum:]]+):([[:digit:]]+)-([[:digit:]]+)$",
#                                       replace="\\1;\\2;\\3"),split=";")))
junc.names=do.call(rbind,strsplit(rownames(junc),split="[:-]"))
#alnum-- Letters and Numbers；digit -- Numbers
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
  print(chr)
  return(c(same.start,same.end[!sapply(same.end,is.null)]))
},mc.cores=10))


juncs <- junc.as

junc.as <- junc.as[sapply(junc.as,nrow)>1]
names(junc.as) <- sapply(junc.as,function(x) paste(x$names,collapse="_"))

## junction just has two alternative starts/ends:
junc.as.2=junc.as[sapply(junc.as,nrow)==2]

table(sapply(junc.as.2,nrow))

junc.as.3=junc.as[sapply(junc.as,nrow)==3]

table(sapply(junc.as.3,nrow))


#### 3.calculate PSI of all junction 2 ----------------------------------------------------------------------------------------------------------



psi_sj2=mclapply(junc.as.2,function(sjs){
  sjs.names=sjs$names[order(as.numeric(sjs$end)-as.numeric(sjs$start))]# find the splicing_in isoform to calculate the psi
  tab=as.matrix(t(junc[sjs.names,]))
  distc=max(c(diff(as.numeric(sjs$start)),diff(as.numeric(sjs$end))))
  rs<-apply(tab, 1, function(a) sum(a))
  tsrs<-tab/rs
  return(tsrs[,1])
},mc.cores=8)

psi.co<-t(do.call(cbind,psi_sj2))

psi.co -> psi.co2

#### 3.calculate PSI of all junction 3 ----------------------------------------------------------------------------------------------------------



psi_sj3=mclapply(junc.as.3,function(sjs){
  sjs.names=sjs$names[order(as.numeric(sjs$end)-as.numeric(sjs$start))]# find the splicing_in isoform to calculate the psi
  tab=as.matrix(t(junc[sjs.names,]))
  distc=max(c(diff(as.numeric(sjs$start)),diff(as.numeric(sjs$end))))
  rs<-apply(tab, 1, function(a) sum(a))
  tsrs<-tab/rs
  return(tsrs)
},mc.cores=8)

psi_sj3 -> psi_sj3_s
colSums(na.omit(psi_sj3[[12]])>0)

colSums(psi_sj3[[12]]>0, na.rm = T)

sj_sum <- mclapply(psi_sj3, function(x) {colSums(x > 0, na.rm = T)}, mc.cores = 6)

psi_sj3[[12]][,colSums(psi_sj3[[12]]>0, na.rm = T)>1]


## We just remain the SJ that appeare in more than 1 cells of this 3 SJs 
psi_sj3_sub <- mclapply(psi_sj3, function(x){x[, colSums(x > 0, na.rm = T) > 1]}, mc.cores = 6)

psi_sj3_sub <- psi_sj3_sub[as.vector(unlist(lapply(psi_sj3_sub, is.matrix)))]
table(as.numeric(unlist(lapply(psi_sj3_sub, ncol))))
#   0    2    3
#1306 9668 4941

psi_sj3_sub_tu <- psi_sj3_sub[as.vector(unlist(lapply(psi_sj3_sub, function(x) ncol(x) == 2)))]

### Except PSI, as.sj should remove the SJs that just appeare in 1 cells too.
psi_sj3_sub <- mclapply(psi_sj3, function(x){colSums(x > 0, na.rm = T) > 1}, mc.cores = 6)

junc.as.3 -> junc.as.3_sub
for(i in 1:length(junc.as.3_sub)){
	junc.as.3_sub[[i]] <- junc.as.3_sub[[i]][junc.as.3_sub[[i]]$names %in% colnames(psi_sj3[[i]])[colSums(psi_sj3[[i]]>0, na.rm = T) > 1],]
}

table(as.vector(unlist(lapply(junc.as.3_sub, nrow))))
#    0    1    2    3
# 1306 6147 9668 4941

junc.as.3_sub_tu <- junc.as.3_sub[as.vector(unlist(lapply(junc.as.3_sub, nrow))) ==2]

names(junc.as.3_sub_tu) <- sapply(junc.as.3_sub_tu,function(x) paste(x$names,collapse="_"))


psi_sj32=mclapply(psi_sj3_sub_tu,function(sjs){
  sjs.names=colnames(sjs)# find the splicing_in isoform to calculate the psi
  tab=as.matrix(t(junc[sjs.names,]))
  #distc=max(c(diff(as.numeric(sjs$start)),diff(as.numeric(sjs$end))))
  rs<-apply(tab, 1, function(a) sum(a))
  tsrs<-tab/rs
  return(tsrs[,1])
},mc.cores=8)

names(psi_sj32) <- names(junc.as.3_sub_tu)


psi_sj3.co<-t(do.call(cbind,psi_sj32))

dim(psi_sj3.co)
#[1] 9668 3574
dim(psi.co)
#[1] 72769  3574

psi.co <- rbind(psi.co,psi_sj3.co)




#### 4.find SE ---------------------------------------------------------------------------------------------------------------------------------




junc.as2.df=do.call(rbind,junc.as.2)
dim(junc.as2.df)
# [1] 145538      4

#junc.as2.df <- rbind(junc.as2.df,do.call(rbind,junc.as.3_sub_tu))
#dim(junc.as2.df)
# [1] 164874      4



junc.as2.df[order(junc.as2.df$chr,junc.as2.df$start,junc.as2.df$end),][1:20,]

junc.as.2.ot <- junc.as2.df[order(junc.as2.df$chr,junc.as2.df$start,junc.as2.df$end),]
junc.as.2.ot[1:20,]
dim(junc.as.2.ot)


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

dim(SE)
# [1] 6848    7



#### 5.validate SE using GTF file ---------------------------------------------------




gtf <- read.table(pfgtf,head=F,sep="\t")

exon <- subset(gtf,gtf$V3=="exon")

exon <- exon[!duplicated(exon[,c(1,4,5)]),]

exon_ca <- dlply(exon, "V1")

exon_c <- exon[,c(1,4,5)]

mode(exon_c)

exon_cn <- split(IRanges(exon_c[,2], exon_c[,3]), exon_c[,1])

unique(SE$chr)
chr <- unique(SE$chr)

SE_type <- dlply(SE, .variables = "chr")



for(j in 1:length(chr)){
  print(chr[j])
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
    print(paste(i," chr", chr[j], ":", test[i,2],"-",test[i,3]," ", AStype," - ",ASlen,sep=""))
  }
}
SE_type <- do.call(rbind, SE_type)
head(SE_type)

table(SE_type$AStype)
#      exon-skipping_exactly     exon-skipping_half-exon
#                       4170                         174
#exon-skipping_multipul-exon    intergenic-splicing_site
#                        734                        1592
#                      Other
#                        178
identical(SE$loci, SE_type$loci)
# TRUE

SE_type -> SE_type_1

#### 5.find SE of sj3 subset SJs ---------------------------------------------------------------------------------------------------------------------------------


junc.as2.df <- do.call(rbind,junc.as.3_sub_tu)
dim(junc.as2.df)
# [1] 19336     4

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

dim(SE)
# [1] 33  7


#### 7.validate SE using GTF file ---------------------------------------------------

chr <- unique(SE$chr)

SE_type <- dlply(SE, .variables = "chr")



for(j in 1:length(chr)){
  print(chr[j])
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
    print(paste(i," chr", chr[j], ":", test[i,2],"-",test[i,3]," ", AStype," - ",ASlen,sep=""))
  }
}
SE_type <- do.call(rbind, SE_type)
head(SE_type)

table(SE_type$AStype)
#      exon-skipping_exactly     exon-skipping_half-exon
#                         23                           3
#exon-skipping_multipul-exon    intergenic-splicing_site
#                          2                           5

identical(SE$loci, SE_type$loci)

SE_type -> SE_type_2


SE_type <- rbind(SE_type_1, SE_type_2)

table(SE_type$AStype)
#      exon-skipping_exactly     exon-skipping_half-exon
#                       4193                         177
#exon-skipping_multipul-exon    intergenic-splicing_site
#                        736                        1597
#                      Other
#                        178



#### 6.Find the skipping exon and original gene ----------------------------------------------------------------------------------------------------------



## all chr


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
})
#     user   system  elapsed 
# 4079.416   40.844 4119.096 

system.time(SE_GTF <- bind_rows(SE_gtf_chr))

dim(SE)
# [1] 6848   10
dim(SE_type)
# [1] 6881   11
dim(SE_GTF)
# [1] 35706    21



## merge PSI and SE information
for(i in colnames(psi.co)){
  SE_type[,paste(i, ".L_psi", sep = "")] <- psi.co[as.character(SE_type$right),i]
  SE_type[,paste(i, ".R_psi", sep = "")] <- psi.co[as.character(SE_type$left),i]
  SE_type[,i] <- rowMeans(cbind(psi.co[as.character(SE_type$left),i], psi.co[as.character(SE_type$right),i]), na.rm = T)
}# as.vector() is necessary


save(SE_type, SE_GTF, file = paste(pffo, "SE/RData/SE_DSU_10reads_result_20180303.RData", sep = ""))




#### 7. find A3SS -------------------------------------------------------------------------------------------------------------------------------


junc.as2.df=do.call(rbind,junc.as.2)
dim(junc.as2.df)
junc.as.2.ot <- junc.as2.df[order(junc.as2.df$chr,junc.as2.df$start,junc.as2.df$end),]


junc.as.2.ot_uniq <- junc.as.2.ot[!duplicated(junc.as.2.ot),]
## 去除重复是因为要将所有的 AS SJ 都用起来，主要目的是搜集那些在 SE 的时候因为是出现 3 次但其实并不是 SE 的那些 SJ, 
## 它们之前出现三次，我们认为有可能是 SE ， 但是之后因为 AS insert 等原因判断出并不是 SE，但也出现了三次，
## 其中还有一次是重复的，所以把重复去掉。

sum(table(junc.as.2.ot_uniq$start)==2)
# [1] 5056
## 把真正的 SE 去掉：
head(junc.as.2.ot_uniq)
junc.as.2.ot_uniq$event <- do.call(rbind, strsplit(rownames(junc.as.2.ot_uniq), split = "\\."))[,1]

a3_raw <- junc.as.2.ot_uniq[junc.as.2.ot_uniq$start %in% names(table(junc.as.2.ot_uniq$start)[table(junc.as.2.ot_uniq$start)==2]), ]

sum(a3_raw$event %in% SE_type$left)
sum(a3_raw$event %in% SE_type$left == FALSE)

a3_raw <- a3_raw[a3_raw$event %in% SE_type$left == FALSE, ]
## 把真正的 SE 去掉


a3_raw <- dlply(a3_raw, c("chr","start"))
sum(unlist(lapply(a3_raw,function(x) length(unique(x$chr)) != 1)))
# if != 0, there are some problem!!!

names(a3_raw) <- sapply(a3_raw,function(x) paste(x$names,collapse="_"))

names(a3_raw)
as.vector(sapply(a3_raw,function(x) diff(as.numeric(x$end))))

cbind(names(a3_raw),as.vector(sapply(a3_raw,function(x) diff(x$end))))

sapply(a3_raw, function(x) mean(as.numeric(x$start)))
sapply(a3_raw, function(x) min(as.numeric(x$end)))
sapply(a3_raw, function(x) max(as.numeric(x$end)))
sapply(a3_raw, function(x) unique(x$chr))



a3 <- data.frame(row.names = names(a3_raw), name = names(a3_raw),
                 distan = as.character(sapply(a3_raw,function(x) diff(as.numeric(x$end)))),
                 chr = as.character(sapply(a3_raw, function(x) unique(x$chr))),
                 start = sapply(a3_raw, function(x) mean(as.numeric(x$start))),
                 end1 = sapply(a3_raw, function(x) min(as.numeric(x$end))),
                 end2 = sapply(a3_raw, function(x) max(as.numeric(x$end))),stringsAsFactors=F)



dim(a3)
head(a3)

A3SS <- a3[rownames(a3) %in% rownames(psi.co),]


identical(row.names(psi.co)[which(row.names(psi.co) %in% row.names(A3SS))],rownames(A3SS))

psi <- psi.co[as.character(row.names(A3SS)),]

A3SS_psi <- cbind(A3SS, psi)

A3SS_psi <- A3SS_psi[order(A3SS_psi$chr,A3SS_psi$start,A3SS_psi$end1,A3SS_psi$end2),]

#head(a3[rownames(a3) %in% rownames(psi.co) == FALSE,])


#### 8.validate A3SS using GTF --------------------------------------------------------




a3_od <- A3SS_psi


unique(a3_od$chr)
chr=unique(a3_od$chr)
A3SS_type <- dlply(a3_od[,1:7], .variables = "chr")

for(j in 1:length(chr)){
  print(chr[j])
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
    print(paste(i," chr", chr[j], ":", test[i,3],"-",test[i,4]," ", AStype," - ",ASlen,sep=""))
  }
}

A3SS_type <- do.call(rbind, A3SS_type)

head(A3SS_type)
dim(A3SS_type)

table(A3SS_type$AStype)
#                                          alt_5/3-splicing_site
#                                                           4112
#                                          exon-skipping_exactly
#                                                           4879
#                                    exon-skipping_multiple-exon
#                                                           7302
#exon-skipping_multiple-exon_and_first_and_last_exon_unannotated
#                                                           2107
#         exon-skipping_multiple-exon_and_first_exon_unannotated
#                                                            949
#                     exon-skipping_skipped_one_unannotated_exon
#                                                           3164
#                               intergenic-alt_5/3-splicing_site
#                                                           2209
#                                                         Others
#                                                           5130




#### 9.Find the A3SS exon and original gene ----------------------------------------------------------------------------------------------------------





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
})
#     user   system  elapsed 
# 4079.416   40.844 4119.096 

system.time(A3SS_GTF <- bind_rows(A3SS_gtf_chr))


dim(A3SS)
#[1] 29612    6
dim(A3SS_psi)
#[1] 29612 3580
dim(A3SS_type)
#[1] 29612   11
dim(A3SS_GTF)
#[1] 312693    21


save(A3SS, A3SS_psi, A3SS_type, A3SS_GTF, file = paste(pffo, "A3SS/RData/A3SS_DSU_10reads_result_20180228.RData", sep = ""))






















