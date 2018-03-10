
## 根据 PSI.R 和 PSI.L 是否都存在，判断某 SE events 是否至少在一个细胞中完整存在。
## 结果表明任一 SE 都至少在一个细胞中完整存在。

load("./SE/RData/SE_DSU_10reads_2cells_result_20180228.RData")
ls()
dim(SE_type)
SE_type[1:5,1:10]
SE_type[1:5,21:30]
SE_type[1,]
colnames(SE_type)[1:30]
SE_type[,12:ncol(SE_type)] -> psi
row.names(psi) <- SE_type[,1]
psi[1:5,1:5]
grep(pattern="psi",colnames(psi),invert=T)
dim(psi[,grep(pattern="psi",colnames(psi),invert=T)])
# [1] 5714 3574

rowSums(!is.na(psi[,grep(pattern="psi",colnames(psi),invert=T)])) -> psi_rowsum_1
summary(psi_rowsum_1)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    2.0    26.0    62.0   150.3   146.0  3554.0


rowSums(!is.na((psi[,grep(pattern="L_psi",colnames(psi))] + psi[,grep(pattern="R_psi",colnames(psi))])/2)) -> psi_rowsum_2
summary(psi_rowsum_2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#      2      20      52     134     130    3540


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
pfsj <- file.path("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ(all_more_than_10_more_than_1_cells)_merged_touse.RData")

## path of figure/file output
pffo <- file.path("/mnt/data5/BGI/UCB/tangchao/DSU/")

## path of RData output
pfro <- file.path("/mnt/data5/BGI/UCB/tangchao/DSU/RData/")

## path of gtf
pfgtf <- file.path("/mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.87.gtf")


load(pfsj)

SJ_tu -> count_use

count_use[is.na(count_use)] <- 0


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
#     2
# 38134


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


junc.as2.df=do.call(rbind,junc.as.2)
dim(junc.as2.df)

junc.as2.df[order(junc.as2.df$chr,junc.as2.df$start,junc.as2.df$end),][1:20,]

junc.as.2.ot <- junc.as2.df[order(junc.as2.df$chr,junc.as2.df$start,junc.as2.df$end),]
junc.as.2.ot[1:20,]
dim(junc.as.2.ot)


temp <- dlply(data.frame(chr = unique(junc.as.2.ot$chr)),"chr")


for(i in 1:length(unique(junc.as.2.ot$chr))){
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
    temp[[i]] <- se
  }
}


SE <- do.call(rbind, temp)
dim(SE)
# 5834 3
rownames(SE) <- do.call(rbind, strsplit(rownames(SE), split = "\\."))[,2]


loci <- do.call(rbind, strsplit(as.vector(SE$loci), split = "[-:]"))

SE$chr <- loci[,1]

SE$dist_left <- as.numeric(loci[,5]) - as.numeric(loci[,3])

SE$dist_right <- as.numeric(loci[,4]) - as.numeric(loci[,2])

SE$dist_insert <- as.numeric(loci[,4]) - as.numeric(loci[,3])

SE <- SE[SE$dist_insert > 0,]
dim(SE)


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
#                       4052                         123 
#exon-skipping_multipul-exon    intergenic-splicing_site 
#                        462                        973 
#                      Other
#                        104

identical(SE$loci, SE_type$loci)


## merge PSI and SE information
for(i in colnames(psi.co)){
  SE_type[,paste(i, ".L_psi", sep = "")] <- psi.co[as.character(SE_type$right),i]
  SE_type[,paste(i, ".R_psi", sep = "")] <- psi.co[as.character(SE_type$left),i]
  SE_type[,paste(i, ".Y_psi", sep = "")] <- rowMeans(cbind(psi.co[as.character(SE_type$left),i], psi.co[as.character(SE_type$right),i]), na.rm = T)
  SE_type[,paste(i, ".N_psi", sep = "")] <- rowMeans(cbind(psi.co[as.character(SE_type$left),i], psi.co[as.character(SE_type$right),i]))
}# as.vector() is necessary


SE_type[,12:ncol(SE_type)] -> psi2
row.names(psi2) <- SE_type[,1]

rowSums(!is.na(psi2[,grep(pattern="Y_psi",colnames(psi2))])) -> psi_rowsum_21
rowSums(!is.na(psi2[,grep(pattern="N_psi",colnames(psi2))])) -> psi_rowsum_22
summary(psi_rowsum_21)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    2.0    26.0    62.0   150.3   146.0  3554.0
summary(psi_rowsum_22)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#      2      20      52     134     130    3540





load("/mnt/data5/BGI/UCB/tangchao/DSU/SE/RData/SE_DSU_10reads_result_20180228.RData")
ls()
dim(SE_type)
SE_type[1:5,1:10]
SE_type[1:5,21:30]
SE_type[1,]
colnames(SE_type)[1:30]
SE_type[,9:ncol(SE_type)] -> psi
row.names(psi) <- SE_type[,1]
psi[1:5,1:5]
grep(pattern="psi",colnames(psi),invert=T)
dim(psi[,grep(pattern="psi",colnames(psi),invert=T)])
# [1] 6848 3574

rowSums(!is.na(psi[,grep(pattern="psi",colnames(psi),invert=T)])) -> psi_rowsum_1
summary(psi_rowsum_1)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   1.00   10.00   31.00   90.39   86.00 3554.00


rowSums(!is.na((psi[,grep(pattern="L_psi",colnames(psi))] + psi[,grep(pattern="R_psi",colnames(psi))])/2)) -> psi_rowsum_2
summary(psi_rowsum_2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   1.00    7.00   25.00   76.87   72.00 3529.00

head(cbind(psi_rowsum_1,psi_rowsum_2))
#                                  psi_rowsum_1 psi_rowsum_2
#1:954083-954463-954524-955922               59           58
#1:1090429-1090875-1090984-1091471            5            5
#1:1374446-1374704-1374791-1375153          488          452
#1:1485172-1485781-1522900-1523510            7            2
#1:1629312-1629384-1629567-1629638           23           17
#1:1636482-1636681-1703625-1703823           36            1























