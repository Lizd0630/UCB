#### date: Fri Feb 23 2018

#### server: The new server

#### **** DSU10 **** ####

#### 0. basic settings ---------------------------------------------------------------------------------------------------------------------------
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
pfsj <- file.path("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_touse.RData")

## path of figure/file output
pffo <- file.path("/mnt/data5/BGI/UCB/tangchao/DSU/")

## path of RData output
pfro <- file.path("/mnt/data5/BGI/UCB/tangchao/DSU/RData/")

## path of gtf
pfgtf <- file.path("/mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.87.gtf")


#### 1.load data --------------------------------------------------------------------------------------------------------------------------------

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


splice.count<-table(sapply(junc.as,nrow))


splice.count.sum<-sum(as.numeric(names(splice.count))*splice.count)


splice.count.percent<-round((as.numeric(names(splice.count))*splice.count)*100/splice.count.sum,2)

#todo pie chart with the percentage and numbers
pie.name <- paste("Spl ",names(splice.count), " ", splice.count," ", splice.count.percent, sep="")

pdf(paste(pffo,"Pie_chart_of_splice_count_and_percentage.pdf",sep = ""))
pie(splice.count,labels = pie.name,col=rainbow(length(pie.name)),main="Splice count and percentage")
dev.off()

juncs <- junc.as

junc.as <- junc.as[sapply(junc.as,nrow)>1]
names(junc.as) <- sapply(junc.as,function(x) paste(x$names,collapse="_"))

## junction just has two alternative starts/ends:
junc.as.2=junc.as[sapply(junc.as,nrow)==2]

table(sapply(junc.as.2,nrow))



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
save(psi.co, file = paste(pfro,"DSU10_psi.co.RData",sep = ""))

pdf(paste(pffo, "PSI_distribution.pdf", sep = ""))
plot(density(as.numeric(melt(psi.co)$value), na.rm = T), main = "Distribution of all PSI", xlab = "PSI")

plot(density(psi.co[,1], na.rm = T), main = "PSI distribution of every cell", xlab = "PSI", ylim = c(0,3.5))
for (i in 1:ncol(psi.co)) {
  lines(density(psi.co[,i], na.rm = T))
}
dev.off()



#### 4.find SE ---------------------------------------------------------------------------------------------------------------------------------



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
# 3877 4
rownames(SE) <- do.call(rbind, strsplit(rownames(SE), split = "\\."))[,2]


loci <- do.call(rbind, strsplit(as.vector(SE$loci), split = "[-:]"))

SE$chr <- loci[,1]

SE$dist_left <- as.numeric(loci[,5]) - as.numeric(loci[,3])

SE$dist_right <- as.numeric(loci[,4]) - as.numeric(loci[,2])

SE$dist_insert <- as.numeric(loci[,4]) - as.numeric(loci[,3])


## merge PSI and SE information
#for(i in colnames(psi.co)){
#  SE[,paste(i, ".L_psi", sep = "")] <- psi.co[as.character(SE$right),i]
#  SE[,paste(i, ".R_psi", sep = "")] <- psi.co[as.character(SE$left),i]
#  SE[,i] <- rowMeans(cbind(psi.co[as.character(SE$left),i], psi.co[as.character(SE$right),i]), na.rm = T)
#}# as.vector() is necessary

SE[,"PSI.R"] <- psi.co[as.character(SE$right),1]
SE[,"PSI.L"] <- psi.co[as.character(SE$left),1]
SE[,"PSI"] <- rowMeans(cbind(psi.co[as.character(SE$left),1], psi.co[as.character(SE$right),1]), na.rm = T)


## because the insert must be positive in SE, so:
## the negative insert is caused by some mixing AS
SE <- SE[SE$dist_insert > 0,]
dim(SE)
# 3789   10

pdf(paste(pffo, "SE/SE_PSI_distribution.pdf", sep = ""))
## density plot
plot(density(as.numeric(SE$PSI), na.rm = T), main = "Distribution of all PSI", xlab = "PSI")
## vioplot
vioplot(as.numeric(SE$PSI)[!is.na(as.numeric(SE$PSI))], drawRect = F)

## bean plot
beanplot(as.numeric(SE$PSI)[!is.na(as.numeric(SE$PSI))], what = c(0,1,0,0))
dev.off()




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
#                       3011                          38 
#exon-skipping_multipul-exon    intergenic-splicing_site 
#                        424                         278 
#                      Other
#                         38

identical(SE$loci, SE_type$loci)
# TRUE



pdf(paste(pffo,"SE/PSI_distribution_of_SE_type.pdf", sep = ""))
par(mar = c(2,12,2,12))
pie(table(SE_type$AStype), labels = paste(names(table(SE_type$AStype)), table(SE_type$AStype), round(table(SE_type$AStype)*100/nrow(SE_type),2), sep=" "))
par(depar)

ggplot(data = SE_type, aes(AStype, PSI))+
  geom_violin()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


#### 6.Find the skipping exon and original gene ----------------------------------------------------------------------------------------------------------


countOverlaps(exon_chr, ASrange, type="within")
queryHits(findOverlaps(exon_chr, ASrange, type="within"))


#temp <- SE_type[SE_type$chr == chr[24],]
#exon_chr <- exon_cn$Y

if(nrow(temp) > 1){
  test <- as.matrix(do.call(rbind, strsplit(as.vector(temp$loci), split = "[-:]"))[,2:5])
}else{
  test <- t(as.matrix(do.call(rbind, strsplit(as.vector(temp$loci), split = "[-:]"))[,2:5]))
}

colnames(test) <- paste("V", 1:4, sep = "")
mode(test) <- "numeric"

se_Add_gtf_chr <- list()
system.time(for(i in 1:nrow(test)){
  print(i)
  ASrange <- IRanges(test[i,2]+1, test[i,3]-1)
  ## because the positions of SJ are intron sets, so we  +- 1 to match exon sets
  
  hit_within <- queryHits(findOverlaps(exon_chr, ASrange, type="within"))
  hit_start <- queryHits(findOverlaps(exon_chr, ASrange, type="start"))
  hit_end <- queryHits(findOverlaps(exon_chr, ASrange, type="end"))
  hit_equal <- queryHits(findOverlaps(exon_chr, ASrange, type="equal"))
  
  if(length(hit_within) > 0){
    within_all <- data.frame(OverType = "within", exon_ca[[1]][hit_within,])
  }else{
    within_all <- data.frame(OverType = "within", t(data.frame(rep(NA,9), row.names = paste("V", 1:9, sep = ""))))
  }
  
  if(length(hit_start) >0 ){
    start_all <- data.frame(OverType = "start", exon_ca[[1]][hit_start,])
  }else{
    start_all <- data.frame(OverType = "start", t(data.frame(rep(NA,9), row.names = paste("V", 1:9, sep = ""))))
  }
  
  if(length(hit_end) > 0){
    end_all <- data.frame(OverType = "end", exon_ca[[1]][hit_end,])
  }else{
    end_all <- data.frame(OverType = "end", t(data.frame(rep(NA,9), row.names = paste("V", 1:9, sep = ""))))
  }
  
  if(length(hit_equal) > 0){
    equal_all <- data.frame(OverType = "equal", exon_ca[[1]][hit_equal,])
  }else{
    equal_all <- data.frame(OverType = "equal", t(data.frame(rep(NA,9), row.names = paste("V", 1:9, sep = ""))))
  }
  
  se_Add_gtf_chr[[i]] <- suppressWarnings(cbind(SE_type[i,], rbind(within_all,start_all,end_all,equal_all)))
})

system.time(se_Add_gtf_chr1 <- do.call(rbind, se_Add_gtf_chr))
system.time(se_gtf_chr1 <- bind_rows(se_Add_gtf_chr))
# [1] 1208   21


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
# [1] 3789   10
dim(SE_type)
# [1] 3789   14
dim(SE_GTF)
# [1] 21104    24


#SE_type <- SE_type[,-c("PSI.R", "PSI.L", "PSI")]
SE_type <- SE_type[,-(8:10)]

## merge PSI and SE information
for(i in colnames(psi.co)){
  SE_type[,paste(i, ".L_psi", sep = "")] <- psi.co[as.character(SE_type$right),i]
  SE_type[,paste(i, ".R_psi", sep = "")] <- psi.co[as.character(SE_type$left),i]
  SE_type[,i] <- rowMeans(cbind(psi.co[as.character(SE_type$left),i], psi.co[as.character(SE_type$right),i]), na.rm = T)
}# as.vector() is necessary


save(SE_type, SE_GTF, file = paste(pffo, "SE/RData/SE_DSU_10cells_result_20180223.RData", sep = ""))



#### 7. find A3SS -------------------------------------------------------------------------------------------------------------------------------



head(junc.as.2)
head(junc.as.2.ot)

junc.as.2.ot_uniq <- junc.as.2.ot[!duplicated(junc.as.2.ot),]
## 去除重复是因为要将所有的 AS SJ 都用起来，主要目的是搜集那些在 SE 的时候因为是出现 3 次但其实并不是 SE 的那些 SJ, 
## 它们之前出现三次，我们认为有可能是 SE ， 但是之后因为 AS insert 等原因判断出并不是 SE，但也出现了三次，
## 其中还有一次是重复的，所以把重复去掉。

sum(table(junc.as.2.ot_uniq$start)==2)

## 把真正的 SE 去掉：
head(junc.as.2.ot_uniq)
junc.as.2.ot_uniq$event <- do.call(rbind, strsplit(rownames(junc.as.2.ot_uniq), split = "\\."))[,1]

a3_raw <- junc.as.2.ot_uniq[junc.as.2.ot_uniq$start %in% names(table(junc.as.2.ot_uniq$start)[table(junc.as.2.ot_uniq$start)==2]), ]

sum(a3_raw$event %in% SE$left)
sum(a3_raw$event %in% SE$left == FALSE)

a3_raw <- a3_raw[a3_raw$event %in% SE$left == FALSE, ]
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
#                                                           1976
#                                          exon-skipping_exactly
#                                                           1998
#                                    exon-skipping_multiple-exon
#                                                           3312
#exon-skipping_multiple-exon_and_first_and_last_exon_unannotated
#                                                            194
#         exon-skipping_multiple-exon_and_first_exon_unannotated
#                                                            130
#                     exon-skipping_skipped_one_unannotated_exon
#                                                            440
#                               intergenic-alt_5/3-splicing_site
#                                                            178
#                                                         Others
#                                                            753


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
#[1] 8981    6
dim(A3SS_psi)
#[1] 8981 3580
dim(A3SS_type)
#[1] 8981   11
dim(A3SS_GTF)
#[1] 69880    21


save(A3SS, A3SS_psi, A3SS_type, A3SS_GTF, file = paste(pffo, "A3SS/RData/A3SS_DSU_10cells_result_20180223.RData", sep = ""))





#### 10. find A5SS -------------------------------------------------------------------------------------------------------------------------------




head(junc.as.2)
head(junc.as.2.ot)

junc.as.2.ot_uniq <- junc.as.2.ot[!duplicated(junc.as.2.ot, fromLast = TRUE),]

sum(table(junc.as.2.ot_uniq$end)==2)

head(junc.as.2.ot_uniq)
junc.as.2.ot_uniq$event <- do.call(rbind, strsplit(rownames(junc.as.2.ot_uniq), split = "\\."))[,1]

## 把真正的 SE 去掉：

a5_raw <- junc.as.2.ot_uniq[junc.as.2.ot_uniq$end %in% names(table(junc.as.2.ot_uniq$end)[table(junc.as.2.ot_uniq$end)==2]), ]

a5_raw <- a5_raw[order(a5_raw$chr, a5_raw$end, a5_raw$start),]


sum(a5_raw$event %in% SE$right)


a5_raw <- a5_raw[a5_raw$event %in% SE$right == FALSE, ]
## 把真正的 SE 去掉

a5_raw <- dlply(a5_raw, c("chr","end"))
sum(unlist(lapply(a5_raw,function(x) length(unique(x$chr)) != 1)))
# if != 0, there are some problem!!!


names(a5_raw) <- sapply(a5_raw,function(x) paste(x$names,collapse="_"))

names(a5_raw)
as.vector(sapply(a5_raw,function(x) diff(as.numeric(x$start))))

cbind(names(a5_raw),as.vector(sapply(a5_raw,function(x) diff(as.numeric(x$start)))))

sapply(a5_raw, function(x) min(as.numeric(x$start)))
sapply(a5_raw, function(x) max(as.numeric(x$start)))
sapply(a5_raw, function(x) mean(as.numeric(x$end)))
sapply(a5_raw, function(x) unique(x$chr))

a5 <- data.frame(name = names(a5_raw), row.names = names(a5_raw),
                 distan = as.integer(sapply(a5_raw,function(x) diff(as.numeric(x$start)))), 
                 chr = as.character(sapply(a5_raw, function(x) unique(x$chr))),
                 start1 = sapply(a5_raw, function(x) min(as.numeric(x$start))),
                 start2 = sapply(a5_raw, function(x) max(as.numeric(x$start))),
                 end = sapply(a5_raw, function(x) mean(as.numeric(x$end))), stringsAsFactors = F)


A5SS <- a5[rownames(a5) %in% rownames(psi.co),]

identical(rownames(psi.co[row.names(A5SS),]),rownames(A5SS))
identical(row.names(psi.co)[which(row.names(psi.co) %in% row.names(A5SS))],rownames(A5SS))

A5SS_psi <- cbind(A5SS, psi.co[as.character(row.names(A5SS)),])

A5SS_psi <- A5SS_psi[order(A5SS_psi$chr,A5SS_psi$end,A5SS_psi$start1,A5SS_psi$start2),]

A5SS <- A5SS[order(A5SS$chr,A5SS$end,A5SS$start1,A5SS$start2),]

plot(density(A5SS$distan))

vioplot::vioplot(A5SS$distan[A5SS$distan<300])

dim(SE)
dim(A3SS)
dim(A5SS)




#### 11.validate A5SS using GTF --------------------------------------------------------




a5_od <- A5SS_psi


unique(a5_od$chr)
chr=unique(a5_od$chr)
A5SS_type <- dlply(a5_od[,1:7], .variables = "chr")


for(j in 1:length(chr)){
  print(chr[j])
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
    print(paste(i," chr", chr[j], ":", test[i,3],"-",test[i,4]," ", AStype," - ",ASlen,sep=""))
  }
}


A5SS_type <- do.call(rbind, A5SS_type)

head(A5SS_type)
dim(A5SS_type)

table(A5SS_type$AStype)
#                                          alt_5/3-splicing_site
#                                                           1916
#                                                  exon-skipping
#                                                           2070
#                                    exon-skipping_multiple-exon
#                                                           3418
#exon-skipping_multiple-exon_and_first_and_last_exon_unannotated
#                                                            218
#          exon-skipping_multiple-exon_and_last_exon_unannotated
#                                                            127
#                     exon-skipping_skipped_one_unannotated_exon
#                                                            496
#                               intergenic-alt_5/3-splicing_site
#                                                            165
#                                                         Others
#                                                            767

dim(A5SS)
dim(A5SS_psi)
dim(A5SS_type)





#### 9.Find the A5SS exon and original gene ----------------------------------------------------------------------------------------------------------





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
})
#     user   system  elapsed 
# 2748.168   22.768 2768.673 

system.time(A5SS_GTF <- bind_rows(A5SS_gtf_chr))

dim(A5SS)
#[1] 9177    6
dim(A5SS_psi)
#[1] 9177 3580
dim(A5SS_type)
#[1] 9177   11
dim(A5SS_GTF)
#[1] 71900    21


save(A5SS, A5SS_psi, A5SS_type, A5SS_GTF, file = paste(pffo, "A5SS/RData/A5SS_DSU_10cells_result_20180223.RData", sep = ""))





#### 13.find MXE ---------------------------------------------------------------------------------





junc.as.2.ot

sum(table(junc.as.2.ot$start)==2)
sum(table(junc.as.2.ot$end)==2)

start_pos <- names(table(junc.as.2.ot$start)[table(junc.as.2.ot$start)==2])

end_pos <- names(table(junc.as.2.ot$end)[table(junc.as.2.ot$end)==2])

junc.as.2.ot[junc.as.2.ot$start %in% start_pos | junc.as.2.ot$end %in% end_pos, ]

## start1 < end2 < start2
## end1 < start1 < end2
## start < end1 < start1 < end2 < start2 < end

head(a3[order(a3$chr, a3$start, a3$end1, a3$end2),])
head(a5[order(a5$chr, a5$start1, a5$start2, a5$end),])


merge.data.frame(x = a5, y = a3, by = "name", all = T)

merge(x = data.table(a3, key = "name"),y = data.table(a5, key = "name"),all = T, by = "name")

a3$name %in% a5$name

a3_for_mxe <- data.frame(a3[order(a3$chr, a3$start, a3$end1, a3$end2),])
a5_for_mxe <- data.frame(a5[order(a5$chr, a5$start1, a5$start2, a5$end),])

## 思路是将 A3 与 A5 每两个都进行两两组合，然后根据位置（start < end1 < start1 < end2 < start2 < end）判断是否为 MXE
library(dplyr)

a5_for_mxe_rep <- bind_rows(replicate(nrow(a3_for_mxe), a5_for_mxe, simplify = FALSE))
dim(a5_for_mxe_rep)
head(a5_for_mxe_rep)

library(tidyr)

a3_for_mxe_rep <- bind_rows(replicate(nrow(a5_for_mxe), a3_for_mxe, simplify = FALSE))

a3_for_mxe_rep <- data.frame(a3_for_mxe_rep[order(a3_for_mxe_rep$chr, a3_for_mxe_rep$start, a3_for_mxe_rep$end1, a3_for_mxe_rep$end2),])




mxe_raw <- cbind(a3_for_mxe_rep,a5_for_mxe_rep)


dim(mxe_raw)
colnames(mxe_raw)



MXE <- mxe_raw[mxe_raw[,"start1"] < mxe_raw[,"end2"] & mxe_raw[,"start1"] > mxe_raw[,"end1"] &
                 mxe_raw[,"end2"] < mxe_raw[,"start2"] & mxe_raw[,"end2"] > mxe_raw[,"start1"], ]


## we just neeed the A3 and A5 from the same chr
MXE <- MXE[apply(MXE[,colnames(MXE) == "chr"], 1, 
                         function(x) sum(duplicated(as.character(x))) > 0),]



#### 14.validate MXE using GTF ---------------------------------------------------------------------------------



unique(MXE$chr)
chr=unique(MXE$chr)

MXE_type <- dlply(MXE, .variables = "chr")

for(j in 1:length(chr)){
  print(chr[j])
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
    print(paste(i," chr", chr[j], ":", test[i,5],"-",test[i,4]," ", AStype," - ",ASlen,sep=""))
  }
}


MXE_type <- do.call(rbind, MXE_type)

head(MXE_type)
dim(MXE_type)

table(MXE_type$AStype)
#         Mixed-splicing Mutually-exclusive-exon
#                    599                     278

dim(MXE)
dim(MXE_type)

sum(MXE_type$ASrange %in% paste("chr", rownames(count_use),sep = "") == FALSE & MXE_type$AStype == "Mutually-exclusive-exon")

## the SJ of start1 and end2 are not exist
## This is important
MXE_type <- MXE_type[MXE_type$ASrange %in% paste("chr", rownames(count_use),sep = "") == FALSE, ]
dim(MXE_type)

table(MXE_type$AStype)
#         Mixed-splicing Mutually-exclusive-exon
#                    368                     103



save(MXE,MXE_type, file = paste(pffo, "MXE/RData/MXE_DSU_10cells_result_20171223.RData", sep = ""))
























