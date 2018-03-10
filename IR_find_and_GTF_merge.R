require(plyr)
require(multicore)
require(VGAM)
library("data.table")

setwd("/mnt/data5/BGI/UCB/tangchao/IR/Result")


####**** 1.merge all cells' data ****-------------------------------------------------------------------------------------

#### median 
#    merge multiple tables from a directory
multmerge = function(mypath){
  #need to change the pattern accordingly !!!
  filenames=list.files(path=mypath, pattern="txt", full.names=TRUE)
  datalist = lapply(filenames, function(x){
    tmp<-fread(x,header=TRUE, select = 4, sep = "\t")
    setnames(tmp, "median",tail(strsplit(x,"/")[[1]],n=1))
    return(tmp)})
  Reduce(function(x,y) {cbind(x,y)}, datalist)
}

path<-file.path("/mnt/data5/BGI/UCB/tangchao/IR/Result/")

system.time(median_merge<-multmerge(path))


coordinates <- as.data.frame(fread(paste(path, "UCB1.01441.txt", sep = ""),select = 1:3, header =T, sep = "\t"))
event_id <- paste(coordinates$chr,":", coordinates$start, "-", coordinates$end, sep="")

IR_median <- data.frame(as.data.frame(median_merge), row.names = event_id)
colnames(IR_median) <- substr(colnames(IR_median),1,10)

save(IR_median, file = "/mnt/data5/BGI/UCB/tangchao/IR/median/IR_median.RData")
write.table(IR_median, "/mnt/data5/BGI/UCB/tangchao/IR/median/IR_median.txt", sep = "\t", quote=F)


#### mad

multmerge = function(mypath){
  #need to change the pattern accordingly !!!
  filenames=list.files(path=mypath, pattern="txt", full.names=TRUE)
  datalist = lapply(filenames, function(x){
    tmp<-fread(x,header=TRUE, select = 5, sep = "\t")
    setnames(tmp, "mad",tail(strsplit(x,"/")[[1]],n=1))
    return(tmp)})
  Reduce(function(x,y) {cbind(x,y)}, datalist)
}


system.time(mad_merge<-multmerge(path))
#     user    system   elapsed
# 4655.576  5546.676 10307.966


IR_mad <- data.frame(as.data.frame(mad_merge), row.names = event_id)
colnames(IR_mad) <- substr(colnames(IR_mad),1,10)

save(IR_mad, file = "/mnt/data5/BGI/UCB/tangchao/IR/mad/IR_mad.RData")
write.table(IR_mad, "/mnt/data5/BGI/UCB/tangchao/IR/mad/IR_mad.txt", sep = "\t", quote=F)


#### coverage

multmerge = function(mypath){
  #need to change the pattern accordingly !!!
  filenames=list.files(path=mypath, pattern="txt", full.names=TRUE)
  datalist = lapply(filenames, function(x){
    tmp<-fread(x,header=TRUE, select = 6, sep = "\t")
    setnames(tmp, "coverage",tail(strsplit(x,"/")[[1]],n=1))
    return(tmp)})
  Reduce(function(x,y) {cbind(x,y)}, datalist)
}


system.time(coverage_merge<-multmerge(path))
#     user    system   elapsed
# 4655.576  5546.676 10307.966


IR_coverage <- data.frame(as.data.frame(coverage_merge), row.names = event_id)
colnames(IR_coverage) <- substr(colnames(IR_coverage),1,10)

save(IR_coverage, file = "/mnt/data5/BGI/UCB/tangchao/IR/coverage/IR_coverage.RData")
write.table(IR_coverage, "/mnt/data5/BGI/UCB/tangchao/IR/coverage/IR_coverage.txt", sep = "\t", quote=F)


#### IRL

multmerge = function(mypath){
  #need to change the pattern accordingly !!!
  filenames=list.files(path=mypath, pattern="txt", full.names=TRUE)
  datalist = lapply(filenames, function(x){
    tmp<-fread(x,header=TRUE, select = 9, sep = "\t")
    setnames(tmp, "irl",tail(strsplit(x,"/")[[1]],n=1))
    return(tmp)})
  Reduce(function(x,y) {cbind(x,y)}, datalist)
}


system.time(irl_merge<-multmerge(path))
#     user    system   elapsed
# 4655.576  5546.676 10307.966


IR_irl <- data.frame(as.data.frame(irl_merge), row.names = event_id)
colnames(IR_irl) <- substr(colnames(IR_irl),1,10)

save(IR_irl, file = "/mnt/data5/BGI/UCB/tangchao/IR/IRL/IR_irl.RData")
write.table(IR_irl, "/mnt/data5/BGI/UCB/tangchao/IR/IRL/IR_irl.txt", sep = "\t", quote=F)



#### merge old data IR and new data IR ####

IR_median -> IR_median_new
IR_mad -> IR_mad_new
IR_coverage -> IR_coverage_new
IR_irl -> IR_irl_new

load("/mnt/data5/BGI/UCB/tangchao/IR/IR_before/median/IR_median.RData")
load("/mnt/data5/BGI/UCB/tangchao/IR/IR_before/mad/IR_mad.RData")
load("/mnt/data5/BGI/UCB/tangchao/IR/IR_before/coverage/IR_coverage.RData")
load("/mnt/data5/BGI/UCB/tangchao/IR/IR_before/IRL/IR_irl.RData")


## coverage
identical(row.names(IR_coverage_new),row.names(IR_coverage))
IR_coverage <- cbind(IR_coverage, IR_coverage_new)

colnames(IR_coverage) <- gsub(colnames(IR_coverage), pattern="UCB3", replacement="UCB4")
colnames(IR_coverage) <- gsub(colnames(IR_coverage), pattern="UCB1", replacement="UCB3")
colnames(IR_coverage) <- gsub(colnames(IR_coverage), pattern="UCB5", replacement="UCB1")

read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Basic_stat.txt", sep = "\t", stringsAsFactors=F, header=T)[,1]
## only use the cells that BGI used
IR_coverage <- IR_coverage[,read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Basic_stat.txt", sep = "\t", stringsAsFactors=F, header=T)[,1]]

save(IR_coverage, file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge/IR_coverage.RData")
write.table(IR_coverage, "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge/IR_coverage.txt", sep = "\t", quote=F, row.names = F)


## irl
IR_irl <- cbind(IR_irl, IR_irl_new)

colnames(IR_irl) <- gsub(colnames(IR_irl), pattern="UCB3", replacement="UCB4")
colnames(IR_irl) <- gsub(colnames(IR_irl), pattern="UCB1", replacement="UCB3")
colnames(IR_irl) <- gsub(colnames(IR_irl), pattern="UCB5", replacement="UCB1")

read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Basic_stat.txt", sep = "\t", stringsAsFactors=F, header=T)[,1]
## only use the cells that BGI used
IR_irl <- IR_irl[,read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Basic_stat.txt", sep = "\t", stringsAsFactors=F, header=T)[,1]]

save(IR_irl, file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge/IR_irl.RData")
write.table(IR_irl, "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge/IR_irl.txt", sep = "\t", quote=F, row.names = F)


## median
IR_median <- cbind(IR_median, IR_median_new)

colnames(IR_median) <- gsub(colnames(IR_median), pattern="UCB3", replacement="UCB4")
colnames(IR_median) <- gsub(colnames(IR_median), pattern="UCB1", replacement="UCB3")
colnames(IR_median) <- gsub(colnames(IR_median), pattern="UCB5", replacement="UCB1")

read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Basic_stat.txt", sep = "\t", stringsAsFactors=F, header=T)[,1]
## only use the cells that BGI used
IR_median <- IR_median[,read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Basic_stat.txt", sep = "\t", stringsAsFactors=F, header=T)[,1]]

save(IR_median, file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge/IR_median.RData")
write.table(IR_median, "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge/IR_median.txt", sep = "\t", quote=F, row.names = F)


## median
IR_mad <- cbind(IR_mad, IR_mad_new)

colnames(IR_mad) <- gsub(colnames(IR_mad), pattern="UCB3", replacement="UCB4")
colnames(IR_mad) <- gsub(colnames(IR_mad), pattern="UCB1", replacement="UCB3")
colnames(IR_mad) <- gsub(colnames(IR_mad), pattern="UCB5", replacement="UCB1")

read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Basic_stat.txt", sep = "\t", stringsAsFactors=F, header=T)[,1]
## only use the cells that BGI used
IR_mad <- IR_mad[,read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Basic_stat.txt", sep = "\t", stringsAsFactors=F, header=T)[,1]]

save(IR_mad, file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge/IR_mad.RData")
write.table(IR_mad, "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge/IR_mad.txt", sep = "\t", quote=F, row.names = F)






#### IRL
irl_rowSum <- rowSums(IR_irl, na.rm=T)

irl_NASum <- rowSums(is.na(IR_irl))

irl_dataSum <- rowSums(!is.na(IR_irl) & IR_irl>0)


sum(irl_dataSum>=10)
# 9880


coverage_Sum <- rowSums(IR_coverage == 100)

sum(coverage_Sum>=10)
# 7010

sub_irl <- IR_irl[intersect(names(coverage_Sum[coverage_Sum>=10]),names(irl_dataSum[irl_dataSum>=10])),]
dim(sub_irl)
#[1] 7010 3574




#### 5.validate IR using GTF file ---------------------------------------------------
library(GenomicRanges)
library(dplyr)
library(plyr)
library(pryr)




gtf <- read.table("/mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.87.gtf",head=F,sep="\t")

exon <- subset(gtf,gtf$V3=="exon")

exon <- exon[!duplicated(exon[,c(1,4,5)]),]

exon_ca <- dlply(exon, "V1")

exon_c <- exon[,c(1,4,5)]

mode(exon_c)

exon_cn <- split(IRanges(exon_c[,2], exon_c[,3]), exon_c[,1])



IR <- data.frame(chr = do.call(rbind,strsplit(row.names(sub_irl), split = "[:-]"))[,1], 
  start = do.call(rbind,strsplit(row.names(sub_irl), split = "[:-]"))[,2],
  end = do.call(rbind,strsplit(row.names(sub_irl), split = "[:-]"))[,3], event_id = row.names(sub_irl), stringsAsFactors = F)

IR <- IR[IR$chr %in% c(1:22,"X","Y"),]


unique(IR$chr)
chr <- unique(IR$chr)

IR_type <- dlply(IR, .variables = "chr")

for(j in 1:length(chr)){
  print(chr[j])
  eval(parse(text = paste("exon_cn$'",chr[j],"'->exon_chr",sep="")))
  temp <- IR[IR$chr == chr[j],]
  if(nrow(temp) > 1){
    test <- as.matrix(temp[,1:3])
  }else{
    test <- t(as.matrix(temp[,1:3]))
  }
  
  mode(test) <- "numeric"
  for(i in 1:nrow(test)){
    ASrange <- IRanges(test[i,2], test[i,3])
    ## because the positions of SJ are intron sets, so we  +- 1 to match exon sets
    ASlen <- diff(c(test[i,2], test[i,3]))
    ov_within<-countOverlaps(exon_chr, ASrange, type="within")
    ov_start<-countOverlaps(exon_chr, ASrange, type="start")
    ov_end<-countOverlaps(exon_chr, ASrange, type="end")
    ov_equal<-countOverlaps(exon_chr, ASrange, type="equal")
    ov_any<-countOverlaps(exon_chr,ASrange, type="any")
    if (sum(ov_any)==0) {
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
#             1099             5910

IR_type_raw <- IR_type

IR_type <- IR_type[IR_type$AStype == "intron_retention", ]





#### orignal gene finding ####


## all chr


chr <- unique(IR_type$chr)

IR_gtf_chr <- list()

system.time(for(j in 1:length(chr)){
  print(paste("chr",chr[j]))
  eval(parse(text = paste("exon_cn$'",chr[j],"'->exon_chr",sep="")))
  temp <- IR_type[IR_type$chr == chr[j],]
  
  eval(parse(text = paste("exon_ca$'",chr[j],"'->exon_ca_chr",sep="")))
  
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
    
    hit_any <- queryHits(findOverlaps(exon_chr, ASrange, type="any"))
    
    if(length(hit_any) > 0){
      within_all <- data.frame(OverType = "any", exon_ca_chr[hit_any,])
    }
    
    ir_gtf_chr[[i]] <- suppressWarnings(cbind(temp[i,], within_all))
    
  }
  
  IR_gtf_chr[[j]] <- suppressMessages(bind_rows(ir_gtf_chr))
})
#     user   system  elapsed 
# 4079.416   40.844 4119.096 

system.time(IR_GTF <- bind_rows(IR_gtf_chr))

dim(IR)
#[1] 7009    4
dim(IR_type)
#[1] 1099    8
dim(IR_GTF)
#[1] 4008   18


save(IR, IR_type, IR_GTF, file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_GTF_match.RData")






