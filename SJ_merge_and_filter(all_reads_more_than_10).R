## SJ filter and merge (all reads >10)

setwd("/mnt/data5/BGI/UCB/tangchao/data/SJ")

library("data.table")

multmerge = function(mypath){
  #need to change the pattern accordingly !!!
  filenames = list.files(path=mypath, pattern="tab", full.names=TRUE)
  datalist = lapply(filenames, function(x){
    tmp <- unique(fread(x,header=FALSE, select = c(1,2,3,7)))
    tmp <- as.data.frame(tmp)
    tmp$name <- paste(tmp[,1],":",tmp[,2],"-",tmp[,3],sep="")
    tmp <- tmp[,c("name","V7")]
    tmp <- tmp[tmp$V7 >= 10,]
    tmp <- as.data.table(tmp)

    setnames(tmp, "V7",tail(strsplit(x,"/")[[1]],n=1))
    setkey(tmp,name)
    return(tmp)})
  Reduce(function(x,y) {merge(x,y,all=T,by="name")}, datalist)
}
path<-file.path("/mnt/data5/BGI/UCB/tangchao/data/SJ/UCB_134")
system.time(te<-multmerge(path))
#      user    system   elapsed
# 149478.08  66937.69  85488.78

dim(te)

save(te, file="./SJ_merged_raw_te(all_more_than_10).RData")


te_table <- as.data.frame(te)

row.names(te_table) <- te_table$name
te_table <- te_table[,-1]

colnames(te_table) <- substr(colnames(te_table), 1, 10)


te_table_rowsum <- apply(te_table,1,function(x) sum(x>0, na.rm=T))
length(te_table_rowsum)
# [1] 461646

sum(te_table_rowsum>=10)
# [1]  94802

sum(te_table_rowsum>=2)
# [1] 197797


## We only care about chr1-22 and X Y
chr_tu <- do.call(rbind, strsplit(rownames(te_table), split = ":"))[,1] %in% c(1:22,"X","Y")
sum(chr_tu)
# [1] 460327

SJ_tu <- te_table[chr_tu,]
dim(SJ_tu)
# [1] 460327   3574

save(SJ_tu, file="./SJ(all_more_than_10)_merged_touse.RData")


#### repeat in at least 2 cells

te_table_rowsum <- apply(te_table,1,function(x) sum(x>0, na.rm=T))
length(te_table_rowsum)
# [1] 461646

sum(te_table_rowsum>=10)
# [1]  94802

sum(te_table_rowsum>=2)
# [1] 197797



te_table <- te_table[te_table_rowsum >= 2, ]

## We only care about chr1-22 and X Y
chr_tu <- do.call(rbind, strsplit(rownames(te_table), split = ":"))[,1] %in% c(1:22,"X","Y")
sum(chr_tu)
# [1] 197427

SJ_tu <- te_table[chr_tu,]
dim(SJ_tu)
# [1] 197427   3574

save(SJ_tu, file="./SJ(all_more_than_10_more_than_1_cells)_merged_touse.RData")


