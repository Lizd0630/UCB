
##  Merge the IR results =================================================================================================================================

setwd("/mnt/data5/BGI/UCB/tangchao/IR/Result")
library("data.table")

mypath<-file.path("/mnt/data5/BGI/UCB/tangchao/IR/Result")
filenames = list.files(path=mypath, pattern="txt", full.names=TRUE)
for(i in 1:length(filenames)){
	if(i%%100 == 0) print(paste(i, "of", length(filenames), "--", date()))
	if(i == 1){
		tmp <- fread(filenames[i],header=TRUE)
    	tmp <- as.data.frame(tmp)
    	tmp$name <- paste(tmp[,1],":",tmp[,2],"-",tmp[,3],sep="")
    	tmp <- tmp[,c("name","median", "mad", "coverage", "median_left", "median_right", "irl")]
    	tmp <- tmp[tmp$coverage == 100 & tmp$median > 2 & tmp$median_right > 2 & tmp$median_left > 2 , ]
    	colnames(tmp)[2:ncol(tmp)] <- paste(substr(tail(strsplit(filenames[i],"/")[[1]],n=1),1,10),colnames(tmp)[2:ncol(tmp)], sep = "_")
    	tmp <- as.data.table(tmp)
    	setkey(tmp,name)
    	te <- tmp
	}else{
		tmp <- fread(filenames[i],header=TRUE)
    	tmp <- as.data.frame(tmp)
    	tmp$name <- paste(tmp[,1],":",tmp[,2],"-",tmp[,3],sep="")
    	tmp <- tmp[,c("name","median", "mad", "coverage", "median_left", "median_right", "irl")]
    	tmp <- tmp[tmp$coverage == 100 & tmp$median > 2 & tmp$median_right > 2 & tmp$median_left > 2 , ]
    	colnames(tmp)[2:ncol(tmp)] <- paste(substr(tail(strsplit(filenames[i],"/")[[1]],n=1),1,10),colnames(tmp)[2:ncol(tmp)], sep = "_")
    	tmp <- as.data.table(tmp)
    	setkey(tmp,name)
    	te <- merge(x = te, y = tmp, all=T, by="name")
	}
	gc()	
}
save(te, file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/te_all_cell.RData")



#### filtering ===========================================================================================================================================




IR <- data.frame(as.data.frame(te)[, 2:ncol(as.data.frame(te))], row.names = as.data.frame(te)[,1])

IR_irl <- IR[, grep("irl", colnames(IR))]
IR_cov <- IR[, grep("coverage", colnames(IR))]
IR_med_l <- IR[, grep("median_left", colnames(IR))]
IR_med_r <- IR[, grep("median_right", colnames(IR))]
IR_mad <- IR[, grep("mad", colnames(IR))]
IR_med <- IR[, grep("median$", colnames(IR))]

colnames(IR_irl) <- substr(colnames(IR_irl),1,10)
colnames(IR_cov) <- substr(colnames(IR_cov),1,10)
colnames(IR_med) <- substr(colnames(IR_med),1,10)

colnames(IR_irl) <- gsub(pattern="UCB3", replacement="UCB4", colnames(IR_irl))
colnames(IR_irl) <- gsub(pattern="UCB1", replacement="UCB3", colnames(IR_irl))
colnames(IR_irl) <- gsub(pattern="UCB5", replacement="UCB1", colnames(IR_irl))

colnames(IR_cov) <- colnames(IR_irl)
colnames(IR_med) <- colnames(IR_irl)
colnames(IR_mad) <- colnames(IR_irl)
colnames(IR_med_l) <- colnames(IR_irl)
colnames(IR_med_r) <- colnames(IR_irl)



cell_tu <- fread("/mnt/data5/BGI/UCB/ExpMat_NewID/Basic_stat.txt", select = 1, sep = "\t")

length(cell_tu[[1]])
# [1] 3574
IR_irl <- IR_irl[,cell_tu[[1]]]
IR_cov <- IR_cov[,cell_tu[[1]]]
IR_med <- IR_med[,cell_tu[[1]]]
IR_mad <- IR_mad[,cell_tu[[1]]]
IR_med_l <- IR_med_l[,cell_tu[[1]]]
IR_med_r <- IR_med_r[,cell_tu[[1]]]




IR_cov_rowsum <- rowSums(!is.na(IR_cov))
sum(IR_cov_rowsum>1)
# [1] 7251
IR_med_rowsum <- apply(IR_med, 1, function(x) sum(x>0, na.rm = T))
sum(IR_med_rowsum>1)
# [1] 7251
IR_med_l_rowsum <- apply(IR_med_l, 1, function(x) sum(x>0, na.rm = T))
sum(IR_med_l_rowsum>1)
# [1] 7159
IR_med_r_rowsum <- apply(IR_med_r, 1, function(x) sum(x>0, na.rm = T))
sum(IR_med_r_rowsum>1)
# [1] 7159

sum(IR_cov_rowsum>1 & IR_med_rowsum>1 & IR_med_l_rowsum>1 & IR_med_r_rowsum>1)
# [1] 7039


IR_irl_tu <- IR_irl[IR_cov_rowsum>1 & IR_med_rowsum>1 & IR_med_l_rowsum>1 & IR_med_r_rowsum>1, ]
dim(IR_irl_tu)
# [1] 7039 3686

IR_irl_rowsum <- rowSums(!is.na(IR_irl_tu))
IR_irl_colsum <- colSums(!is.na(IR_irl_tu))
summary(IR_irl_rowsum)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   2.00    2.00    3.00   15.66    8.00 3574.00
summary(IR_irl_colsum)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   2.00   20.00   27.00   30.84   38.00  143.00



library(reshape2)
IR_med_melt <-  melt(IR_med)
IR_med_melt <-  na.omit(IR_med_melt)
dim(IR_med_melt)
# [1] 127327      2


IR_irl_melt <-  melt(IR_irl)
IR_irl_melt <- na.omit(IR_irl_melt)
summary(IR_irl_melt$value)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#   0.0000    0.5229    0.9430    0.8652    1.0935 1123.1250
save(IR_irl_tu, file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/IR_irl_tu.RData")

write.table(IR_irl_tu, "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/IR_irl_tu.txt", sep = "\t")


IR_cov_tu <- IR_cov[row.names(IR_irl_tu), colnames(IR_irl_tu)]
IR_med_tu <- IR_med[row.names(IR_irl_tu), colnames(IR_irl_tu)]
IR_mad_tu <- IR_mad[row.names(IR_irl_tu), colnames(IR_irl_tu)]
save(IR_cov_tu, file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/IR_cov_tu.RData")
write.table(IR_cov_tu, "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/IR_cov_tu.txt", sep = "\t")
save(IR_med_tu, file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/IR_med_tu.RData")
write.table(IR_med_tu, "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/IR_med_tu.txt", sep = "\t")
save(IR_mad_tu, file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/IR_mad_tu.RData")
write.table(IR_mad_tu, "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/IR_mad_tu.txt", sep = "\t")





