
##  Merge the IR results
setwd("/mnt/data5/BGI/UCB/tangchao/IR/Result")
library("data.table")

multmerge = function(mypath){
  #need to change the pattern accordingly !!!
  filenames = list.files(path=mypath, pattern="txt", full.names=TRUE)
  datalist = lapply(filenames, function(x){
    tmp <- fread(x,header=TRUE)
    tmp <- as.data.frame(tmp)
    tmp$name <- paste(tmp[,1],":",tmp[,2],"-",tmp[,3],sep="")
    tmp <- tmp[,c("name","median", "mad", "coverage", "median_left", "median_right", "irl")]
    colnames(tmp)[2:ncol(tmp)] <- paste(substr(tail(strsplit(x,"/")[[1]],n=1),1,10),colnames(tmp)[2:ncol(tmp)], sep = "_")
    tmp <- as.data.table(tmp)
    setkey(tmp,name)
    return(tmp)})
  Reduce(function(x,y) {merge(x,y,all=T,by="name")}, datalist)
}
mypath<-file.path("/mnt/data5/BGI/UCB/tangchao/IR/Result")

# system.time(te<-multmerge(mypath)) # killed after 
# 

filenames = list.files(path=mypath, pattern="txt", full.names=TRUE)
for(i in 1:length(filenames)){
	if(i%%100 == 0) print(paste(i, "of", length(filenames), "--", date()))
	if(i == 1){
		tmp <- fread(filenames[i],header=TRUE)
    	tmp <- as.data.frame(tmp)
    	tmp$name <- paste(tmp[,1],":",tmp[,2],"-",tmp[,3],sep="")
    	tmp <- tmp[,c("name","median", "mad", "coverage", "median_left", "median_right", "irl")]
    	colnames(tmp)[2:ncol(tmp)] <- paste(substr(tail(strsplit(filenames[i],"/")[[1]],n=1),1,10),colnames(tmp)[2:ncol(tmp)], sep = "_")
    	tmp <- as.data.table(tmp)
    	setkey(tmp,name)
    	te <- tmp
	}else{
		tmp <- fread(filenames[i],header=TRUE)
    	tmp <- as.data.frame(tmp)
    	tmp$name <- paste(tmp[,1],":",tmp[,2],"-",tmp[,3],sep="")
    	tmp <- tmp[,c("name","median", "mad", "coverage", "median_left", "median_right", "irl")]
    	colnames(tmp)[2:ncol(tmp)] <- paste(substr(tail(strsplit(filenames[i],"/")[[1]],n=1),1,10),colnames(tmp)[2:ncol(tmp)], sep = "_")
    	tmp <- as.data.table(tmp)
    	setkey(tmp,name)
    	te <- merge(x = te, y = tmp, all=T, by="name")
	}
	gc()	
}


save(te, file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/te_1557.RData")


length(filenames)
# [1] 3686
filenames2 = list.files(path=mypath, pattern="txt", full.names=TRUE)[1558:2557]
for(i in 1:length(filenames2)){
	if(i%%100 == 0) print(paste(i, "of", length(filenames2), "--", date()))
	if(i == 1){
		tmp <- fread(filenames2[i],header=TRUE)
    	tmp <- as.data.frame(tmp)
    	tmp$name <- paste(tmp[,1],":",tmp[,2],"-",tmp[,3],sep="")
    	tmp <- tmp[,c("name","median", "mad", "coverage", "median_left", "median_right", "irl")]
    	colnames(tmp)[2:ncol(tmp)] <- paste(substr(tail(strsplit(filenames2[i],"/")[[1]],n=1),1,10),colnames(tmp)[2:ncol(tmp)], sep = "_")
    	tmp <- as.data.table(tmp)
    	setkey(tmp,name)
    	te2 <- tmp
	}else{
		tmp <- fread(filenames2[i],header=TRUE)
    	tmp <- as.data.frame(tmp)
    	tmp$name <- paste(tmp[,1],":",tmp[,2],"-",tmp[,3],sep="")
    	tmp <- tmp[,c("name","median", "mad", "coverage", "median_left", "median_right", "irl")]
    	colnames(tmp)[2:ncol(tmp)] <- paste(substr(tail(strsplit(filenames2[i],"/")[[1]],n=1),1,10),colnames(tmp)[2:ncol(tmp)], sep = "_")
    	tmp <- as.data.table(tmp)
    	setkey(tmp,name)
    	te2 <- merge(x = te2, y = tmp, all=T, by="name")
	}
	gc()	
}


save(te2, file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/te_1558_2557.RData")




filenames3 = list.files(path=mypath, pattern="txt", full.names=TRUE)[2558:3686]
for(i in 1:length(filenames3)){
	if(i%%100 == 0) print(paste(i, "of", length(filenames3), "--", date()))
	if(i == 1){
		tmp <- fread(filenames3[i],header=TRUE)
    	tmp <- as.data.frame(tmp)
    	tmp$name <- paste(tmp[,1],":",tmp[,2],"-",tmp[,3],sep="")
    	tmp <- tmp[,c("name","median", "mad", "coverage", "median_left", "median_right", "irl")]
    	colnames(tmp)[2:ncol(tmp)] <- paste(substr(tail(strsplit(filenames3[i],"/")[[1]],n=1),1,10),colnames(tmp)[2:ncol(tmp)], sep = "_")
    	tmp <- as.data.table(tmp)
    	setkey(tmp,name)
    	te3 <- tmp
	}else{
		tmp <- fread(filenames3[i],header=TRUE)
    	tmp <- as.data.frame(tmp)
    	tmp$name <- paste(tmp[,1],":",tmp[,2],"-",tmp[,3],sep="")
    	tmp <- tmp[,c("name","median", "mad", "coverage", "median_left", "median_right", "irl")]
    	colnames(tmp)[2:ncol(tmp)] <- paste(substr(tail(strsplit(filenames3[i],"/")[[1]],n=1),1,10),colnames(tmp)[2:ncol(tmp)], sep = "_")
    	tmp <- as.data.table(tmp)
    	setkey(tmp,name)
    	te3 <- merge(x = te3, y = tmp, all=T, by="name")
	}
	gc()	
}

save(te3, file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/te_2558_3686.RData")

load("/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/te_1557.RData")
load("/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/te_1558_2557.RData")
load("/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/te_2558_3686.RData")


dim(te)
#[1] 562452   9343
dim(te2)
#[1] 354484   6001
dim(te3)
#[1] 359206   6775

## I have used too many ways to merge the 3 tables, but all ways have failed.
## So, I make a cutoff before merge the tables
load("/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/te_1557.RData")
te <- data.frame(te[,2:ncol(te)], row.names=te[,1])

te_intron_cov <- apply(te[, grep("coverage", colnames(te))],1, function(x) sum(x == 100, na.rm = T))
sum(te_intron_cov>=1)
# [1] 22093
te <- te[te_intron_cov>=1, ]

load("/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/te_1558_2557.RData")
te2 <- data.frame(te2[,2:ncol(te2)], row.names=te2[,1])
te2_intron_cov <- apply(te2[, grep("coverage", colnames(te2))],1, function(x) sum(x == 100, na.rm = T))
sum(te2_intron_cov>=1)
# [1] 10725
te2 <- te2[te2_intron_cov>=1, ]


load("/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/te_2558_3686.RData")
te3 <- data.frame(te3[,2:ncol(te3)], row.names=te3[,1])
te3_intron_cov <- apply(te3[, grep("coverage", colnames(te3))],1, function(x) sum(x == 100, na.rm = T))
sum(te3_intron_cov>=1)
# [1] 13412
te3 <- te3[te3_intron_cov>=1, ]


save(te,te2,te3, file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/te_for_merge.RData")
load("/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/te_for_merge.RData")
#te12 <- merge(te,te2,all=T)
dim(te)
dim(te2)
dim(te3)
length(unique(c(row.names(te),row.names(te2),row.names(te3))))
# [1] 34125
#te_irl <- te[, grep("irl", colnames(te))]
#te2_irl <- te2[, grep("irl", colnames(te2))]
#te12_irl <- merge(te_irl,te2_irl,all=T)

## Even though I did the filtering before the merge, I still couldn't.
## So, I used a stupid way:


test1 <- te[row.names(te) %in% row.names(te2),]
test2 <- te2[row.names(te2) %in% row.names(te),]
identical(row.names(test1), row.names(test2))
# [1] TRUE
test12 <- cbind(test1,test2)
dim(test12)
# [1] 5313 15342

test1_only <- te[!row.names(te) %in% row.names(te2),]
dim(test1_only)
# [1] 16780  9342
test2_only <- te2[!row.names(te2) %in% row.names(te),]
dim(test2_only)
# [1] 5412 6000



test1_only_in_test2 <- matrix(NA, nrow = nrow(test1_only), ncol = ncol(test2))
colnames(test1_only_in_test2) <- colnames(test2)
row.names(test1_only_in_test2) <- row.names(test1_only)
dim(test1_only_in_test2) # 第二批次中第一批次所特有的
# [1] 16780  6000
test1_only_in_all <- cbind(test1_only, test1_only_in_test2)


test2_only_in_test1 <- matrix(NA, nrow = nrow(test2_only), ncol = ncol(test1))
colnames(test2_only_in_test1) <- colnames(test1)
row.names(test2_only_in_test1) <- row.names(test2_only)
dim(test2_only_in_test1) # 
# [1] 5412 9342
test2_only_in_all <- cbind(test2_only_in_test1, test2_only)

dim(test12)
dim(test1_only_in_all)
dim(test2_only_in_all)
te <- rbind(test12,test1_only_in_all,test2_only_in_all)
dim(te)
# [1] 27505 15342

te2 <- te3 

test1 <- te[row.names(te) %in% row.names(te2),]
test2 <- te2[row.names(te2) %in% row.names(te),]
identical(row.names(test1), row.names(test2))
# [1] FALSE
test2 <- test2[row.names(test1),]
identical(row.names(test1), row.names(test2))
# [1] TRUE
test12 <- cbind(test1,test2)
dim(test12)
# [1]  6792 22116

test1_only <- te[!row.names(te) %in% row.names(te2),]
dim(test1_only)
# [1] 20713 15342
test2_only <- te2[!row.names(te2) %in% row.names(te),]
dim(test2_only)
# [1] 6620 6774



test1_only_in_test2 <- matrix(NA, nrow = nrow(test1_only), ncol = ncol(test2))
colnames(test1_only_in_test2) <- colnames(test2)
row.names(test1_only_in_test2) <- row.names(test1_only)
dim(test1_only_in_test2) # 第二批次中第一批次所特有的
# [1] 20713  6774
test1_only_in_all <- cbind(test1_only, test1_only_in_test2)


test2_only_in_test1 <- matrix(NA, nrow = nrow(test2_only), ncol = ncol(test1))
colnames(test2_only_in_test1) <- colnames(test1)
row.names(test2_only_in_test1) <- row.names(test2_only)
dim(test2_only_in_test1) # 
# [1] 6620 15342
test2_only_in_all <- cbind(test2_only_in_test1, test2_only)

dim(test12)
dim(test1_only_in_all)
dim(test2_only_in_all)

te123 <- rbind(test12,test1_only_in_all,test2_only_in_all)
dim(te123)
# [1] 34125 22116

te123 <- te123[order(row.names(te123)),]
te123 <- te123[order(row.names(te123)),]
save(te123, file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/IR_merged.RData")




