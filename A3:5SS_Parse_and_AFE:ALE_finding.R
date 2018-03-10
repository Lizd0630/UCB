#### modify the GTF to find alternative first/last exon



library(plyr)

setwd("/mnt/data1/projects/UCB/results_ucb/tangchao/IR")
data <- read.table("table_for_intron.txt",sep="\t",header=T)

data_split <- dlply(data, .variables = "transcript_id")

data_split_rev <- lapply(data_split, function(x) {
	x$exon_numble_rev <- -rev(x$exon_numble)
	return(x)
	})

library(dplyr)
data_rev <- do.call(bind_rows, data_split_rev)

save(data_rev, file = '/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/RData/GTF_exon_info_for_alternative_firstlast_exon.RData')


load("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/A3SS/RData/all_cells_A3SS_GTF.RData")

A3SS_GTF[!is.na( A3SS_GTF$V1),] -> test

test$exon_name <- paste(test$V1,":",test$V4,"-",test$V5,sep="")
data_rev$exon_name <- paste(data_rev$chr,":",data_rev$start,"-",data_rev$end,sep="")
test_merge <- merge(x	= test, y = data_rev, by = "exon_name", all.x = T)

dim(test_merge[test_merge$strand == "+" & test_merge$exon_numble_rev == -1,])
# [1] 4568   33
identical(test_merge$strand, test_merge$V7)
# [1] FALSE
dim(test_merge[test_merge$V7 == "+" & test_merge$exon_numble_rev == -1,])
# [1] 4568   33
identical(as.character(test_merge$strand), test_merge$V7)
# [1] TRUE
# dim(test_merge[test_merge$V7 == "+" & test_merge$exon_numble_rev == 1,])
# [1]  0 33

test_merge[test_merge$AStype == "exon-skipping_exactly" & test_merge$V7 == "+" & test_merge$exon_numble_rev == -1,] -> ALE_merge

test_merge[test_merge$AStype == "exon-skipping_exactly" & test_merge$V7 == "-" & test_merge$exon_numble == 1,] -> AFE_merge



load("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/A5SS/RData/all_cells_A5SS_GTF.RData")

A5SS_GTF[!is.na( A5SS_GTF$V1),] -> test

test$exon_name <- paste(test$V1,":",test$V4,"-",test$V5,sep="")
#data_rev$exon_name <- paste(data_rev$chr,":",data_rev$start,"-",data_rev$end,sep="")
test_merge <- merge(x	= test, y = data_rev, by = "exon_name", all.x = T)

dim(test_merge[test_merge$AStype == "exon-skipping" & test_merge$V7 == "+" & test_merge$exon_numble == 1,])
# [1] 2727   33
identical(test_merge$strand, test_merge$V7)
# [1] FALSE
dim(test_merge[test_merge$AStype == "exon-skipping" & test_merge$V7 == "-" & test_merge$exon_numble_rev == -1,])
# [1] 1467   33
identical(as.character(test_merge$strand), test_merge$V7)


test_merge[test_merge$AStype == "exon-skipping" & test_merge$V7 == "+" & test_merge$exon_numble == 1,] -> AFE_merge_2

test_merge[test_merge$AStype == "exon-skipping" & test_merge$V7 == "-" & test_merge$exon_numble_rev == -1,] -> ALE_merge_2


length(unique(ALE_merge$name))
# [1] 624
length(unique(ALE_merge_2$name))
# [1] 587
length(unique(AFE_merge_2$name))
# [1] 1191
length(unique(AFE_merge$name))
# [1] 1098




load("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/A3SS/RData/all_cells_A3SS_type_and_PSI.RData")

sum(unique(ALE_merge$name) %in% all_cells_A3SS_type_and_PSI$name)
# [1] 624

all_cells_A3SS_type_and_PSI[all_cells_A3SS_type_and_PSI$name %in% unique(ALE_merge$name), ] -> all_cells_ALE_type_and_PSI_from_A3SS

all_cells_A3SS_type_and_PSI[all_cells_A3SS_type_and_PSI$name %in% unique(AFE_merge$name), ] -> all_cells_AFE_type_and_PSI_from_A3SS


load("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/A5SS/RData/all_cells_A5SS_type_and_PSI.RData")

all_cells_A5SS_type_and_PSI[all_cells_A5SS_type_and_PSI$name %in% unique(ALE_merge_2$name), ] -> all_cells_ALE_type_and_PSI_from_A5SS

all_cells_A5SS_type_and_PSI[all_cells_A5SS_type_and_PSI$name %in% unique(AFE_merge_2$name), ] -> all_cells_AFE_type_and_PSI_from_A5SS


colnames(all_cells_AFE_type_and_PSI_from_A3SS) <- colnames(all_cells_AFE_type_and_PSI_from_A5SS)
rbind(all_cells_AFE_type_and_PSI_from_A3SS, all_cells_AFE_type_and_PSI_from_A5SS) -> all_cells_AFE_PSI
all_cells_AFE_PSI$AStype <- "AFE"
dim(all_cells_AFE_PSI)
#[1] 2289 3584

colnames(AFE_merge) <- colnames(AFE_merge_2)
rbind(AFE_merge, AFE_merge_2) -> AFE_GTF
AFE_GTF$AStype <- "AFE"

save(all_cells_AFE_PSI, AFE_GTF, file = "/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/AFE/RData/all_cells_AFE_PSI_GTF.RData")


colnames(all_cells_ALE_type_and_PSI_from_A5SS) <- colnames(all_cells_ALE_type_and_PSI_from_A3SS)
rbind(all_cells_ALE_type_and_PSI_from_A3SS, all_cells_ALE_type_and_PSI_from_A5SS) -> all_cells_ALE_PSI
all_cells_ALE_PSI$AStype <- "ALE"
dim(all_cells_ALE_PSI)
#[1] 1211 3584

colnames(ALE_merge_2) <- colnames(ALE_merge)
rbind(ALE_merge, ALE_merge_2) -> ALE_GTF
ALE_GTF$AStype <- "ALE"

save(all_cells_ALE_PSI, ALE_GTF, file = "/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/ALE/RData/all_cells_ALE_PSI_GTF.RData")



#### Parse A3SS/A5SS ==================================================================================================================================================


dim(A3SS_GTF)

A3SS_GTF[A3SS_GTF$AStype == "alt_5/3-splicing_site", ] -> test
test[!is.na(test$V1), ] -> test
unique(test[test$V7 == "+", ]$name) -> A3SS_name_1
unique(test[test$V7 == "-", ]$name) -> A5SS_name_1


A5SS_GTF[A5SS_GTF$AStype == "alt_5/3-splicing_site", ] -> test
test[!is.na(test$V1), ] -> test
unique(test[test$V7 == "+", ]$name) -> A5SS_name_2
unique(test[test$V7 == "-", ]$name) -> A3SS_name_2

all_cells_A3SS_type_and_PSI[all_cells_A3SS_type_and_PSI$name %in% A3SS_name_1, ] -> A3SS_type_and_PSI_1
all_cells_A3SS_type_and_PSI[all_cells_A3SS_type_and_PSI$name %in% A5SS_name_1, ] -> A5SS_type_and_PSI_1


all_cells_A5SS_type_and_PSI[all_cells_A5SS_type_and_PSI$name %in% A3SS_name_2, ] -> A3SS_type_and_PSI_2
all_cells_A5SS_type_and_PSI[all_cells_A5SS_type_and_PSI$name %in% A5SS_name_2, ] -> A5SS_type_and_PSI_2


colnames(A3SS_type_and_PSI_2) <- colnames(A3SS_type_and_PSI_1)
rbind(A3SS_type_and_PSI_1, A3SS_type_and_PSI_2) -> A3SS_type_and_PSI
dim(A3SS_type_and_PSI)
# [1] 5874 3584


colnames(A5SS_type_and_PSI_1) <- colnames(A5SS_type_and_PSI_2)
rbind(A5SS_type_and_PSI_1, A5SS_type_and_PSI_2) -> A5SS_type_and_PSI
dim(A5SS_type_and_PSI)
# [1] 3908 3584


A3SS_GTF[A3SS_GTF$name %in% A3SS_name_1, ] -> A3SS_GTF_1
A3SS_GTF[A3SS_GTF$name %in% A5SS_name_1, ] -> A5SS_GTF_1

A5SS_GTF[A5SS_GTF$name %in% A3SS_name_2, ] -> A3SS_GTF_2
A5SS_GTF[A5SS_GTF$name %in% A5SS_name_2, ] -> A5SS_GTF_2

colnames(A3SS_GTF_2) <- colnames(A3SS_GTF_1)
rbind(A3SS_GTF_2, A3SS_GTF_1) -> A3SS_GTF

colnames(A5SS_GTF_1) <- colnames(A5SS_GTF_2)
rbind(A5SS_GTF_1, A5SS_GTF_2) -> A5SS_GTF

A3SS_type_and_PSI$AStype <- "A3SS"
A3SS_GTF$AStype <- "A3SS"
save(A3SS_type_and_PSI, A3SS_GTF, file = "/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/A3SS/RData/Parsed_A3SS_PSI_GTF.RData")


A5SS_type_and_PSI$AStype <- "A3SS"
A5SS_GTF$AStype <- "A3SS"
save(A5SS_type_and_PSI, A5SS_GTF, file = "/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/A5SS/RData/Parsed_A5SS_PSI_GTF.RData")





