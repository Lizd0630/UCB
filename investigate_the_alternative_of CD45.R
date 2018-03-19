
## investigate the alternative of CD45
## CD45: 1:198638671-198757283_+    gene_id "ENSG00000081237"; gene_name "PTPRC"

#### IR ========================================================================================================================================================

load(file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/IR_GTF_merge/IR_known_introns_match.RData")
load(file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/IR_irl_tu.RData")

sum(IR_GTF$chr == 1 & IR_GTF$start >= 198638671 & IR_GTF$end <= 198757283)
# [1] 50
table(IR_GTF[IR_GTF$chr == 1 & IR_GTF$start >= 198638671 & IR_GTF$end <= 198757283, ]$gene_name)
# PTPRC
#    50

IR_GTF[IR_GTF$chr == 1 & IR_GTF$start >= 198638671 & IR_GTF$end <= 198757283, ] -> IR_CD45
length(unique(IR_CD45$event_id))
# [1] 13
length(unique(IR_CD45$transcript_id))
# [1] 13

IR_irl_CD45 <- IR_irl_tu[unique(IR_CD45$event_id),]
rowSums(!is.na(IR_irl_CD45))
#1:198639138-198639225 1:198696910-198699563 1:198699705-198702386
#                    2                    23                     7
#1:198699705-198703297 1:198702531-198703297 1:198706953-198708132
#                    6                    33                    33
#1:198709825-198712952 1:198731727-198732299 1:198734426-198735126
#                   17                     4                     3
#1:198742027-198742231 1:198752372-198752593 1:198752773-198754268
#                   11                     2                     4
#1:198754405-198755905
#                   15

library(reshape2)
summary(melt(IR_irl_CD45)$value)
summary(na.omit(melt(IR_irl_CD45)$value))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.004451 0.083175 0.201347 0.361122 0.600254 1.579303


apply(IR_irl_CD45,1,function(x) which(!is.na(x)))
apply(IR_irl_CD45,1,function(x) colnames(IR_irl_CD45)[which(!is.na(x))])
table(table(as.character(unlist(apply(IR_irl_CD45,1,function(x) colnames(IR_irl_CD45)[which(!is.na(x))])))))
#   1   2   3
# 133   9   3
## Only 3 cells have 3 IR events, and 133 cells only have one IR events.

apply(IR_irl_CD45,1,function(x) x[!is.na(x)])



CD45_IR <- data.frame(event_id = rep(names(rowSums(!is.na(IR_irl_CD45))), as.numeric(rowSums(!is.na(IR_irl_CD45)))),
	cell = as.character(unlist(apply(IR_irl_CD45,1,function(x) colnames(IR_irl_CD45)[which(!is.na(x))]))),
	psi = as.character(unlist(apply(IR_irl_CD45,1,function(x) x[!is.na(x)]))))


gtf <- read.table("/mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.87.gtf",sep="\t",header=F, stringsAsFactors = F)
transcript <- subset(gtf,gtf$V3=="transcript")
test <- as.data.frame(do.call(rbind, strsplit(transcript$V9, split = "[; ]")))

valueFind <- function(x, type = "gene_id"){
  for(i in 1:length(x)){
    if(sum(x == type) == 0){
      return(NA)
    }else{
      return(x[which(x == type)[1]+1])
    }
  }
}

transcript_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
           transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
           gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
           transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
           gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
           transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

merge(x = IR_CD45, y = transcript_infor, by = "transcript_id", all)
table(IR_CD45$transcript_biotype)
#      non_stop_decay processed_transcript       protein_coding
#                   7                    3                   37
#     retained_intron
#                   3
table(IR_CD45[!duplicated(IR_CD45$transcript_id), c("transcript_id", "transcript_biotype")]$transcript_biotype)
#      non_stop_decay processed_transcript       protein_coding
#                   1                    2                    8
#     retained_intron
#                   2




#### SE ========================================================================================================================================================

load("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/SE/RData/all_cells_SE_GTF.RData")
load("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/SE/RData/all_cells_SE_type_and_PSI.RData")

sum(do.call(rbind,strsplit(all_cells_SE_type_and_PSI$ASrange, split="[:-]"))[,1] == 1 & 
	do.call(rbind,strsplit(all_cells_SE_type_and_PSI$ASrange, split="[:-]"))[,2] >= 198638671 & 
	do.call(rbind,strsplit(all_cells_SE_type_and_PSI$ASrange, split="[:-]"))[,3] <= 198757283)
# [1] 0

sum(SE_GTF$V1 == 1 & SE_GTF$V4 >= 198638671 & SE_GTF$V5 <= 198757283, na.rm=T)
# [1] 72

SE_GTF_CD45 <- SE_GTF[which(SE_GTF$V1 == 1 & SE_GTF$V4 >= 198638671 & SE_GTF$V5 <= 198757283), ]
dim(SE_GTF_CD45)
# [1] 72 21
table(SE_GTF_CD45$AStype)
#      exon-skipping_exactly     exon-skipping_half-exon
#                         29                           2
#exon-skipping_multipul-exon
#                         41


length(unique(SE_GTF_CD45$loci))
# [1] 16


test <- as.data.frame(do.call(rbind, strsplit(SE_GTF_CD45$V9, split = "[; ]")))

valueFind <- function(x, type = "gene_id"){
  for(i in 1:length(x)){
    if(sum(x == type) == 0){
      return(NA)
    }else{
      return(x[which(x == type)[1]+1])
    }
  }
}

gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
           transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
           gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
           transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
           gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
           transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

dim(gtf_infor)
# [1] 72  6

SE_GTF_CD45 <- cbind(SE_GTF_CD45,gtf_infor)
table(SE_GTF_CD45[!duplicated(SE_GTF_CD45$transcript_id),]$transcript_biotype)
# protein_coding retained_intron
#              5               1

SE_PSI_CD45 <- all_cells_SE_type_and_PSI[all_cells_SE_type_and_PSI$loci %in% SE_GTF_CD45$loci, ]
dim(SE_PSI_CD45)
# [1]   16 3585

row.names(SE_PSI_CD45) <- SE_PSI_CD45$loci
SE_PSI_CD45 <- SE_PSI_CD45[,12:ncol(SE_PSI_CD45)]
rowSums(!is.na(SE_PSI_CD45))
#1:198692374-198696711-198696910-198699563
#                                      600
#1:198699705-198702386-198702531-198703297
#                                      735
#1:198708262-198709686-198709825-198712952
#                                       28
#1:198692374-198699563-198699705-198703297
#                                       19
#1:198692374-198694017-198696910-198699563
#                                       58
#1:198692374-198699563-198702531-198703297
#                                       15
#1:198692374-198694017-198694141-198696711
#                                       14
#1:198692374-198696711-198699705-198703297
#                                        2
#1:198692374-198696711-198698932-198699563
#                                        2
#1:198692374-198696711-198702531-198703297
#                                       12
#1:198748200-198749415-198749550-198750491
#                                        2
#1:198752773-198754268-198754405-198755905
#                                        1
#1:198699705-198702349-198702531-198703297
#                                        4
#1:198706953-198708132-198709825-198712952
#                                        1
#1:198692374-198694017-198694141-198699563
#                                        1
#1:198699705-198701774-198702531-198703297
#                                        1

summary(na.omit(melt(SE_PSI_CD45)$value))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0006127 0.1656313 0.5384609 0.5303473 0.9459642 0.9995992

table(table(as.character(unlist(apply(SE_PSI_CD45,1,function(x) colnames(SE_PSI_CD45)[which(!is.na(x))])))))
#   1   2   3
# 732 365  11
## Only 11 cells have 3 IR events, and 732 cells only have one IR events.


(SE_PSI_CD45["1:198692374-198696711-198696910-198699563",!is.na()])
intersect(colnames(SE_PSI_CD45)[!is.na(SE_PSI_CD45["1:198692374-198696711-198696910-198699563",])],colnames(SE_PSI_CD45)[!is.na(SE_PSI_CD45["1:198699705-198702386-198702531-198703297",])])
length(intersect(colnames(SE_PSI_CD45)[!is.na(SE_PSI_CD45["1:198692374-198696711-198696910-198699563",])],colnames(SE_PSI_CD45)[!is.na(SE_PSI_CD45["1:198699705-198702386-198702531-198703297",])]))
# [1] 324
## 600 intersect 735 = 324

SE_GTF_CD45[SE_GTF_CD45$loci == "1:198692374-198696711-198696910-198699563", ]$AStype
#  "exon-skipping_exactly"
SE_GTF_CD45[SE_GTF_CD45$loci == "1:198699705-198702386-198702531-198703297", ]$AStype
#  "exon-skipping_exactly"






#### A3SS ========================================================================================================================================================




load("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/A3SS/RData/Parsed_A3SS_PSI_GTF.RData")

A3SS_GTF <- na.omit(A3SS_GTF)
sum(A3SS_GTF$chr == 1 & A3SS_GTF$start >= 198638671 & A3SS_GTF$end <= 198757283, na.rm=T)
# [1] 0





#### A5SS ========================================================================================================================================================




load("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/A5SS/RData/Parsed_A5SS_PSI_GTF.RData")

A5SS_GTF <- na.omit(A5SS_GTF)
sum(A5SS_GTF$chr == 1 & A5SS_GTF$start >= 198638671 & A5SS_GTF$end <= 198757283, na.rm=T)
# [1] 0









#### MXE ========================================================================================================================================================




load("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/MXE/RData/all_cells_MXE_type_and_existing.RData")

A5SS_GTF <- na.omit(A5SS_GTF)
sum(all_cells_MXE_type_and_existing$chr == 1 & all_cells_MXE_type_and_existing$start >= 198638671 & all_cells_MXE_type_and_existing$end <= 198757283, na.rm=T)
# [1] 2
## There are two MXE in CD45

MXE_CD45 <- all_cells_MXE_type_and_existing[all_cells_MXE_type_and_existing$chr == 1 & all_cells_MXE_type_and_existing$start >= 198638671 & all_cells_MXE_type_and_existing$end <= 198757283,]
dim(MXE_CD45)
# [1]    2 3590

rowSums(MXE_CD45[,17:ncol(MXE_CD45)] == 1)
# 1206 1701
#    2    1
## It makes no sense








#### plot in tsne ======================================================================================================================




load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")
tsne <- as.data.frame(Combine@dr$tsne@cell.embeddings)

dim(IR_irl_CD45)
dim(SE_PSI_CD45)



identical(colnames(SE_PSI_CD45), colnames(IR_irl_CD45))
rbind(SE_PSI_CD45,IR_irl_CD45) -> PSI
PSI <- PSI[,row.names(tsne)]
dim(PSI)

data.frame(rowSums(!is.na(PSI)))
data.frame(rowSums(!is.na(PSI)), row.names = 1:nrow(PSI), name = row.names(PSI))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/top_AS_in_cells.pdf")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

plot(x = tsne$tSNE_1, y = tsne$tSNE_2, pch = 19, cex = .6, 
     col = rainbow(1, start =.1, alpha = .1), xlab = "tSNE_1", ylab = "tSNE_2")

points(x = tsne[colnames(PSI)[!is.na(PSI[1,])],]$tSNE_1, 
       y = tsne[colnames(PSI)[!is.na(PSI[1,])],]$tSNE_2, 
       col = rainbow(1), pch = 19, cex = .6)
points(x = tsne[colnames(PSI)[!is.na(PSI[2,])],]$tSNE_1, 
       y = tsne[colnames(PSI)[!is.na(PSI[2,])],]$tSNE_2, 
       col = rainbow(1,start = .2), pch = 19, cex = .6)
# Add legend to top right, outside plot region
legend("topright", inset=c(-0.2,.4), legend=c("SE1","SE2"), pch=19, cex = 1, bty = "n",
       col = c(rainbow(1),rainbow(1,start = .2)))
dev.off()



read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", sep = "\t", header = F, row.names = 1) -> celltype
colnames(celltype) <- "Cell type"
PSI[1:2,] -> psitu
psitu[is.na(psitu)] <- 0

heatmap(as.matrix(psitu), ColSideColors = as.character(celltype$V2), labRow = "")

identical(row.names(celltype), colnames(psitu))

library(pheatmap)

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/top_AS_in_cells_heatmap.pdf")

pheatmap(mat = as.matrix(psitu), annotation_col = celltype, show_colnames = F, show_rownames = F)

test <- psitu[,colSums(psitu) > 0]

pheatmap(mat = as.matrix(test), show_colnames = F, show_rownames = F,
         annotation_col = subset.data.frame(celltype, row.names(celltype) %in% colnames(test)))
dev.off()


celltest <- subset.data.frame(celltype, row.names(celltype) %in% colnames(test))
celltest$YN <- 1
celltest <- celltest[order(celltest$`Cell type`),]

celltest -> test2
test2$cell <- row.names(test2)
as.data.frame(t(test)) -> test1
test1$cell <- row.names(test1)

test3 <- merge(test1, test2, by = "cell")
head(test3)
colnames(test3) <- c("cell","SE1","SE2","Cell","YN")

test3->test4

test4[test4 == 0] <- NA


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/top_AS_in_cells_boxplot.pdf")

par(mar=c(8.1, 4.1, 4.1, 8.1), xpd=TRUE)

boxplot(SE1~Cell, test4, boxwex = 0.25, at = 1:6 - 0.2,col = "yellow", xaxt='n')
boxplot(SE2~Cell, test4, boxwex = 0.25, at = 1:6 + 0.2,col = "orange", add = T, xaxt='n')
axis(side = 1, at = 1:6,line = .3, tick = F, cex.axis=1, las = 2, srt = 45,
     labels = paste(names(table(celltype$`Cell type`)),"\n", table(celltype$`Cell type`)))

table(test4[!is.na(test4$SE1),]$`Cell`)
table(test4[!is.na(test4$SE2),]$`Cell`)

axis(side = 3, at = 1:6 - .2, tick = F, labels = table(test4[!is.na(test4$SE1),]$`Cell`), line = -1)
axis(side = 3, at = 1:6 + .2, tick = F, labels = table(test4[!is.na(test4$SE2),]$`Cell`), line = -1)

legend("topright", inset=c(-0.2,.4), legend=c("SE1","SE2"), pch=19, cex = 1, bty = "n",
       col = c("yellow","orange"))

dev.off()







#### intron centric PSI for CD45 =================================================================================================================





load(file = "/mnt/data5/BGI/UCB/tangchao/DSU/RData/psi_list_left_right_table.RData")

data.frame(do.call(rbind,strsplit(row.names(psi_sj_same_start_table), split="[:-]")), stringsAsFactors = F) -> sjinfo
sum(sjinfo$X1 == 1 & sjinfo$X2 >= 198638671 & sjinfo$X3 <= 198757283, na.rm=T)
# [1] 204

data.frame(do.call(rbind,strsplit(row.names(psi_sj_same_end_table), split="[:-]")), stringsAsFactors = F) -> sjinfo2
sum(sjinfo2$X1 == 1 & sjinfo2$X2 >= 198638671 & sjinfo2$X3 <= 198757283, na.rm=T)
# [1] 76


sj_CD45 <- rbind(psi_sj_same_start_table[sjinfo$X1 == 1 & sjinfo$X2 >= 198638671 & sjinfo$X3 <= 198757283,],
				psi_sj_same_end_table[sjinfo2$X1 == 1 & sjinfo2$X2 >= 198638671 & sjinfo2$X3 <= 198757283,])

sj_CD45[is.na(sj_CD45)] <- 0

summary(colSums(sj_CD45))
sum(colSums(sj_CD45) == 0)
# [1] 502

sj_CD45 <- sj_CD45[,row.names(celltype)] 


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/All_SJ_in_cells_heatmap.pdf")

pheatmap(mat = as.matrix(sj_CD45), annotation_col = celltype, show_colnames = F, show_rownames = F)

pheatmap(mat = as.matrix(sj_CD45[apply(sj_CD45,1,sd)>=.1,]), annotation_col = celltype, show_colnames = F, show_rownames = F)

dev.off()



identical(colnames(sj_CD45), row.names(celltype))
# [1] TRUE

test <- data.frame(cell = row.names(celltype), CellType = celltype$`Cell type`,
					SJ = sj_CD45[1,])

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/All_SJ_in_cells_boxplot.pdf")

par(mar=c(8.1, 4.1, 4.1, 8.1), xpd=TRUE)

boxplot(SJ~CellType, test, boxwex = 0.25, col = "yellow", xaxt='n', main = row.names(sj_CD45)[1])
axis(side = 1, at = 1:6,line = .3, tick = F, cex.axis=1, las = 2, 
     labels = paste(names(table(celltype$`Cell type`)),"\n", table(celltype$`Cell type`)))

axis(side = 3, at = 1:6, tick = F, labels = table(test[test$SJ>0,]$`CellType`), line = -1)

dev.off()


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/All_SJ_in_cells_boxplot.pdf")
par(mar=c(8.1, 4.1, 4.1, 8.1), xpd=TRUE)
for(i in 1:nrow(sj_CD45)){
	print(i)
	test <- data.frame(cell = row.names(celltype), CellType = celltype$`Cell type`,SJ = sj_CD45[i,])
	test[test==0] <- NA

	if(sum(test$SJ, na.rm = T) > 0){
		boxplot(SJ~CellType, test, boxwex = 0.25, col = "yellow", xaxt='n', main = row.names(sj_CD45)[i])
	axis(side = 1, at = 1:6,line = .3, tick = F, cex.axis=1, las = 2, 
     labels = paste(names(table(celltype$`Cell type`)),"\n", table(celltype$`Cell type`)))

	axis(side = 3, at = 1:6, tick = F, labels = table(test[test$SJ>0,]$`CellType`), line = -1)

	}

	
}

dev.off()




pdf("/mnt/data5/BGI/UCB/tangchao/CD45/Meaningful_SJ_in_cells_boxplot.pdf")
par(mar=c(8.1, 4.1, 4.1, 8.1), xpd=TRUE)
for(i in 1:nrow(sj_CD45)){
	print(i)
	test <- data.frame(cell = row.names(celltype), CellType = celltype$`Cell type`,SJ = sj_CD45[i,])
	test[test==0] <- NA

	if(sum(test$SJ, na.rm = T) > 0 & sd(test$SJ, na.rm = T) > 0 & sum(!is.na(test$SJ)) > 500){
		boxplot(SJ~CellType, test, boxwex = 0.25, col = "yellow", xaxt='n', main = row.names(sj_CD45)[i])
	axis(side = 1, at = 1:6,line = .3, tick = F, cex.axis=1, las = 2, 
     labels = paste(names(table(celltype$`Cell type`)),"\n", table(celltype$`Cell type`)))

	axis(side = 3, at = 1:6, tick = F, labels = table(test[test$SJ>0,]$`CellType`), line = -1)

	}

	
}

dev.off()




pdf("/mnt/data5/BGI/UCB/tangchao/CD45/Meaningful2_SJ_in_cells_boxplot.pdf")
par(mar=c(8.1, 4.1, 4.1, 8.1), xpd=TRUE)
for(i in 1:nrow(sj_CD45)){
	print(i)
	test <- data.frame(cell = row.names(celltype), CellType = celltype$`Cell type`,SJ = sj_CD45[i,])
	test[test==0] <- NA

	if(sum(test$SJ, na.rm = T) > 0 & sd(test$SJ, na.rm = T) > 0 & sum(!is.na(test$SJ)) > 500 & median(test$SJ, na.rm = T) != 1){
		boxplot(SJ~CellType, test, boxwex = 0.25, col = "yellow", xaxt='n', main = row.names(sj_CD45)[i])
	axis(side = 1, at = 1:6,line = .3, tick = F, cex.axis=1, las = 2, 
     labels = paste(names(table(celltype$`Cell type`)),"\n", table(celltype$`Cell type`)))

	axis(side = 3, at = 1:6, tick = F, labels = table(test[test$SJ>0,]$`CellType`), line = -1)

	}

	
}

dev.off()

