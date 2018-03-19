
#### Number of every AS type ======================================================================================================================
## SE
setwd("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/")
load("./SE/RData/all_cells_SE_type_and_PSI.RData")
table(all_cells_SE_type_and_PSI$AStype)
#      exon-skipping_exactly     exon-skipping_half-exon
#                       8311                         520
#exon-skipping_multipul-exon    intergenic-splicing_site
#                       3943                        2747
#                      Other
#                        324

rowSums(!is.na(all_cells_SE_type_and_PSI[,12:ncol(all_cells_SE_type_and_PSI)])) -> se_psi_rowsum
table(all_cells_SE_type_and_PSI[se_psi_rowsum>1,]$AStype)
#      exon-skipping_exactly     exon-skipping_half-exon
#                       4590                         146
#exon-skipping_multipul-exon    intergenic-splicing_site
#                       1469                         802
#                      Other
#                         89


## A3SS
load("./A3SS/RData/Parsed_A3SS_PSI_GTF.RData")
dim(A3SS_type_and_PSI)
# [1] 5874 3584
rowSums(!is.na(A3SS_type_and_PSI[,11:ncol(A3SS_type_and_PSI)])) -> a3ss_psi_rowsum
sum(a3ss_psi_rowsum>1)
# [1] 3179


## A5SS
load("./A5SS/RData/Parsed_A5SS_PSI_GTF.RData")
dim(A5SS_type_and_PSI)
# [1] 3908 3584
rowSums(!is.na(A5SS_type_and_PSI[,11:ncol(A5SS_type_and_PSI)])) -> a5ss_psi_rowsum
sum(a5ss_psi_rowsum>1)
# [1] 2256


## AFE
load("./AFE/RData/all_cells_AFE_PSI_GTF.RData")
dim(all_cells_AFE_PSI)
# [1] 2289 3584
rowSums(!is.na(all_cells_AFE_PSI[,11:ncol(all_cells_AFE_PSI)])) -> afe_psi_rowsum
sum(afe_psi_rowsum>1)
# [1] 1361


## ALE
load("./ALE/RData/all_cells_ALE_PSI_GTF.RData")
dim(all_cells_ALE_PSI)
# [1] 1211 3584
rowSums(!is.na(all_cells_ALE_PSI[,11:ncol(all_cells_ALE_PSI)])) -> ale_psi_rowsum
sum(ale_psi_rowsum>1)
# [1] 624


load("./MXE/RData/all_cells_MXE_type_and_existing.RData")
dim(all_cells_MXE_type_and_existing)
# [1]  376 3590
rowSums(all_cells_MXE_type_and_existing[,17:ncol(all_cells_MXE_type_and_existing)]) -> mxe_psi_rowsum
sum(mxe_psi_rowsum>1)
#[1] 154


#### Gene function of every AS type ======================================================================================================================

## SE -----------------------------------------------------------------------------------------------------------------------------------------------------

all_cells_SE_type_and_PSI[se_psi_rowsum > 1 & all_cells_SE_type_and_PSI$AStype == "exon-skipping_exactly",]$loci -> se_loci_more1
all_cells_SE_type_and_PSI[se_psi_rowsum == 1 & all_cells_SE_type_and_PSI$AStype == "exon-skipping_exactly",]$loci -> se_loci_equal1
load("./SE/RData/all_cells_SE_GTF.RData")

length(SE_GTF[SE_GTF$loci %in% se_loci_more1, ]$V9)

test <- as.data.frame(do.call(rbind, strsplit(SE_GTF[SE_GTF$loci %in% se_loci_more1, ]$V9, split = "[; ]")))


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
# [1] 25575    6
gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 7877    6
length(unique(gtf_infor$gene_id))
# [1] 3004
length(unique(gtf_infor$transcript_id))
# [1] 7877

table(gtf_infor$gene_biotype)
#                         antisense                          IG_C_gene
#                                27                                  1
#                           lincRNA             polymorphic_pseudogene
#                                16                                  1
#              processed_transcript                     protein_coding
#                                23                               3949
#                 sense_overlapping   transcribed_processed_pseudogene
#                                 1                                  2
#transcribed_unprocessed_pseudogene                 unitary_pseudogene
#                                27                                  1
#            unprocessed_pseudogene
#                                 4
#                         antisense                          IG_C_gene
#                                34                                  2
#                           lincRNA                           misc_RNA
#                                49                                  3
#              processed_transcript                     protein_coding
#                                85                               7647
#  transcribed_processed_pseudogene     transcribed_unitary_pseudogene
#                                 1                                  5
#transcribed_unprocessed_pseudogene                          TR_C_gene
#                                35                                  4
#                unitary_pseudogene             unprocessed_pseudogene
#                                 9                                  3
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype)
#                         antisense                          IG_C_gene
#                                24                                  2
#                           lincRNA                           misc_RNA
#                                24                                  3
#              processed_transcript                     protein_coding
#                                30                               2891
#  transcribed_processed_pseudogene     transcribed_unitary_pseudogene
#                                 1                                  1
#transcribed_unprocessed_pseudogene                          TR_C_gene
#                                19                                  3
#                unitary_pseudogene             unprocessed_pseudogene
#                                 3                                  3

gtf_infor -> SE_gene_info_more1
unique(as.character(gtf_infor$gene_id)) -> SE_gene_list_more1

unique(gtf_infor[gtf_infor$gene_biotype == "protein_coding", ]$gene_id) -> SE_pcgene_list_more1

## GO enrichment


## equal
length(se_loci_equal1)
# [1] 3721
test <- as.data.frame(do.call(rbind, strsplit(SE_GTF[SE_GTF$loci %in% se_loci_equal1, ]$V9, split = "[; ]")))


gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
           transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
           gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
           transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
           gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
           transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

dim(gtf_infor)
# [1] 20225    6
gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 6591    6
length(unique(gtf_infor$gene_id))
# [1] 2864
length(unique(gtf_infor$transcript_id))
# [1] 6591

table(gtf_infor$gene_biotype)
#                         antisense                          IG_C_gene
#                                45                                  1
#                         IG_J_gene                            lincRNA
#                                 1                                 55
#                          misc_RNA               processed_transcript
#                                 1                                107
#                    protein_coding                             snoRNA
#                              6310                                  1
#  transcribed_processed_pseudogene     transcribed_unitary_pseudogene
#                                 2                                  5
#transcribed_unprocessed_pseudogene                          TR_C_gene
#                                51                                  5
#                unitary_pseudogene             unprocessed_pseudogene
#                                 4                                  3
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype)
#                         antisense                          IG_C_gene
#                                28                                  1
#                         IG_J_gene                            lincRNA
#                                 1                                 29
#                          misc_RNA               processed_transcript
#                                 1                                 33
#                    protein_coding                             snoRNA
#                              2731                                  1
#  transcribed_processed_pseudogene     transcribed_unitary_pseudogene
#                                 2                                  1
#transcribed_unprocessed_pseudogene                          TR_C_gene
#                                26                                  4
#                unitary_pseudogene             unprocessed_pseudogene
#                                 3                                  3

gtf_infor -> SE_gene_info_equal1
unique(as.character(gtf_infor$gene_id)) -> SE_gene_list_equal1

unique(gtf_infor[gtf_infor$gene_biotype == "protein_coding", ]$gene_id) -> SE_pcgene_list_equal1



## intersect of more 1 and equal 1
length(SE_gene_list_more1)
# [1] 3004
length(SE_gene_list_equal1)
# [1] 2864
length(intersect(SE_gene_list_more1,SE_gene_list_equal1))
# [1] 1259

intersect(SE_gene_list_more1,SE_gene_list_equal1) -> SE_gene_list_com

SE_gene_list_equal1[!SE_gene_list_equal1 %in% SE_gene_list_more1] -> SE_gene_list_only_1cell

SE_gene_list_more1[!SE_gene_list_more1 %in% SE_gene_list_equal1] -> SE_gene_list_exactly_more_1cell
length(SE_gene_list_exactly_more_1cell)
# [1] 1745

#### GO
library(org.Hs.eg.db)#####datdbase target to human-geneIDs
library(AnnotationDbi)
biocLite("clusterProfiler")
library(DOSE)
library(GO.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)

for(i in c("SE_gene_list_more1","SE_pcgene_list_more1","SE_gene_list_equal1","SE_pcgene_list_equal1","SE_gene_list_com","SE_gene_list_only_1cell","SE_gene_list_exactly_more_1cell")){
  ## enrichment
  eval(parse(text = paste("genelist <- ", i, sep ="")))
  ego_BP <- enrichGO(gene = genelist,
                     OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                     ont  = "BP", pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
  ego_MF <- enrichGO(gene = genelist,
                     OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                     ont  = "MF", pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
  ego_CC <- enrichGO(gene = genelist,
                     OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                     ont  = "CC", pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
  
  pdf(paste("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/SE/", i, "_BP.pdf", sep = ""),height = 10,width = 15)
  print(dotplot(ego_BP,showCategory = 10,font.size = 16,title = "BP"))
  plotGOgraph(ego_BP,firstSigNodes = 15, useInfo = "all")
  enrichMap(ego_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
  dev.off()
  pdf(paste("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/SE/", i, "_MF.pdf", sep = ""),height = 10,width = 15)
  print(dotplot(ego_MF,showCategory = 10,font.size = 16,title = "MF"))
  plotGOgraph(ego_MF,firstSigNodes = 15, useInfo = "all")
  enrichMap(ego_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
  dev.off()
  pdf(paste("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/SE/", i, "_CC.pdf", sep = ""),height = 10,width = 15)
  print(dotplot(ego_CC,showCategory = 10,font.size = 16,title = "CC"))
  plotGOgraph(ego_CC,firstSigNodes = 15, useInfo = "all")
  enrichMap(ego_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = .4)
  dev.off()
}







test <- all_cells_SE_type_and_PSI[se_psi_rowsum == 1 & all_cells_SE_type_and_PSI$AStype == "exon-skipping_exactly",]
dim(test)
# [1] 3721 3585
test2 <- colSums(!is.na(test[,12:ncol(test)]))
summary(test2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.000   0.000   1.000   1.041   1.000  53.000
table(test2)
#   0    1    2    3    4    5    6    7    8    9   10   11   12   16   20   32
#1638 1053  483  213   92   47   13   11    9    5    1    3    2    1    1    1
#  53
#   1
pdf("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/SE/Histogram_of_cells_has_specific_SE.pdf")
hist(test2, freq=F, breaks=20, col = rainbow(1,start = .2, end = .3,alpha = .4), 
	main = "Histogram of cells has specific SE(only in one cell)", xlab = "No. specific SE")
dev.off()



## A3SS -----------------------------------------------------------------------------------------------------------------------------------------------------

load("./A3SS/RData/Parsed_A3SS_PSI_GTF.RData")
dim(A3SS_type_and_PSI)
# [1] 5874 3584
rowSums(!is.na(A3SS_type_and_PSI[,11:ncol(A3SS_type_and_PSI)])) -> a3ss_psi_rowsum
sum(a3ss_psi_rowsum>1)
# [1] 3179

length(a3ss_psi_rowsum)
sum(a3ss_psi_rowsum>1)
A3SS_type_and_PSI$name[a3ss_psi_rowsum>1] -> a3ss_loci_more1
A3SS_type_and_PSI$name[a3ss_psi_rowsum==1] -> a3ss_loci_equal1


test <- as.data.frame(do.call(rbind, strsplit(A3SS_GTF[!is.na(A3SS_GTF$V9),]$V9, split = "[; ]")))

gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
           transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
           gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
           transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
           gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
           transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 8073    6
length(unique(gtf_infor$gene_id))
# [1] 3603
length(unique(gtf_infor$transcript_id))
# [1] 8073

table(gtf_infor$gene_biotype)
#                         antisense                          IG_C_gene
#                                92                                  4
#                         IG_V_gene                    IG_V_pseudogene
#                                10                                  1
#                           lincRNA             polymorphic_pseudogene
#                                91                                  3
#              processed_transcript                     protein_coding
#                                65                               7719
#  transcribed_processed_pseudogene     transcribed_unitary_pseudogene
#                                 1                                  4
#transcribed_unprocessed_pseudogene                          TR_C_gene
#                                48                                  7
#                         TR_V_gene                    TR_V_pseudogene
#                                12                                  1
#                unitary_pseudogene             unprocessed_pseudogene
#                                 7                                  8
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype)
#                         antisense                          IG_C_gene
#                                60                                  3
#                         IG_V_gene                    IG_V_pseudogene
#                                10                                  1
#                           lincRNA             polymorphic_pseudogene
#                                55                                  1
#              processed_transcript                     protein_coding
#                                31                               3376
#  transcribed_processed_pseudogene     transcribed_unitary_pseudogene
#                                 1                                  2
#transcribed_unprocessed_pseudogene                          TR_C_gene
#                                33                                  5
#                         TR_V_gene                    TR_V_pseudogene
#                                12                                  1
#                unitary_pseudogene             unprocessed_pseudogene
#                                 4                                  8


test <- as.data.frame(do.call(rbind, strsplit(na.omit(A3SS_GTF[A3SS_GTF$name %in% a3ss_loci_more1,]$V9), split = "[; ]")))

gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
           transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
           gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
           transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
           gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
           transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 4812    6
length(unique(gtf_infor$gene_id))
# [1] 2218
length(unique(gtf_infor$transcript_id))
# [1] 4812

table(gtf_infor$gene_biotype)
#                         antisense                          IG_C_gene
#                                45                                  4
#                         IG_V_gene                    IG_V_pseudogene
#                                 2                                  1
#                           lincRNA               processed_transcript
#                                35                                 43
#                    protein_coding transcribed_unprocessed_pseudogene
#                              4643                                 19
#                         TR_C_gene                          TR_V_gene
#                                 5                                  6
#                   TR_V_pseudogene                 unitary_pseudogene
#                                 1                                  5
#            unprocessed_pseudogene
#                                 3
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype)
#                         antisense                          IG_C_gene
#                                26                                  3
#                         IG_V_gene                    IG_V_pseudogene
#                                 2                                  1
#                           lincRNA               processed_transcript
#                                25                                 22
#                    protein_coding transcribed_unprocessed_pseudogene
#                              2110                                 14
#                         TR_C_gene                          TR_V_gene
#                                 3                                  6
#                   TR_V_pseudogene                 unitary_pseudogene
#                                 1                                  2
#            unprocessed_pseudogene
#                                 3
unique(gtf_infor$gene_id) -> A3SS_gene_list_more1
unique(gtf_infor[gtf_infor$gene_biotype == "protein_coding",]$gene_id) -> A3SS_pcgene_list_more1


test <- as.data.frame(do.call(rbind, strsplit(na.omit(A3SS_GTF[A3SS_GTF$name %in% a3ss_loci_equal1,]$V9), split = "[; ]")))

gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
           transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
           gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
           transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
           gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
           transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 4044    6
length(unique(gtf_infor$gene_id))
# [1] 2183
length(unique(gtf_infor$transcript_id))
# [1] 4044

table(gtf_infor$gene_biotype)
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype)
#                         antisense                          IG_C_gene
#                                38                                  2
#                         IG_V_gene                            lincRNA
#                                 8                                 34
#            polymorphic_pseudogene               processed_transcript
#                                 1                                 18
#                    protein_coding   transcribed_processed_pseudogene
#                              2040                                  1
#    transcribed_unitary_pseudogene transcribed_unprocessed_pseudogene
#                                 2                                 21
#                         TR_C_gene                          TR_V_gene
#                                 5                                  6
#                unitary_pseudogene             unprocessed_pseudogene
#                                 2                                  5
unique(gtf_infor$gene_id) -> A3SS_gene_list_equal1
unique(gtf_infor[gtf_infor$gene_biotype == "protein_coding",]$gene_id) -> A3SS_pcgene_list_equal1

intersect(A3SS_gene_list_more1, A3SS_gene_list_equal1) -> A3SS_gene_list_com
length(A3SS_gene_list_com)
# [1] 798
A3SS_gene_list_equal1[!A3SS_gene_list_equal1 %in% A3SS_gene_list_more1] -> A3SS_gene_list_only_1cell

A3SS_gene_list_more1[!A3SS_gene_list_more1 %in% A3SS_gene_list_equal1] -> A3SS_gene_list_exactly_more_1cell
length(A3SS_gene_list_only_1cell)
# [1] 1385
length(A3SS_gene_list_exactly_more_1cell)
# [1] 1420



#A3SS_gene_list_more1,A3SS_pcgene_list_more1,A3SS_gene_list_equal1,A3SS_pcgene_list_equal1,A3SS_gene_list_com,A3SS_gene_list_only_1cell,A3SS_gene_list_exactly_more_1cell

## GO enrichment

#######
## all gene

for(i in c("A3SS_gene_list_more1","A3SS_pcgene_list_more1","A3SS_gene_list_equal1","A3SS_pcgene_list_equal1","A3SS_gene_list_com","A3SS_gene_list_only_1cell","A3SS_gene_list_exactly_more_1cell")){
  ## enrichment
  eval(parse(text = paste("genelist <- ", i, sep ="")))
  ego_BP <- enrichGO(gene = genelist,
                     OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                     ont  = "BP", pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
  ego_MF <- enrichGO(gene = genelist,
                     OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                     ont  = "MF", pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
  ego_CC <- enrichGO(gene = genelist,
                     OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                     ont  = "CC", pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
  
  pdf(paste("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/A3SS/", i, "_BP.pdf", sep = ""),height = 10,width = 15)
  print(dotplot(ego_BP,showCategory = 10,font.size = 16,title = "BP"))
  suppressMessages(plotGOgraph(ego_BP,firstSigNodes = 15, useInfo = "all"))
  enrichMap(ego_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
  dev.off()
  pdf(paste("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/A3SS/", i, "_MF.pdf", sep = ""),height = 10,width = 15)
  print(dotplot(ego_MF,showCategory = 10,font.size = 16,title = "MF"))
  suppressMessages(plotGOgraph(ego_MF,firstSigNodes = 15, useInfo = "all"))
  enrichMap(ego_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
  dev.off()
  pdf(paste("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/A3SS/", i, "_CC.pdf", sep = ""),height = 10,width = 15)
  print(dotplot(ego_CC,showCategory = 10,font.size = 16,title = "CC"))
  suppressMessages(plotGOgraph(ego_CC,firstSigNodes = 15, useInfo = "all"))
  enrichMap(ego_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = .4)
  dev.off()
}











## A5SS -----------------------------------------------------------------------------------------------------------------------------------------------------

load("./A5SS/RData/Parsed_A5SS_PSI_GTF.RData")
dim(A5SS_type_and_PSI)
# [1] 3908 3584
rowSums(!is.na(A5SS_type_and_PSI[,11:ncol(A5SS_type_and_PSI)])) -> a5ss_psi_rowsum
sum(a5ss_psi_rowsum>1)
# [1] 2256

A5SS_type_and_PSI$name[a5ss_psi_rowsum>1] -> a5ss_loci_more1
A5SS_type_and_PSI$name[a5ss_psi_rowsum==1] -> a5ss_loci_equal1
length(a5ss_loci_more1)
# [1] 2256
length(a5ss_loci_equal1)
# [1] 1652


test <- as.data.frame(do.call(rbind, strsplit(na.omit(A5SS_GTF[A5SS_GTF$name %in% a5ss_loci_more1,]$V9), split = "[; ]")))

gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
           transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
           gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
           transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
           gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
           transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 3476    6
length(unique(gtf_infor$gene_id))
# [1] 1735
length(unique(gtf_infor$transcript_id))
# [1] 3476

table(gtf_infor$gene_biotype)
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype)
#                         antisense                          IG_C_gene
#                                33                                  1
#                           lincRNA               processed_transcript
#                                26                                 22
#                    protein_coding     transcribed_unitary_pseudogene
#                              1623                                  1
#transcribed_unprocessed_pseudogene                          TR_C_gene
#                                17                                  2
#                         TR_J_gene                          TR_V_gene
#                                 1                                  6
#            unprocessed_pseudogene
#                                 3

unique(gtf_infor$gene_id) -> A5SS_gene_list_more1
unique(gtf_infor[gtf_infor$gene_biotype == "protein_coding",]$gene_id) -> A5SS_pcgene_list_more1


test <- as.data.frame(do.call(rbind, strsplit(na.omit(A5SS_GTF[A5SS_GTF$name %in% a5ss_loci_equal1,]$V9), split = "[; ]")))

gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
           transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
           gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
           transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
           gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
           transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 2588    6
length(unique(gtf_infor$gene_id))
# [1] 1476
length(unique(gtf_infor$transcript_id))
# [1] 2588

table(gtf_infor$gene_biotype)
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype)
#                         antisense                          IG_C_gene
#                                29                                  2
#                           lincRNA               processed_transcript
#                                29                                 24
#                    protein_coding                  sense_overlapping
#                              1360                                  1
#  transcribed_processed_pseudogene     transcribed_unitary_pseudogene
#                                 1                                  3
#transcribed_unprocessed_pseudogene                          TR_C_gene
#                                17                                  2
#                         TR_J_gene                          TR_V_gene
#                                 1                                  4
#                unitary_pseudogene             unprocessed_pseudogene
#                                 1                                  2

unique(gtf_infor$gene_id) -> A5SS_gene_list_equal1
unique(gtf_infor[gtf_infor$gene_biotype == "protein_coding",]$gene_id) -> A5SS_pcgene_list_equal1

intersect(A5SS_gene_list_more1, A5SS_gene_list_equal1) -> A5SS_gene_list_com
length(A5SS_gene_list_com)
# [1] 412
A5SS_gene_list_equal1[!A5SS_gene_list_equal1 %in% A5SS_gene_list_more1] -> A5SS_gene_list_only_1cell

A5SS_gene_list_more1[!A5SS_gene_list_more1 %in% A5SS_gene_list_equal1] -> A5SS_gene_list_exactly_more_1cell
length(A5SS_gene_list_only_1cell)
# [1] 1064
length(A5SS_gene_list_exactly_more_1cell)
# [1] 1323


for(i in c("A5SS_gene_list_more1","A5SS_pcgene_list_more1","A5SS_gene_list_equal1","A5SS_pcgene_list_equal1","A5SS_gene_list_com","A5SS_gene_list_only_1cell","A5SS_gene_list_exactly_more_1cell")){
  ## enrichment
  eval(parse(text = paste("genelist <- ", i, sep ="")))
  ego_BP <- enrichGO(gene = genelist,
                     OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                     ont  = "BP", pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
  ego_MF <- enrichGO(gene = genelist,
                     OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                     ont  = "MF", pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
  ego_CC <- enrichGO(gene = genelist,
                     OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                     ont  = "CC", pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
  
  pdf(paste("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/A5SS/", i, "_BP.pdf", sep = ""),height = 10,width = 15)
  print(dotplot(ego_BP,showCategory = 10,font.size = 16,title = "BP"))
  plotGOgraph(ego_BP,firstSigNodes = 15, useInfo = "all")
  enrichMap(ego_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
  dev.off()
  pdf(paste("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/A5SS/", i, "_MF.pdf", sep = ""),height = 10,width = 15)
  print(dotplot(ego_MF,showCategory = 10,font.size = 16,title = "MF"))
  plotGOgraph(ego_MF,firstSigNodes = 15, useInfo = "all")
  enrichMap(ego_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
  dev.off()
  pdf(paste("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/A5SS/", i, "_CC.pdf", sep = ""),height = 10,width = 15)
  print(dotplot(ego_CC,showCategory = 10,font.size = 16,title = "CC"))
  plotGOgraph(ego_CC,firstSigNodes = 15, useInfo = "all")
  enrichMap(ego_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = .4)
  dev.off()
}





## AFE -----------------------------------------------------------------------------------------------------------------------------------------------------


## AFE
load("./AFE/RData/all_cells_AFE_PSI_GTF.RData")
dim(all_cells_AFE_PSI)
# [1] 2289 3584
rowSums(!is.na(all_cells_AFE_PSI[,11:ncol(all_cells_AFE_PSI)])) -> afe_psi_rowsum
sum(afe_psi_rowsum>1)
# [1] 1361

all_cells_AFE_PSI$name[afe_psi_rowsum>1] -> afe_loci_more1
all_cells_AFE_PSI$name[afe_psi_rowsum==1] -> afe_loci_equal1
length(afe_loci_more1)
# [1] 1361
length(afe_loci_equal1)
# [1] 928


test <- as.data.frame(do.call(rbind, strsplit(na.omit(AFE_GTF[AFE_GTF$name %in% afe_loci_more1,]$V9), split = "[; ]")))

gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
           transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
           gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
           transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
           gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
           transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 1649    6
length(unique(gtf_infor$gene_id))
# [1] 1204
length(unique(gtf_infor$transcript_id))
# [1] 1649

table(gtf_infor$gene_biotype)
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype)
#                         antisense                          IG_J_gene
#                                25                                  6
#                         IG_V_gene                            lincRNA
#                                 1                                 24
#              processed_transcript                     protein_coding
#                                15                               1102
#                            snoRNA   transcribed_processed_pseudogene
#                                 3                                  1
#transcribed_unprocessed_pseudogene                          TR_J_gene
#                                 7                                 16
#                         TR_V_gene
#                                 4
unique(gtf_infor$gene_id) -> AFE_gene_list_more1
unique(gtf_infor[gtf_infor$gene_biotype == "protein_coding",]$gene_id) -> AFE_pcgene_list_more1


test <- as.data.frame(do.call(rbind, strsplit(na.omit(AFE_GTF[AFE_GTF$name %in% afe_loci_equal1,]$V9), split = "[; ]")))

gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
           transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
           gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
           transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
           gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
           transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 1128    6
length(unique(gtf_infor$gene_id))
# [1] 871
length(unique(gtf_infor$transcript_id))
# [1] 1128

table(gtf_infor$gene_biotype)
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype)
#                         antisense                          IG_J_gene
#                                35                                  2
#                         IG_V_gene                            lincRNA
#                                 3                                 22
#              processed_transcript                     protein_coding
#                                 4                                768
#                    sense_intronic                             snoRNA
#                                 2                                  1
#    transcribed_unitary_pseudogene transcribed_unprocessed_pseudogene
#                                 1                                  7
#                         TR_C_gene                          TR_J_gene
#                                 1                                 14
#                   TR_J_pseudogene                          TR_V_gene
#                                 1                                  9
#                   TR_V_pseudogene
#                                 1
unique(gtf_infor$gene_id) -> AFE_gene_list_equal1
unique(gtf_infor[gtf_infor$gene_biotype == "protein_coding",]$gene_id) -> AFE_pcgene_list_equal1

intersect(AFE_gene_list_more1, AFE_gene_list_equal1) -> AFE_gene_list_com
length(AFE_gene_list_com)
# [1] 155
AFE_gene_list_equal1[!AFE_gene_list_equal1 %in% AFE_gene_list_more1] -> AFE_gene_list_only_1cell

AFE_gene_list_more1[!AFE_gene_list_more1 %in% AFE_gene_list_equal1] -> AFE_gene_list_exactly_more_1cell
length(AFE_gene_list_only_1cell)
# [1] 716
length(AFE_gene_list_exactly_more_1cell)
# [1] 1049


for(i in c("AFE_gene_list_more1","AFE_pcgene_list_more1","AFE_gene_list_equal1","AFE_pcgene_list_equal1","AFE_gene_list_com","AFE_gene_list_only_1cell","AFE_gene_list_exactly_more_1cell")){
  ## enrichment
  eval(parse(text = paste("genelist <- ", i, sep ="")))
  ego_BP <- enrichGO(gene = genelist,
                     OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                     ont  = "BP", pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
  ego_MF <- enrichGO(gene = genelist,
                     OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                     ont  = "MF", pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
  ego_CC <- enrichGO(gene = genelist,
                     OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                     ont  = "CC", pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
  if(length(ego_BP$geneID) > 0 ){
  	pdf(paste("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/AFE/", i, "_BP.pdf", sep = ""),height = 10,width = 15)
  	print(dotplot(ego_BP,showCategory = 10,font.size = 16,title = "BP"))
  	plotGOgraph(ego_BP,firstSigNodes = 15, useInfo = "all")
  	enrichMap(ego_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
  	dev.off()
  }
  if(length(ego_MF$geneID) > 0 ){
  	pdf(paste("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/AFE/", i, "_MF.pdf", sep = ""),height = 10,width = 15)
  	print(dotplot(ego_MF,showCategory = 10,font.size = 16,title = "MF"))
  	plotGOgraph(ego_MF,firstSigNodes = 15, useInfo = "all")
  	enrichMap(ego_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
  	dev.off()
  }
  if(length(ego_CC$geneID) > 0 ){
  	pdf(paste("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/AFE/", i, "_CC.pdf", sep = ""),height = 10,width = 15)
  	print(dotplot(ego_CC,showCategory = 10,font.size = 16,title = "CC"))
  	plotGOgraph(ego_CC,firstSigNodes = 15, useInfo = "all")
  	enrichMap(ego_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = .4)
  	dev.off()
  }  
}


pdf(paste("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/AFE/", i, "_MF.pdf", sep = ""),height = 10,width = 15)
  print(dotplot(ego_MF,showCategory = 10,font.size = 16,title = "MF"))
  plotGOgraph(ego_MF,firstSigNodes = 15, useInfo = "all")
  enrichMap(ego_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
  dev.off()



## ALE -----------------------------------------------------------------------------------------------------------------------------------------------------


## ALE
load("./ALE/RData/all_cells_ALE_PSI_GTF.RData")
dim(all_cells_ALE_PSI)
# [1] 1211 3584
rowSums(!is.na(all_cells_ALE_PSI[,11:ncol(all_cells_ALE_PSI)])) -> ALE_psi_rowsum
sum(ALE_psi_rowsum>1)
# [1] 624

all_cells_ALE_PSI$name[ALE_psi_rowsum>1] -> ALE_loci_more1
all_cells_ALE_PSI$name[ALE_psi_rowsum==1] -> ALE_loci_equal1
length(ALE_loci_more1)
# [1] 624
length(ALE_loci_equal1)
# [1] 587


test <- as.data.frame(do.call(rbind, strsplit(na.omit(ALE_GTF[ALE_GTF$name %in% ALE_loci_more1,]$V9), split = "[; ]")))

gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
                        transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
                        gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
                        transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
                        gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
                        transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 857    6
length(unique(gtf_infor$gene_id))
# [1] 566
length(unique(gtf_infor$transcript_id))
# [1] 857

table(gtf_infor$gene_biotype)
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype)
#                         antisense                          IG_C_gene
#                                11                                  1
#                           lincRNA               processed_transcript
#                                17                                  5
#                    protein_coding transcribed_unprocessed_pseudogene
#                               525                                  6
#                unitary_pseudogene
#                                 1
unique(gtf_infor$gene_id) -> ALE_gene_list_more1
unique(gtf_infor[gtf_infor$gene_biotype == "protein_coding",]$gene_id) -> ALE_pcgene_list_more1


test <- as.data.frame(do.call(rbind, strsplit(na.omit(ALE_GTF[ALE_GTF$name %in% ALE_loci_equal1,]$V9), split = "[; ]")))

gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
                        transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
                        gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
                        transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
                        gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
                        transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 846    6
length(unique(gtf_infor$gene_id))
# [1] 545
length(unique(gtf_infor$transcript_id))
# [1] 846

table(gtf_infor$gene_biotype)
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype)
#          3prime_overlapping_ncRNA                          antisense
#                                 1                                 17
#                         IG_C_gene                            lincRNA
#                                 1                                 26
#              processed_transcript                     protein_coding
#                                 9                                474
#                    sense_intronic   transcribed_processed_pseudogene
#                                 2                                  2
#transcribed_unprocessed_pseudogene                          TR_C_gene
#                                 7                                  1
#                         TR_V_gene             unprocessed_pseudogene
#                                 1                                  4

unique(gtf_infor$gene_id) -> ALE_gene_list_equal1
unique(gtf_infor[gtf_infor$gene_biotype == "protein_coding",]$gene_id) -> ALE_pcgene_list_equal1

intersect(ALE_gene_list_more1, ALE_gene_list_equal1) -> ALE_gene_list_com
length(ALE_gene_list_com)
# [1] 78
ALE_gene_list_equal1[!ALE_gene_list_equal1 %in% ALE_gene_list_more1] -> ALE_gene_list_only_1cell

ALE_gene_list_more1[!ALE_gene_list_more1 %in% ALE_gene_list_equal1] -> ALE_gene_list_exactly_more_1cell
length(ALE_gene_list_only_1cell)
# [1] 467
length(ALE_gene_list_exactly_more_1cell)
# [1] 488


for(i in c("ALE_gene_list_more1","ALE_pcgene_list_more1","ALE_gene_list_equal1","ALE_pcgene_list_equal1","ALE_gene_list_com","ALE_gene_list_only_1cell","ALE_gene_list_exactly_more_1cell")){
  ## enrichment
  eval(parse(text = paste("genelist <- ", i, sep ="")))
  ego_BP <- enrichGO(gene = genelist,
                     OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                     ont  = "BP", pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
  ego_MF <- enrichGO(gene = genelist,
                     OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                     ont  = "MF", pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
  ego_CC <- enrichGO(gene = genelist,
                     OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                     ont  = "CC", pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
  
  if(length(ego_BP$geneID) > 0){
  	pdf(paste("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/ALE/", i, "_BP.pdf", sep = ""),height = 10,width = 15)
  	print(dotplot(ego_BP,showCategory = 10,font.size = 16,title = "BP"))
  	plotGOgraph(ego_BP,firstSigNodes = 15, useInfo = "all")
  	enrichMap(ego_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
  	dev.off()
  }
  if(length(ego_MF$geneID) > 0){
  	pdf(paste("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/ALE/", i, "_MF.pdf", sep = ""),height = 10,width = 15)
  	print(dotplot(ego_MF,showCategory = 10,font.size = 16,title = "MF"))
  	plotGOgraph(ego_MF,firstSigNodes = 15, useInfo = "all")
  	enrichMap(ego_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
  	dev.off()
  }
  if(length(ego_CC$geneID) > 0){
  	pdf(paste("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/ALE/", i, "_CC.pdf", sep = ""),height = 10,width = 15)
  	print(dotplot(ego_CC,showCategory = 10,font.size = 16,title = "CC"))
  	plotGOgraph(ego_CC,firstSigNodes = 15, useInfo = "all")
  	enrichMap(ego_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = .4)
  	dev.off()
  }  
}






#### gene biotype statistic **********************************************************************************************************************************************


setwd("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/")
load("./SE/RData/all_cells_SE_type_and_PSI.RData")
load("./SE/RData/all_cells_SE_GTF.RData")
table(all_cells_SE_type_and_PSI$AStype)
#      exon-skipping_exactly     exon-skipping_half-exon
#                       8311                         520
#exon-skipping_multipul-exon    intergenic-splicing_site
#                       3943                        2747
#                      Other
#                        324

all_cells_SE_type_and_PSI[all_cells_SE_type_and_PSI$AStype == "exon-skipping_exactly", ]$loci -> se_loci
dim(SE_GTF[SE_GTF$loci %in% se_loci & !is.na(SE_GTF$V9),])

test <- as.data.frame(do.call(rbind, strsplit(na.omit(SE_GTF[SE_GTF$loci %in% se_loci,]$V9), split = "[; ]")))

gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
                        transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
                        gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
                        transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
                        gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
                        transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 12981    6
length(unique(gtf_infor$gene_id))
# [1] 4609
length(unique(gtf_infor$transcript_id))
# [1] 12981

table(gtf_infor$gene_biotype)
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype)
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype) -> SE_gene_biotype_table
gtf_infor -> se_gtf_infor


## A3SS
load("./A3SS/RData/Parsed_A3SS_PSI_GTF.RData")
dim(A3SS_type_and_PSI)
# [1] 5874 3584
rowSums(!is.na(A3SS_type_and_PSI[,11:ncol(A3SS_type_and_PSI)])) -> a3ss_psi_rowsum
sum(a3ss_psi_rowsum>1)
# [1] 3179

test <- as.data.frame(do.call(rbind, strsplit(na.omit(A3SS_GTF$V9), split = "[; ]")))

gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
                        transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
                        gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
                        transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
                        gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
                        transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 8073    6
length(unique(gtf_infor$gene_id))
# [1] 3603
length(unique(gtf_infor$transcript_id))
# [1] 8073

table(gtf_infor$gene_biotype)
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype)
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype) -> A3SS_gene_biotype_table
gtf_infor -> a3ss_gtf_infor

## A5SS
load("./A5SS/RData/Parsed_A5SS_PSI_GTF.RData")
dim(A5SS_type_and_PSI)
# [1] 3908 3584
rowSums(!is.na(A5SS_type_and_PSI[,11:ncol(A5SS_type_and_PSI)])) -> a5ss_psi_rowsum
sum(a5ss_psi_rowsum>1)
# [1] 2256

test <- as.data.frame(do.call(rbind, strsplit(na.omit(A5SS_GTF$V9), split = "[; ]")))
gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
                        transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
                        gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
                        transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
                        gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
                        transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 5685    6
length(unique(gtf_infor$gene_id))
# [1] 2799
length(unique(gtf_infor$transcript_id))
# [1] 5685

table(gtf_infor$gene_biotype)
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype)
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype) -> A5SS_gene_biotype_table
gtf_infor -> a5ss_gtf_infor

## AFE
load("./AFE/RData/all_cells_AFE_PSI_GTF.RData")
dim(all_cells_AFE_PSI)
# [1] 2289 3584
rowSums(!is.na(all_cells_AFE_PSI[,11:ncol(all_cells_AFE_PSI)])) -> afe_psi_rowsum
sum(afe_psi_rowsum>1)
# [1] 1361

test <- as.data.frame(do.call(rbind, strsplit(na.omit(AFE_GTF$V9), split = "[; ]")))
gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
                        transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
                        gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
                        transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
                        gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
                        transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 2668    6
length(unique(gtf_infor$gene_id))
# [1] 1920
length(unique(gtf_infor$transcript_id))
# [1] 2668

table(gtf_infor$gene_biotype)
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype)
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype) -> AFE_gene_biotype_table
gtf_infor -> afe_gtf_infor



## ALE
load("./ALE/RData/all_cells_ALE_PSI_GTF.RData")
dim(all_cells_ALE_PSI)
# [1] 1211 3584
rowSums(!is.na(all_cells_ALE_PSI[,11:ncol(all_cells_ALE_PSI)])) -> ale_psi_rowsum
sum(ale_psi_rowsum>1)
# [1] 624

test <- as.data.frame(do.call(rbind, strsplit(na.omit(ALE_GTF$V9), split = "[; ]")))
gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
                        transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
                        gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
                        transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
                        gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
                        transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 1634    6
length(unique(gtf_infor$gene_id))
# [1] 1033
length(unique(gtf_infor$transcript_id))
# [1] 1634

table(gtf_infor$gene_biotype)
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype)
table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype) -> AFE_gene_biotype_table
gtf_infor -> ale_gtf_infor

#### IR
load("/mnt/data5/BGI/UCB/tangchao/IR/IR_GTF_match.RData")
test <- as.data.frame(do.call(rbind, strsplit(as.character(IR_GTF[,"V9"]), split = "[; ]")))
IR_gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
           transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
           gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
           transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
           gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
           transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})))


dim(IR_gtf_infor)
# [1] 4008    6
length(unique(IR_gtf_infor$gene_id))
# [1] 868
length(unique(IR_gtf_infor$transcript_id))
# [1] 2629

dim(IR_gtf_infor[!duplicated(IR_gtf_infor$gene_id),])
# [1] 868    6

table(IR_gtf_infor[!duplicated(IR_gtf_infor$gene_id),]$gene_biotype) -> IR_gene_biotype_table



AS_biotype <- rbind(se_gtf_infor[!duplicated(se_gtf_infor$gene_id), c("gene_id", "gene_biotype")],
	  a3ss_gtf_infor[!duplicated(a3ss_gtf_infor$gene_id), c("gene_id", "gene_biotype")],
	  a5ss_gtf_infor[!duplicated(a5ss_gtf_infor$gene_id), c("gene_id", "gene_biotype")],
	  afe_gtf_infor[!duplicated(afe_gtf_infor$gene_id), c("gene_id", "gene_biotype")],
	  ale_gtf_infor[!duplicated(ale_gtf_infor$gene_id), c("gene_id", "gene_biotype")],
	  IR_gtf_infor[!duplicated(IR_gtf_infor$gene_id), c("gene_id", "gene_biotype")])

AS_biotype$AStype <- rep(c("SE", "A3SS", "A5SS","AFE","ALE","IR"), 
	c(length(unique(se_gtf_infor$gene_id)),length(unique(a3ss_gtf_infor$gene_id)),
		length(unique(a5ss_gtf_infor$gene_id)),length(unique(afe_gtf_infor$gene_id)),
		length(unique(ale_gtf_infor$gene_id)),length(unique(IR_gtf_infor$gene_id))))


dim(AS_biotype)
table(AS_biotype$AStype, AS_biotype$gene_biotype)
table(AS_biotype$gene_biotype, AS_biotype$AStype)
# write.csv(as.matrix(table(AS_biotype$gene_biotype, AS_biotype$AStype)),"~/test.csv")




#### Pie chart of AS type ==============================================================================================================================
labs <- paste(c(4052,2685,2583,1099,60),"(",round(c(4052,2685,2583,1099,60)/sum(c(4052,2685,2583,1099,60))*100,2),"%)", sep = "")

pdf("/mnt/data5/BGI/UCB/tangchao/AS_stat/AS_type_pie_chart.pdf")
pie(c(4052,2685,2583,1099,60), labels = labs, main = "pie chart of AS types", col = rainbow(5))
legend("topright", c("SE","A3SS","A5SS","IR","MXE"), cex = 0.8, fill = rainbow(5), bty = "n")
dev.off()

n <- c(8311,5874,3908,2289,1211,376,1099)
labs <- paste(n, "(", round(n/sum(n)*100, 2),"%)", sep = "")

pdf("/mnt/data5/BGI/UCB/tangchao/AS_stat/AS_type(cell_by_cell)_pie_chart.pdf")
pie(n, labels = labs, main = "pie chart of AS types", col = rainbow(7))
legend("topright", c("SE","A3SS","A5SS","AFE","ALE","MXE","IR"), cex = 0.8, fill = rainbow(7), bty = "n")
dev.off()







#### Fisher/prop test of AS type in ---------------------------------------------------------------------------------------------------------------------------




pfgtf <- file.path("/mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.87.gtf")
gtf <- read.table(pfgtf,head=F,sep="\t")

gene <- subset(gtf,gtf$V3=="gene")

test <- as.data.frame(do.call(rbind, strsplit(as.character(gene$V9), split = "[; ]")))


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
           gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
           gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})), stringsAsFactors=F)


gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 58051     3


ReadCount <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/ReadCount_mat.noERCC.txt", sep = "\t", header = T, row.names = 1)
dim(ReadCount)
# [1] 34581  3574
TPM <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/TPM_mat.txt", sep = "\t", header = T, row.names = 1)
dim(TPM)
# [1] 58051  3574

gtf_infor[gtf_infor$gene_id %in% row.names(ReadCount),] -> gtf_infor_sub
dim(gtf_infor_sub)
# [1] 34581     3
table(gtf_infor_sub$gene_biotype)

          3prime_overlapping_ncRNA                          antisense
                                24                               3348
     bidirectional_promoter_lncRNA                          IG_C_gene
                                 1                                 14
                   IG_C_pseudogene                          IG_J_gene
                                 8                                  3
                     IG_pseudogene                          IG_V_gene
                                 1                                144
                   IG_V_pseudogene                            lincRNA
                               101                               3196
                      macro_lncRNA                              miRNA
                                 1                                 20
                          misc_RNA                            Mt_rRNA
                               563                                  2
                        non_coding             polymorphic_pseudogene
                                 3                                 10
              processed_pseudogene               processed_transcript
                              6287                                409
                    protein_coding                         pseudogene
                             16592                                 15
                          ribozyme                               rRNA
                                 2                                 77
                            scaRNA                              scRNA
                                18                                  1
                    sense_intronic                  sense_overlapping
                               504                                132
                            snoRNA                              snRNA
                               188                                316
                               TEC   transcribed_processed_pseudogene
                               613                                342
    transcribed_unitary_pseudogene transcribed_unprocessed_pseudogene
                                35                                544
                         TR_C_gene                          TR_V_gene
                                 6                                105
                   TR_V_pseudogene                 unitary_pseudogene
                                15                                 43
            unprocessed_pseudogene                           vaultRNA
                               897                                  1




## A3SS
fisher.test(matrix(c(60,5526,3551-60,34270-5526),nrow = 2))
fisher.test(matrix(c(55,7533,3551-55,34270-7533),nrow = 2))
fisher.test(matrix(c(31,515,3551-31,34270-515),nrow = 2))
fisher.test(matrix(c(3376,19961,3551-3376,34270-19961),nrow = 2))
fisher.test(matrix(c(33,735,3551-33,34270-735),nrow = 2))
## A5SS
fisher.test(matrix(c(56,5526,2771-56,34270-5526),nrow = 2))
fisher.test(matrix(c(51,7533,2771-51,34270-7533),nrow = 2))
fisher.test(matrix(c(41,515,2771-41,34270-515),nrow = 2))
fisher.test(matrix(c(2592,19961,2771-2592,34270-19961),nrow = 2))
fisher.test(matrix(c(31,735,2771-31,34270-735),nrow = 2))


## A5SS
fisher.test(matrix(c(56,   3348,  2771-56,   24089-3348),nrow = 2))
fisher.test(matrix(c(51,   3196,  2771-51,   24089-3196),nrow = 2))
fisher.test(matrix(c(41,   409,   2771-41,   24089-409),nrow = 2))
fisher.test(matrix(c(2592, 16592, 2771-2592, 24089-16592),nrow = 2))
fisher.test(matrix(c(31,   544,   2771-31,   24089-544),nrow = 2))


prop.test(matrix(c(56,   3348,  2771-56,   24089-3348),nrow = 2))
prop.test(matrix(c(51,   3196,  2771-51,   24089-3196),nrow = 2))
prop.test(matrix(c(41,   409,   2771-41,   24089-409),nrow = 2))
prop.test(matrix(c(2592, 16592, 2771-2592, 24089-16592),nrow = 2))
prop.test(matrix(c(31,   544,   2771-31,   24089-544),nrow = 2))


###################### prop.test of gene biotype in 6 AS types
c(60,56,60,27,26,48) -> antisense
c(55,51,41,41,5,47) -> lincRNA
c(31,41,17,12,10,51) -> processed
c(3376,2592,1728,926,781,4399) -> protein
c(33,31,14,13,12,38) -> unprocessed
c(3555,2771,1860,1019,834,4583) -> sum

prop.test(antisense, sum)
prop.test(lincRNA, sum)
prop.test(processed, sum)
prop.test(protein, sum)
prop.test(unprocessed, sum)










#### IR GO ===================================================================================================================================================



## IR

load("/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/IR_GTF_merge/IR_known_introns_match.RData")
length(unique(IR_GTF$gene_id))
# [1] 2320
table(IR_GTF[!duplicated(IR_GTF$gene_id),]$gene_biotype)
#          3prime_overlapping_ncRNA                          antisense
#                                 1                                134
#                         IG_C_gene                          IG_V_gene
#                                 2                                 10
#                   IG_V_pseudogene                            lincRNA
#                                 1                                 57
#              processed_pseudogene               processed_transcript
#                                 1                                 30
#                    protein_coding                     sense_intronic
#                              2029                                 12
#                 sense_overlapping   transcribed_processed_pseudogene
#                                 2                                  1
#transcribed_unprocessed_pseudogene                          TR_C_gene
#                                15                                  3
#                         TR_V_gene                    TR_V_pseudogene
#                                18                                  1
#            unprocessed_pseudogene
#                                 3

unique(IR_GTF$gene_id) -> IR_gene_list

for(i in c("IR_gene_list")){
  ## enrichment
  eval(parse(text = paste("genelist <- ", i, sep ="")))
  ego_BP <- enrichGO(gene = genelist,
                     OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                     ont  = "BP", pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
  ego_MF <- enrichGO(gene = genelist,
                     OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                     ont  = "MF", pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
  ego_CC <- enrichGO(gene = genelist,
                     OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                     ont  = "CC", pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
  
  if(length(ego_BP$geneID) > 0){
  	pdf(paste("/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/figure/", i, "_BP.pdf", sep = ""),height = 10,width = 15)
  	print(dotplot(ego_BP,showCategory = 10,font.size = 16,title = "BP"))
  	plotGOgraph(ego_BP,firstSigNodes = 15, useInfo = "all")
  	enrichMap(ego_BP,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
  	dev.off()
  }
  if(length(ego_MF$geneID) > 0){
  	pdf(paste("/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/figure/", i, "_MF.pdf", sep = ""),height = 10,width = 15)
  	print(dotplot(ego_MF,showCategory = 10,font.size = 16,title = "MF"))
  	plotGOgraph(ego_MF,firstSigNodes = 15, useInfo = "all")
  	enrichMap(ego_MF,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
  	dev.off()
  }
  if(length(ego_CC$geneID) > 0){
  	pdf(paste("/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/figure/", i, "_CC.pdf", sep = ""),height = 10,width = 15)
  	print(dotplot(ego_CC,showCategory = 10,font.size = 16,title = "CC"))
  	plotGOgraph(ego_CC,firstSigNodes = 15, useInfo = "all")
  	enrichMap(ego_CC,n = 30, vertex.label.font = 1,cex = 1, font.size = .4)
  	dev.off()
  }  
}








