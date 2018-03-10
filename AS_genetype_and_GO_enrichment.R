setwd("/mnt/data5/BGI/UCB/tangchao/DSU/SE/RData")
load("./SE_DSU_10reads_2cells_result_20180228.RData")
load("../../A3SS/RData/A3SS_DSU_10reads_2cells_result_20180228.RData")
load("../../A5SS/RData/A5SS_DSU_10reads_2cells_result_20180228.RData")
load("../../MXE/RData/MXE_DSU_10reads_2cells_result_20171228.RData")

load("../../../IR/IR_GTF_match.RData")

table(SE_type$AStype)
#      exon-skipping_exactly     exon-skipping_half-exon
#                       4052                         123
#exon-skipping_multipul-exon    intergenic-splicing_site
#                        462                         973
#                      Other
#                        104


table(A3SS_type$AStype)
#
#                                          alt_5/3-splicing_site
#                                                           2685
#                                          exon-skipping_exactly
#                                                           2685
#                                    exon-skipping_multiple-exon
#                                                           4134
#exon-skipping_multiple-exon_and_first_and_last_exon_unannotated
#                                                            169
#         exon-skipping_multiple-exon_and_first_exon_unannotated
#                                                            280
#                     exon-skipping_skipped_one_unannotated_exon
#                                                           1228
#                               intergenic-alt_5/3-splicing_site
#                                                            391
#                                                         Others
#                                                           1905

table(A5SS_type$AStype)
#
#                                          alt_5/3-splicing_site
#                                                           2583
#                                                  exon-skipping
#                                                           2688
#                                    exon-skipping_multiple-exon
#                                                           4161
#exon-skipping_multiple-exon_and_first_and_last_exon_unannotated
#                                                            194
#          exon-skipping_multiple-exon_and_last_exon_unannotated
#                                                            281
#                     exon-skipping_skipped_one_unannotated_exon
#                                                           1147
#                               intergenic-alt_5/3-splicing_site
#                                                            347
#                                                         Others
#                                                           1819

table(MXE_type$AStype)
#         Mixed-splicing Mutually-exclusive-exon
#                    198                      60

table(IR_type$AStype)
#intron_retention
#             1099


labs <- paste(c(4052,2685,2583,1099,60),"(",round(c(4052,2685,2583,1099,60)/sum(c(4052,2685,2583,1099,60))*100,2),"%)", sep = "")

pdf("/mnt/data5/BGI/UCB/tangchao/AS_stat/AS_type_pie_chart.pdf")
pie(c(4052,2685,2583,1099,60), labels = labs, main = "pie chart of AS types", col = rainbow(5))
legend("topright", c("SE","A3SS","A5SS","IR","MXE"), cex = 0.8, fill = rainbow(5), bty = "n")
dev.off()



exactly_SE_GTF <- SE_GTF[SE_GTF$AStype == "exon-skipping_exactly" & SE_GTF$OverType == "equal",]

test <- as.data.frame(do.call(rbind, strsplit(exactly_SE_GTF[,24], split = "[; ]")))


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
           transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})))


dim(gtf_infor)
# [1] 4052    6
length(unique(gtf_infor$gene_id))
# [1] 3178
length(unique(gtf_infor$transcript_id))
# [1] 3557

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

table(gtf_infor[!duplicated(gtf_infor[,"gene_id"]),]$gene_biotype)
#                         antisense                          IG_C_gene 
#                                27                                  1 
#                           lincRNA             polymorphic_pseudogene 
#                                15                                  1 
#              processed_transcript                     protein_coding 
#                                22                               3084 
#                 sense_overlapping   transcribed_processed_pseudogene 
#                                 1                                  2 
#transcribed_unprocessed_pseudogene                 unitary_pseudogene 
#                                20                                  1 
#            unprocessed_pseudogene 
#                                 4 


save(gtf_infor, file = "/mnt/data5/BGI/UCB/tangchao/AS_stat/SE_exactly_gene_info.RData")



############################ SE GO

library(org.Hs.eg.db)#####datdbase target to human-geneIDs
library(AnnotationDbi)
biocLite("clusterProfiler")
library(DOSE)
library(GO.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)
#######for ENSEMBLID of genelist
unique(as.character(gtf_infor[gtf_infor$gene_biotype == "protein_coding","gene_id"]))

ego_pc_BP <- enrichGO(gene = unique(as.character(gtf_infor[gtf_infor$gene_biotype == "protein_coding","gene_id"])),
                OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                ont  = "BP", pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
ego_pc_MF <- enrichGO(gene = unique(as.character(gtf_infor[gtf_infor$gene_biotype == "protein_coding","gene_id"])),
                      OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                      ont  = "MF", pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
ego_pc_CC <- enrichGO(gene = unique(as.character(gtf_infor[gtf_infor$gene_biotype == "protein_coding","gene_id"])),
                      OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                      ont  = "CC", pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/AS_stat/SE_exactly_genes_ego_pc_BP.pdf",height = 10,width = 15)
dotplot(ego_pc_BP,showCategory = 10,font.size = 16,title = "ego_cell")
plotGOgraph(ego_pc_BP,firstSigNodes = 15, useInfo = "all")
enrichMap(ego_pc_BP,n = 30, vertex.label.font = 1,cex = 6, font.size = 1)
dev.off()
pdf("/mnt/data5/BGI/UCB/tangchao/AS_stat/SE_exactly_genes_ego_pc_MF.pdf",height = 10,width = 15)
dotplot(ego_pc_MF,showCategory = 10,font.size = 16,title = "ego_cell")
plotGOgraph(ego_pc_MF,firstSigNodes = 15, useInfo = "all")
enrichMap(ego_pc_MF,n = 30, vertex.label.font = 1,cex = 6, font.size = 1)
dev.off()
pdf("/mnt/data5/BGI/UCB/tangchao/AS_stat/SE_exactly_genes_ego_pc_CC.pdf",height = 10,width = 15)
dotplot(ego_pc_CC,showCategory = 10,font.size = 16,title = "ego_cell")
plotGOgraph(ego_pc_CC,firstSigNodes = 15, useInfo = "all")
enrichMap(ego_pc_CC,n = 30, vertex.label.font = 1,cex = 6, font.size = 1)
dev.off()



#######for ENSEMBLID of genelist
names(table(as.character(gtf_infor$gene_id)))[table(as.character(gtf_infor$gene_id)) > 1]

ego_MS_BP <- enrichGO(gene = names(table(as.character(gtf_infor$gene_id)))[table(as.character(gtf_infor$gene_id)) > 1],
                      OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                      ont  = "BP", pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
ego_MS_MF <- enrichGO(gene = names(table(as.character(gtf_infor$gene_id)))[table(as.character(gtf_infor$gene_id)) > 1],
                      OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                      ont  = "MF", pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
ego_MS_CC <- enrichGO(gene = names(table(as.character(gtf_infor$gene_id)))[table(as.character(gtf_infor$gene_id)) > 1],
                      OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                      ont  = "CC", pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)


pdf("/mnt/data5/BGI/UCB/tangchao/AS_stat/SE_exactly_genes_ego_MS_BP.pdf",height = 10,width = 15)
dotplot(ego_MS_BP,showCategory = 10,font.size = 16,title = "ego_cell")
plotGOgraph(ego_MS_BP,firstSigNodes = 15, useInfo = "all")
enrichMap(ego_MS_BP,n = 30, vertex.label.font = 1,cex = 6, font.size = 1)
dev.off()
pdf("/mnt/data5/BGI/UCB/tangchao/AS_stat/SE_exactly_genes_ego_MS_MF.pdf",height = 10,width = 15)
dotplot(ego_MS_MF,showCategory = 10,font.size = 16,title = "ego_cell")
plotGOgraph(ego_MS_MF,firstSigNodes = 15, useInfo = "all")
enrichMap(ego_MS_MF,n = 30, vertex.label.font = 1,cex = 6, font.size = 1)
dev.off()
pdf("/mnt/data5/BGI/UCB/tangchao/AS_stat/SE_exactly_genes_ego_MS_CC.pdf",height = 10,width = 15)
dotplot(ego_MS_CC,showCategory = 10,font.size = 16,title = "ego_cell")
plotGOgraph(ego_MS_CC,firstSigNodes = 15, useInfo = "all")
enrichMap(ego_MS_CC,n = 30, vertex.label.font = 1,cex = 6, font.size = 1)
dev.off()




############################ A3SS GO


table(A3SS_type$AStype)
exactly_A3SS_GTF <- A3SS_GTF[A3SS_GTF$AStype == "alt_5/3-splicing_site" & A3SS_GTF$OverType == "start",]

test <- as.data.frame(do.call(rbind, strsplit(exactly_A3SS_GTF[,"V9"], split = "[; ]")))


valueFind <- function(x, type = "gene_id"){
  for(i in 1:length(x)){
    if(sum(x == type) == 0){
      return(NA)
    }else{
      return(x[which(x == type)[1]+1])
    }
  }
}



a3_gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
           transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
           gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
           transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
           gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
           transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})))


dim(a3_gtf_infor)
# [1] 4393    6
length(unique(a3_gtf_infor$gene_id))
# [1] 2245
length(unique(a3_gtf_infor$transcript_id))
# [1] 4050

dim(a3_gtf_infor[!duplicated(a3_gtf_infor$gene_id),])
# [1] 2245    6

table(a3_gtf_infor[!duplicated(a3_gtf_infor$gene_id),]$gene_biotype)
#                          antisense                          IG_C_gene
#                                43                                  1
#                   IG_V_pseudogene                            lincRNA
#                                 1                                 31
#              processed_transcript                     protein_coding
#                                26                               2110
#                    sense_intronic                  sense_overlapping
#                                 1                                  1
#  transcribed_processed_pseudogene     transcribed_unitary_pseudogene
#                                 3                                  2
#transcribed_unprocessed_pseudogene                          TR_C_gene
#                                13                                  2
#                         TR_V_gene                    TR_V_pseudogene
#                                 5                                  1
#                unitary_pseudogene             unprocessed_pseudogene
#                                 1                                  4


dim(a3_gtf_infor[!duplicated(a3_gtf_infor$gene_id) & a3_gtf_infor$gene_biotype == "protein_coding",])

A3_pc_gene_list <- as.character(a3_gtf_infor[!duplicated(a3_gtf_infor$gene_id) & a3_gtf_infor$gene_biotype == "protein_coding",]$gene_id)


#######for ENSEMBLID of genelist

ego_PC_BP <- enrichGO(gene = A3_pc_gene_list,
                      OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                      ont  = "BP", pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
ego_PC_MF <- enrichGO(gene = A3_pc_gene_list,
                      OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                      ont  = "MF", pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
ego_PC_CC <- enrichGO(gene = A3_pc_gene_list,
                      OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                      ont  = "CC", pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)


pdf("/mnt/data5/BGI/UCB/tangchao/AS_stat/A3SS_PC_genes_ego_BP.pdf",height = 10,width = 15)
dotplot(ego_PC_BP,showCategory = 10,font.size = 16,title = "ego_cell")
plotGOgraph(ego_PC_BP,firstSigNodes = 15, useInfo = "all")
enrichMap(ego_PC_BP,n = 30, vertex.label.font = 1,cex = 6, font.size = 1)
dev.off()
pdf("/mnt/data5/BGI/UCB/tangchao/AS_stat/A3SS_PC_genes_ego_MF.pdf",height = 10,width = 15)
dotplot(ego_PC_MF,showCategory = 10,font.size = 16,title = "ego_cell")
plotGOgraph(ego_PC_MF,firstSigNodes = 15, useInfo = "all")
enrichMap(ego_PC_MF,n = 30, vertex.label.font = 1,cex = 6, font.size = 1)
dev.off()
pdf("/mnt/data5/BGI/UCB/tangchao/AS_stat/A3SS_PC_genes_ego_CC.pdf",height = 10,width = 15)
dotplot(ego_PC_CC,showCategory = 10,font.size = 16,title = "ego_cell")
plotGOgraph(ego_PC_CC,firstSigNodes = 15, useInfo = "all")
enrichMap(ego_PC_CC,n = 30, vertex.label.font = 1,cex = 6, font.size = 1)
dev.off()






############################ A5SS GO




exactly_A5SS_GTF <- A5SS_GTF[A5SS_GTF$AStype == "alt_5/3-splicing_site" & A5SS_GTF$OverType == "end",]

test <- as.data.frame(do.call(rbind, strsplit(exactly_A5SS_GTF[,"V9"], split = "[; ]")))


valueFind <- function(x, type = "gene_id"){
  for(i in 1:length(x)){
    if(sum(x == type) == 0){
      return(NA)
    }else{
      return(x[which(x == type)[1]+1])
    }
  }
}



a5_gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
           transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
           gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
           transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
           gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
           transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})))


dim(a5_gtf_infor)
# [1] 4321    6
length(unique(a5_gtf_infor$gene_id))
# [1] 2174
length(unique(a5_gtf_infor$transcript_id))
# [1] 3939

dim(a5_gtf_infor[!duplicated(a5_gtf_infor$gene_id),])
# [1] 2174    6

table(a5_gtf_infor[!duplicated(a5_gtf_infor$gene_id),]$gene_biotype)
#                        antisense                            lincRNA
#                                29                                 31
#                          misc_RNA               processed_transcript
#                                 1                                 22
#                    protein_coding   transcribed_processed_pseudogene
#                              2061                                  2
#transcribed_unprocessed_pseudogene                          TR_C_gene
#                                19                                  1
#                         TR_V_gene                 unitary_pseudogene
#                                 4                                  2
#            unprocessed_pseudogene
#                                 2

dim(a5_gtf_infor[!duplicated(a5_gtf_infor$gene_id) & a5_gtf_infor$gene_biotype == "protein_coding",])

A5_pc_gene_list <- as.character(a5_gtf_infor[!duplicated(a5_gtf_infor$gene_id) & a5_gtf_infor$gene_biotype == "protein_coding",]$gene_id)


#######for ENSEMBLID of genelist

ego_PC_BP <- enrichGO(gene = A5_pc_gene_list,
                      OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                      ont  = "BP", pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
ego_PC_MF <- enrichGO(gene = A5_pc_gene_list,
                      OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                      ont  = "MF", pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
ego_PC_CC <- enrichGO(gene = A5_pc_gene_list,
                      OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                      ont  = "CC", pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)


pdf("/mnt/data5/BGI/UCB/tangchao/AS_stat/A5SS_PC_genes_ego_BP.pdf",height = 10,width = 15)
dotplot(ego_PC_BP,showCategory = 10,font.size = 16,title = "ego_cell")
plotGOgraph(ego_PC_BP,firstSigNodes = 15, useInfo = "all")
enrichMap(ego_PC_BP,n = 30, vertex.label.font = 1,cex = 6, font.size = 1)
dev.off()
pdf("/mnt/data5/BGI/UCB/tangchao/AS_stat/A5SS_PC_genes_ego_MF.pdf",height = 10,width = 15)
dotplot(ego_PC_MF,showCategory = 10,font.size = 16,title = "ego_cell")
plotGOgraph(ego_PC_MF,firstSigNodes = 15, useInfo = "all")
enrichMap(ego_PC_MF,n = 30, vertex.label.font = 1,cex = 6, font.size = 1)
dev.off()
pdf("/mnt/data5/BGI/UCB/tangchao/AS_stat/A5SS_PC_genes_ego_CC.pdf",height = 10,width = 15)
dotplot(ego_PC_CC,showCategory = 10,font.size = 16,title = "ego_cell")
plotGOgraph(ego_PC_CC,firstSigNodes = 15, useInfo = "all")
enrichMap(ego_PC_CC,n = 30, vertex.label.font = 1,cex = 6, font.size = 1)
dev.off()






############################ IR GO



table(IR_type$AStype)
# intron_retention
#            1099

table(IR_GTF$AStype)
# intron_retention
#             4008

test <- as.data.frame(do.call(rbind, strsplit(as.character(IR_GTF[,"V9"]), split = "[; ]")))


valueFind <- function(x, type = "gene_id"){
  for(i in 1:length(x)){
    if(sum(x == type) == 0){
      return(NA)
    }else{
      return(x[which(x == type)[1]+1])
    }
  }
}



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

table(IR_gtf_infor[!duplicated(IR_gtf_infor$gene_id),]$gene_biotype)
#                         antisense                          IG_C_gene
#                                26                                  2
#                         IG_V_gene                            lincRNA
#                                 4                                  5
#              processed_pseudogene               processed_transcript
#                                 6                                 10
#                    protein_coding                         pseudogene
#                               781                                  1
#                    sense_intronic                  sense_overlapping
#                                 6                                  1
#transcribed_unprocessed_pseudogene                          TR_C_gene
#                                12                                  4
#                         TR_V_gene             unprocessed_pseudogene
#                                 6                                  4

dim(IR_gtf_infor[!duplicated(IR_gtf_infor$gene_id) & IR_gtf_infor$gene_biotype == "protein_coding",])

IR_pc_gene_list <- as.character(IR_gtf_infor[!duplicated(IR_gtf_infor$gene_id) & IR_gtf_infor$gene_biotype == "protein_coding",]$gene_id)



#######for ENSEMBLID of genelist

ego_PC_BP <- enrichGO(gene = IR_pc_gene_list,
                      OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                      ont  = "BP", pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
ego_PC_MF <- enrichGO(gene = IR_pc_gene_list,
                      OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                      ont  = "MF", pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
ego_PC_CC <- enrichGO(gene = IR_pc_gene_list,
                      OrgDb  = org.Hs.eg.db,keyType= 'ENSEMBL',
                      ont  = "CC", pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)


pdf("/mnt/data5/BGI/UCB/tangchao/AS_stat/IR_PC_genes_ego_BP.pdf",height = 10,width = 15)
dotplot(ego_PC_BP,showCategory = 10,font.size = 16,title = "ego_cell")
plotGOgraph(ego_PC_BP,firstSigNodes = 15, useInfo = "all")
enrichMap(ego_PC_BP,n = 30, vertex.label.font = 1,cex = 6, font.size = 1)
dev.off()
pdf("/mnt/data5/BGI/UCB/tangchao/AS_stat/IR_PC_genes_ego_MF.pdf",height = 10,width = 15)
dotplot(ego_PC_MF,showCategory = 10,font.size = 16,title = "ego_cell")
plotGOgraph(ego_PC_MF,firstSigNodes = 15, useInfo = "all")
enrichMap(ego_PC_MF,n = 30, vertex.label.font = 1,cex = 6, font.size = 1)
dev.off()
pdf("/mnt/data5/BGI/UCB/tangchao/AS_stat/IR_PC_genes_ego_CC.pdf",height = 10,width = 15)
dotplot(ego_PC_CC,showCategory = 10,font.size = 16,title = "ego_cell")
plotGOgraph(ego_PC_CC,firstSigNodes = 15, useInfo = "all")
enrichMap(ego_PC_CC,n = 30, vertex.label.font = 1,cex = 6, font.size = 1)
dev.off()




length(unique(gtf_infor$gene_id))
# 3178
length(unique(a3_gtf_infor$gene_id))
# 2245
length(unique(a5_gtf_infor$gene_id))
# 2174
length(unique(IR_gtf_infor$gene_id))
# 868


dim(gtf_infor[!duplicated(gtf_infor$gene_id), c("gene_id", "gene_biotype")])
# [1] 3178    2
dim(a3_gtf_infor[!duplicated(a3_gtf_infor$gene_id), c("gene_id", "gene_biotype")])
# [1] 2245    2
dim(a5_gtf_infor[!duplicated(a5_gtf_infor$gene_id), c("gene_id", "gene_biotype")])
# [1] 2174    2
dim(IR_gtf_infor[!duplicated(IR_gtf_infor$gene_id), c("gene_id", "gene_biotype")])
# [1] 868    2

AS_biotype <- rbind(gtf_infor[!duplicated(gtf_infor$gene_id), c("gene_id", "gene_biotype")],
	  a3_gtf_infor[!duplicated(a3_gtf_infor$gene_id), c("gene_id", "gene_biotype")],
	  a5_gtf_infor[!duplicated(a5_gtf_infor$gene_id), c("gene_id", "gene_biotype")],
	  IR_gtf_infor[!duplicated(IR_gtf_infor$gene_id), c("gene_id", "gene_biotype")])

AS_biotype$AStype <- rep(c("SE", "A3SS", "A5SS","IR"), c(length(unique(gtf_infor$gene_id)),length(unique(a3_gtf_infor$gene_id)),length(unique(a5_gtf_infor$gene_id)),length(unique(IR_gtf_infor$gene_id))))


dim(AS_biotype)
table(AS_biotype$AStype, AS_biotype$gene_biotype)
table(AS_biotype$gene_biotype, AS_biotype$AStype)
#                                     A3SS A5SS   IR   SE
#  antisense                            43   29   26   27
#  IG_C_gene                             1    0    2    1
#  lincRNA                              31   31    5   15
#  polymorphic_pseudogene                0    0    0    1
#  processed_transcript                 26   22   10   22
#  protein_coding                     2110 2061  781 3084
#  sense_overlapping                     1    0    1    1
#  transcribed_processed_pseudogene      3    2    0    2
#  transcribed_unprocessed_pseudogene   13   19   12   20
#  unitary_pseudogene                    1    2    0    1
#  unprocessed_pseudogene                4    2    4    4
#  IG_V_pseudogene                       1    0    0    0
#  sense_intronic                        1    0    6    0
#  transcribed_unitary_pseudogene        2    0    0    0
#  TR_C_gene                             2    1    4    0
#  TR_V_gene                             5    4    6    0
#  TR_V_pseudogene                       1    0    0    0
#  misc_RNA                              0    1    0    0
#  IG_V_gene                             0    0    4    0
#  processed_pseudogene                  0    0    6    0
#  pseudogene                            0    0    1    0

AS_biotype_table <- table(AS_biotype$gene_biotype, AS_biotype$AStype)

AS_biotype_table["protein_coding",]/colSums(AS_biotype_table)
#      A3SS      A5SS        IR        SE
# 0.9398664 0.9480221 0.8997696 0.9704216

pdf("/mnt/data5/BGI/UCB/tangchao/AS_stat/Percentage_of_PC_genes.pdf")
p = AS_biotype_table["protein_coding",]/colSums(AS_biotype_table)*100
x = barplot(p, ylim = c(0,105), main = "Number and percentage of protein coding gene\nin different AS types", col = rainbow(4))
labs <- paste(AS_biotype_table["protein_coding",],"\n",round(p,2),"%", sep = "")
text(x = x,y = p + 4, labels = labs)
dev.off()























