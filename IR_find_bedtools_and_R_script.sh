#!/bin/sh

#  bedtools for supplementary data of BGI UCB
#
#
#  Created by TangChao on 23/02/2018.
#
cd /mnt/data5/BGI/UCB/Star
#for f in $(ls )
#do
#bedtools genomecov -ibam /mnt/data5/BGI/UCB/Star/$f/*bam -bga -split | bedtools intersect -a /mnt/data1/projects/UCB/results_ucb/tangchao/IR/BEDtools/intron/intron.bed -b - -wa -wb > /mnt/data5/BGI/UCB/tangchao/IR/bedtools/${f}_introncov
#done
#
#exit 0

## generate genome coverage of every BAM
for f in $(ls )
do
bedtools genomecov -ibam /mnt/data5/BGI/UCB/Star/$f/*bam -bga -split > /mnt/data5/BGI/UCB/tangchao/IR/bedtools/genomecov/${f}_genomecov



## call all introns' range of whole genome

bedtools intersect -a /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/introns_touse.bed -b /mnt/data5/BGI/UCB/tangchao/IR/bedtools/genomecov/${f}_genomecov -wa -wb > /mnt/data5/BGI/UCB/tangchao/IR/bedtools/intron/${f}_introncov



## call all left exons' range of whole genome

bedtools intersect -a /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/left_exon.bed -b /mnt/data5/BGI/UCB/tangchao/IR/bedtools/genomecov/${f}_genomecov -wa -wb > /mnt/data5/BGI/UCB/tangchao/IR/bedtools/leftexon/${f}_leftexoncov



## call all right exons' range of whole genome

bedtools intersect -a /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/right_exon.bed -b /mnt/data5/BGI/UCB/tangchao/IR/bedtools/genomecov/${f}_genomecov -wa -wb > /mnt/data5/BGI/UCB/tangchao/IR/bedtools/rightexon/${f}_rightexoncov
done

exit 0





## R script
cd /mnt/data5/BGI/UCB/Star

for f in $(ls | sed -n '1,100p')
do
echo $f
Rscript --vanilla /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/MMCLGenomeRegion.R /mnt/data5/BGI/UCB/tangchao/IR/bedtools/intron/${f}_introncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/leftexon/${f}_leftexoncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/rightexon/${f}_rightexoncov /mnt/data5/BGI/UCB/tangchao/IR/Result/${f}.txt
done


for f in $(ls | sed -n '101,200p')
do
echo $f
Rscript --vanilla /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/MMCLGenomeRegion.R /mnt/data5/BGI/UCB/tangchao/IR/bedtools/intron/${f}_introncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/leftexon/${f}_leftexoncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/rightexon/${f}_rightexoncov /mnt/data5/BGI/UCB/tangchao/IR/Result/${f}.txt
done


for f in $(ls | sed -n '201,300p')
do
echo $f
Rscript --vanilla /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/MMCLGenomeRegion.R /mnt/data5/BGI/UCB/tangchao/IR/bedtools/intron/${f}_introncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/leftexon/${f}_leftexoncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/rightexon/${f}_rightexoncov /mnt/data5/BGI/UCB/tangchao/IR/Result/${f}.txt
done


for f in $(ls | sed -n '301,400p')
do
echo $f
Rscript --vanilla /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/MMCLGenomeRegion.R /mnt/data5/BGI/UCB/tangchao/IR/bedtools/intron/${f}_introncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/leftexon/${f}_leftexoncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/rightexon/${f}_rightexoncov /mnt/data5/BGI/UCB/tangchao/IR/Result/${f}.txt
done


for f in $(ls | sed -n '401,500p')
do
echo $f
Rscript --vanilla /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/MMCLGenomeRegion.R /mnt/data5/BGI/UCB/tangchao/IR/bedtools/intron/${f}_introncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/leftexon/${f}_leftexoncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/rightexon/${f}_rightexoncov /mnt/data5/BGI/UCB/tangchao/IR/Result/${f}.txt
done


for f in $(ls | sed -n '501,600p')
do
echo $f
Rscript --vanilla /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/MMCLGenomeRegion.R /mnt/data5/BGI/UCB/tangchao/IR/bedtools/intron/${f}_introncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/leftexon/${f}_leftexoncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/rightexon/${f}_rightexoncov /mnt/data5/BGI/UCB/tangchao/IR/Result/${f}.txt
done


for f in $(ls | sed -n '601,700p')
do
echo $f
Rscript --vanilla /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/MMCLGenomeRegion.R /mnt/data5/BGI/UCB/tangchao/IR/bedtools/intron/${f}_introncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/leftexon/${f}_leftexoncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/rightexon/${f}_rightexoncov /mnt/data5/BGI/UCB/tangchao/IR/Result/${f}.txt
done


for f in $(ls | sed -n '701,800p')
do
echo $f
Rscript --vanilla /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/MMCLGenomeRegion.R /mnt/data5/BGI/UCB/tangchao/IR/bedtools/intron/${f}_introncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/leftexon/${f}_leftexoncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/rightexon/${f}_rightexoncov /mnt/data5/BGI/UCB/tangchao/IR/Result/${f}.txt
done


for f in $(ls | sed -n '801,859p')
do
echo $f
Rscript --vanilla /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/MMCLGenomeRegion.R /mnt/data5/BGI/UCB/tangchao/IR/bedtools/intron/${f}_introncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/leftexon/${f}_leftexoncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/rightexon/${f}_rightexoncov /mnt/data5/BGI/UCB/tangchao/IR/Result/${f}.txt
done


for f in $(ls | sed -n '860,918p')
do
echo $f
Rscript --vanilla /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/MMCLGenomeRegion.R /mnt/data5/BGI/UCB/tangchao/IR/bedtools/intron/${f}_introncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/leftexon/${f}_leftexoncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/rightexon/${f}_rightexoncov /mnt/data5/BGI/UCB/tangchao/IR/Result/${f}.txt
done









