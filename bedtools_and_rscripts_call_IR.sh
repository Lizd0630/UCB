#!/bin/sh

#  bedtools for supplementary data of BGI UCB
#  We used the SJ.out.tab file of STAR to identify the introns range.
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
#bedtools genomecov -ibam /mnt/data5/BGI/UCB/Star/$f/*bam -bga -split > /mnt/data5/BGI/UCB/tangchao/IR/bedtools/genomecov/${f}_genomecov
# I have completed this step before.


## call all introns' range of whole genome

awk '$7>1' /mnt/data5/BGI/UCB/Star/$f/*SJ.out.tab | awk '{print $1"\t"$2"\t"$3}' | bedtools intersect -a - -b /mnt/data5/BGI/UCB/tangchao/IR/bedtools/genomecov/${f}_genomecov -wa -wb > /mnt/data5/BGI/UCB/tangchao/IR/bedtools/intron/${f}_introncov



## call all left exons' range of whole genome

awk '$7>1' /mnt/data5/BGI/UCB/Star/$f/*SJ.out.tab | awk '{print $1"\t"$2-100"\t"$2-1"\t"$1"\t"$2"\t"$3}' | bedtools intersect -a - -b /mnt/data5/BGI/UCB/tangchao/IR/bedtools/genomecov/${f}_genomecov -wa -wb > /mnt/data5/BGI/UCB/tangchao/IR/bedtools/leftexon/${f}_leftexoncov



## call all right exons' range of whole genome

awk '$7>1' /mnt/data5/BGI/UCB/Star/$f/*SJ.out.tab | awk '{print $1"\t"$3+1"\t"$3+100"\t"$1"\t"$2"\t"$3}' | bedtools intersect -a - -b /mnt/data5/BGI/UCB/tangchao/IR/bedtools/genomecov/${f}_genomecov -wa -wb > /mnt/data5/BGI/UCB/tangchao/IR/bedtools/rightexon/${f}_rightexoncov

done

exit 0




## R script
cd /mnt/data5/BGI/UCB/Star

for f in $(ls)
do
echo $f
Rscript --vanilla /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/MMCLGenomeRegion.R /mnt/data5/BGI/UCB/tangchao/IR/bedtools/intron/${f}_introncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/leftexon/${f}_leftexoncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/rightexon/${f}_rightexoncov /mnt/data5/BGI/UCB/tangchao/IR/Result/${f}.txt
done






## data before
cd /mnt/data1/projects/UCB/01.Star


for f in $(ls )
do
#bedtools genomecov -ibam /mnt/data1/projects/UCB/01.Star/$f/*bam -bga -split > /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/genomecov/${f}_genomecov
# I have completed this step before.


## call all introns' range of whole genome

awk '$7>1' /mnt/data1/projects/UCB/01.Star/$f/*SJ.out.tab | awk '{print $1"\t"$2"\t"$3}' | bedtools intersect -a - -b /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/genomecov/${f}_genomecov -wa -wb > /mnt/data5/BGI/UCB/tangchao/IR/bedtools/intron/${f}_introncov



## call all left exons' range of whole genome

awk '$7>1' /mnt/data1/projects/UCB/01.Star/$f/*SJ.out.tab | awk '{print $1"\t"$2-100"\t"$2-1"\t"$1"\t"$2"\t"$3}' | bedtools intersect -a - -b /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/genomecov/${f}_genomecov -wa -wb > /mnt/data5/BGI/UCB/tangchao/IR/bedtools/leftexon/${f}_leftexoncov



## call all right exons' range of whole genome

awk '$7>1' /mnt/data1/projects/UCB/01.Star/$f/*SJ.out.tab | awk '{print $1"\t"$3+1"\t"$3+100"\t"$1"\t"$2"\t"$3}' | bedtools intersect -a - -b /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/genomecov/${f}_genomecov -wa -wb > /mnt/data5/BGI/UCB/tangchao/IR/bedtools/rightexon/${f}_rightexoncov

done


## R script
cd /mnt/data1/projects/UCB/01.Star

for f in $(ls)
do
echo $f
Rscript --vanilla /mnt/data1/projects/UCB/results_ucb/tangchao/IR_v2/MMCLGenomeRegion.R /mnt/data5/BGI/UCB/tangchao/IR/bedtools/intron/${f}_introncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/leftexon/${f}_leftexoncov /mnt/data5/BGI/UCB/tangchao/IR/bedtools/rightexon/${f}_rightexoncov /mnt/data5/BGI/UCB/tangchao/IR/Result/${f}.txt
done


















