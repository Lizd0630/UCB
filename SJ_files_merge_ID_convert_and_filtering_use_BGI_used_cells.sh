## SJ.out.tab files prepare and ID convert

## UCB1 -> UCB3
cd /mnt/data5/BGI/UCB/tangchao/data/SJ/UCB3

# copy old data to the directory
cp /mnt/data1/projects/UCB/data/STAR/SJ/UCB1* ./

# copy new data to the directory
for f in `ls /mnt/data5/BGI/UCB/Star/ | grep UCB1`
do cp /mnt/data5/BGI/UCB/Star/$f/*SJ.out.tab ./
done

# convert the name 
for f in `ls `
do mv $f UCB3${f:4:22}
done


## UCB5 -> UCB1
cd /mnt/data5/BGI/UCB/tangchao/data/SJ/UCB1

# copy old data to the directory
cp /mnt/data1/projects/UCB/data/STAR/SJ/UCB5* ./

# convert the name UCB1 -> ucb3
for f in `ls `
do mv $f UCB1${f:4:22}
done


## UCB3 -> UCB4
cd /mnt/data5/BGI/UCB/tangchao/data/SJ/UCB4

# copy old data to the directory
cp /mnt/data1/projects/UCB/data/STAR/SJ/UCB3* ./

# copy new data to the directory
for f in `ls /mnt/data5/BGI/UCB/Star/ | grep UCB3`
do cp /mnt/data5/BGI/UCB/Star/$f/*SJ.out.tab ./
done

# convert the name 
for f in `ls `
do mv $f UCB4${f:4:22}
done



## only use the cells that BGI used
mkdir /mnt/data5/BGI/UCB/tangchao/data/SJ/UCB_134
for f in `awk -F "\t" '{print $1}' /mnt/data5/BGI/UCB/ExpMat_NewID/Basic_stat.txt | grep UCB1`
do cp /mnt/data5/BGI/UCB/tangchao/data/SJ/UCB1/$f* /mnt/data5/BGI/UCB/tangchao/data/SJ/UCB_134
done


for f in `awk -F "\t" '{print $1}' /mnt/data5/BGI/UCB/ExpMat_NewID/Basic_stat.txt | grep UCB3`
do cp /mnt/data5/BGI/UCB/tangchao/data/SJ/UCB3/$f* /mnt/data5/BGI/UCB/tangchao/data/SJ/UCB_134
done


for f in `awk -F "\t" '{print $1}' /mnt/data5/BGI/UCB/ExpMat_NewID/Basic_stat.txt | grep UCB4`
do cp /mnt/data5/BGI/UCB/tangchao/data/SJ/UCB4/$f* /mnt/data5/BGI/UCB/tangchao/data/SJ/UCB_134
done