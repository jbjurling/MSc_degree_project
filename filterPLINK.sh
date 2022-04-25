#!/bin/bash -l
#
# Script to filter out genotypes for variants in genes from genepanel

while read line1 line2
do 
first=$line1
second=$line2

path=/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/plink_files
plink=/proj/sens2017538/nobackup/ExomeSeq2ndRel
txtfile=/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/vep_files

echo "#! /bin/bash -l" > filterPLINK.sc;
echo "#SBATCH -A sens2017538" >> filterPLINK.sc;
echo "#SBATCH -p core -n 1" >> filterPLINK.sc;
echo "#SBATCH -J $second.filterPLINK" >> filterPLINK.sc;
echo "#SBATCH -t 5:00" >> filterPLINK.sc;
echo "#SBATCH -e /proj/sens2017538/nobackup/Exjobb/Josefin/ErrorAndOut/$second.filter.error" >> filterPLINK.sc;

echo "module load bioinfo-tools" >> filterPLINK.sc;
echo "module load plink2" >> filterPLINK.sc;

echo "cd $path" >> filterPLINK.sc;

echo "plink2 --bed $plink/ukb23155_${first}_b0_v1.bed --bim $plink/UKBexomeOQFE_$second.bim --fam $plink/ukb41143.fam --extract $txtfile/BC_variants.vep --keep-females --freq --make-bed --out BC_UKBexome_$second" >> filterPLINK.sc;
echo "plink2 --bed $plink/ukb23155_${first}_b0_v1.bed --bim $plink/UKBexomeOQFE_$second.bim --fam $plink/ukb41143.fam --extract $txtfile/OC_variants.vep --keep-females --freq --make-bed --out OC_UKBexome_$second" >> filterPLINK.sc;

sbatch filterPLINK.sc;

done <chr.txt
