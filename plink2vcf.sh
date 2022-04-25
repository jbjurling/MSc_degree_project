#!/bin/bash -l
#
# Script to create vcf files from plink .bim .bed .fam files
# First step for annotation with VEP

while read line1 line2
do
first=$line1
second=$line2

path=/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/vcf_files
plinkfiles=/proj/sens2017538/nobackup/ExomeSeq2ndRel

echo "#! /bin/bash -l" > plink2vcf.sc;
echo "#SBATCH -A sens2017538" >> plink2vcf.sc;
echo "#SBATCH -p core -n 1" >> plink2vcf.sc;
echo "#SBATCH -J $second.plink2vcf" >> plink2vcf.sc;
echo "#SBATCH -t 5:00:00" >> plink2vcf.sc;
echo "#SBATCH -o /proj/sens2017538/nobackup/Exjobb/Josefin/ErrorAndOut/$second.plink2vcf.out" >> plink2vcf.sc;
echo "#SBATCH -e /proj/sens2017538/nobackup/Exjobb/Josefin/ErrorAndOut/$second.plink2vcf.error" >> plink2vcf.sc;

#Load modules
echo "module load bioinfo-tools" >> plink2vcf.sc;
echo "module load plink2" >> plink2vcf.sc;

echo "cd $plinkfiles" >> plink2vcf.sc;

echo "plink2 --bed ukb23155_${first}_b0_v1.bed --bim UKBexomeOQFE_$second.bim --fam ukb41143.fam --recode vcf --freq --out $path/UKBBexome_$second" >> plink2vcf.sc;

echo "cd $path" >> plink2vcf.sc;
echo "cat UKBBexome_$second.vcf | awk -vOFS=\"\t\" '{print \$1,\$2,\$3,\$4,\$5}' > tmp_$second" >> plink2vcf.sc;
echo "mv tmp_$second UKBBexome_$second.vcf" >> plink2vcf.sc;

sbatch plink2vcf.sc;

done <chr.txt
