# VEP.sh

#!/bin/bash -l
#
# Script to annotate variants from exome data using VEP. Run on Rackham, internet access needed.

while read line
do
chr=$line

path=/proj/snic2021-22-624/nobackup/temp/Josefin/vcf_files
output=/proj/snic2021-22-624/nobackup/temp/Josefin/vep_files

echo "#! /bin/bash -l" > VEP.sc;
echo "#SBATCH -A snic2021-22-624" >> VEP.sc;
echo "#SBATCH -p core -n 1" >> VEP.sc;
echo "#SBATCH -J $chr.VEP" >> VEP.sc;
echo "#SBATCH -t 6:00:00" >> VEP.sc; #Took ~5 h for chr1 and 2 to run, 20 min for chr21
echo "#SBATCH -o /proj/snic2021-22-624/nobackup/temp/Josefin/Error/$chr.VEP.out" >> VEP.sc;
echo "#SBATCH -e /proj/snic2021-22-624/nobackup/temp/Josefin/Error/$chr.VEP.error" >> VEP.sc;

#Load modules
echo "module load bioinfo-tools" >> VEP.sc;
echo "module load vep" >> VEP.sc;

echo "cd $output" >> VEP.sc;

echo "vep --cache --dir /sw/data/uppnex/vep/99/ -i $path/UKBBexome_$chr.vcf -o UKBBexome_$chr.vep" >> VEP.sc;

sbatch VEP.sc;

done <chr.txt
