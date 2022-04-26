#!/bin/bash -l
#
# Script to filter out variants in genes from genepanel

while read -r line
do 
chr=$line

path=/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/vep_files

echo "#! /bin/bash -l" > filterVEP.sc;
echo "#SBATCH -A sens2017538" >> filterVEP.sc;
echo "#SBATCH -p core -n 1" >> filterVEP.sc;
echo "#SBATCH -J $chr.filter" >> filterVEP.sc;
echo "#SBATCH -t 15:00" >> filterVEP.sc;
echo "#SBATCH -e /proj/sens2017538/nobackup/Exjobb/Josefin/ErrorAndOut/$chr.filter.error" >> filterVEP.sc;

echo "cd $path || exit 1" >> filterVEP.sc;

# Filter out genes associated with breast cancer
echo "grep \"^#\" UKBBexome_$chr.vep > head_$chr" >> filterVEP.sc;
echo "cat genepanel_BC.txt | while read line; do cat UKBBexome_$chr.vep | grep \$line;done > bc_$chr" >> filterVEP.sc;
echo "cat bc_$chr | grep \"HIGH\" > bc_h_$chr" >> filterVEP.sc;
echo "cat head_$chr bc_h_$chr > UKBBexome_${chr}_filtered_BC.vep" >> filterVEP.sc;

#Filter out genes associated with ovarian cancer
echo "cat genepanel_OC.txt | while read line; do cat UKBBexome_$chr.vep | grep \$line;done > oc_$chr" >> filterVEP.sc;
echo "cat oc_$chr | grep \"HIGH\" > oc_h_$chr" >> filterVEP.sc;
echo "cat head_$chr oc_h_$chr > UKBBexome_${chr}_filtered_OC.vep" >> filterVEP.sc;

echo "rm head_$chr bc_$chr bc_h_$chr oc_$chr oc_h_$chr" >> filterVEP.sc;

sbatch filterVEP.sc;

done <chr_v2.txt
