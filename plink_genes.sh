#!/bin/bash -l
#
# Script to get raw data per gene for ttest in R. 

while read line
do
name=$line

path=/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/plink_files
genes=/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/R_analysis/genes

echo "#!/bin/bash -l" > plink_gene.sc;
echo "#SBATCH -A sens2017538" >> plink_gene.sc;
echo "#SBATCH -p core -n 1" >> plink_gene.sc;
echo "#SBATCH -J $name.plink_gene" >> plink_gene.sc;
echo "#SBATCH -t 10:00" >> plink_gene.sc;
echo "#SBATCH -o /proj/sens2017538/nobackup/Exjobb/Josefin/ErrorAndOut/$name.plink_gene.out" >> plink_gene.sc;
echo "#SBATCH -e /proj/sens2017538/nobackup/Exjobb/Josefin/ErrorAndOut/$name.plink_gene.error" >> plink_gene.sc;

echo "module load bioinfo-tools" >> plink_gene.sc;
echo "module load plink" >> plink_gene.sc;

echo "cd $genes" >> plink_gene.sc;
echo "plink --bfile $path/BC_variants --extract $name.txt --export A --out ${name}_BC" >> plink_gene.sc;

sbatch plink_gene.sc;

done <BC_genes.txt
