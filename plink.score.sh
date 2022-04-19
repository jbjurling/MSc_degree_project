#! /bin/bash -l
#SBATCH -A sens2017538
#SBATCH -p core -n 10
#SBATCH -t 48:00:00

module load bioinfo-tools
module load plink

path=/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/breast_cancer #/ovarian_cancer for OC
genotypes=/proj/sens2017538/nobackup/UKBB_IMP_DOSAGE_V3_bim_bed
out=/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/breast_cancer/PRS_score #/ovarian_cancer/PRS_score for OC

cd $out

for chr in {1..22}

do

plink --bed $genotypes/chr${chr}.bed \
--bim $genotypes/chr${chr}.dedup.bim \
--fam $path/ukb41143.BC.fam \ #/ukb41143.OC.fam for OC
--score $path/beta_grid.score 1 2 ${score} header sum \
--out PRS_score_${chr}_${score}

done
