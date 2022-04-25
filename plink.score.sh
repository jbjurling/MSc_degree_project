#! /bin/bash -l 
#SBATCH -A sens2017538
#SBATCH -p core -n 8
#SBATCH -t 48:00:00

genotype=/proj/sens2017538/nobackup/UKBB_IMP_DOSAGE_V3_bim_bed
path=/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/breast_cancer #ovarian_cancer for OC

module load bioinfo-tools
module load plink2

cd /proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/breast_cancer/PRS_score #ovarian_cancer for OC

for chr in {1..22}

do

plink2 --pgen $genotype/chr${chr}.pgen \
--psam $path/ukb41143.BC.fam \ #.OC. for OC
--pvar $genotype/chr${chr}.dedup.bim \
--score $path/beta_grid.score.a1 1 2 ${score} header cols=+scoresums \
--out model_${score}_chr${chr}

done
