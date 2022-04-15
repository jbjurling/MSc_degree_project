# PRS for breast and ovarian cancer
Aim for the project was to construct polygenic risk scores constituted of polygenic factors identified in GWAS for breast cancer (BC) and ovarian cancer (OC) separately. The PRS was created in the UK Biobank (UKB) cohort. To estimate the weights for the PRS, GWAS summary statistics from the Breast Cancer Association Consortium (BCAC) and Ovarian Cancer Association Consortium (OCAC) were used. 

PRS was computed using LDpred2:

Privé F, Arbel J, Vilhjálmsson BJ. LDpred2: better, faster, stronger. Bioinformatics. 2020;36.

# Prepare GWAS summary statistics files and fam files
Downloaded from:

Breast cancer: BCAC - https://bcac.ccge.medschl.cam.ac.uk/

Ovarian cancer: OCAC - https://www.ebi.ac.uk/gwas/studies/GCST004417

The downloaded summary statistics files were rearranged to be used as input for LDpred2. Column order: temporary variant name (var_name), chromosome number (chr), genetic position (pos), reference allele (a0), alternate allele (a1), effect size estimate (beta), standard error of beta (beta_se), p-value for effect estimate (p). The temporary variant name (var_name) is the location of the SNP, so which chromosome it is present on and at which position.


```
#BREAST CANCER
cat icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt | tr " " "\t" | awk '{print $1,$1,$42,$44,$46,$4}' | tr "_" "\t" | awk '{print $1"_"$2,$5,$6,$7,$8,$9,$10,$11,$12}' | grep -v "var_name" | tr " " "," > BC_gwas_summary.txt

#OVARIAN CANCER
unzip Phelan_Archive.zip
rm Summary_chr23.txt # remove chrX (coded as 23) since we do not have chrX from BCAC
for i in Summary_chr*;do cat $i | grep -v "OrigSNPname" | tr "," "\t" | awk '{print $3"_"$4,$3,$4,$6,$5,$9,$10,$12,$7}'| tr " " ",";done > OC_gwas_sum.txt
cat OC_gwas_sum.txt | grep -v "\-99" > OC_gwas_summary.txt #to remove variants with unknown effect size and p-value

#Header added using nano
```
Next the summary files were filtered to only include variants that can also be found in the UKB cohort. Variants with a maf <0.01 were also filtered out for the OC GWAS summary. The downloaded BC GWAS summary statistics was already filtered for maf. 

First created a file containing rsID from the UKB SNPs, also added a var_name column to match agains the GWAS summary statistics. The "dedup.bim" files contain only unique rsID from UKB. 

```
cat *dedup.bim | cut -f 1-2,4 | sort -k1,1V -k3,3n | awk '{print $1"_"$3,$2}' > /proj/sens2017538/nobackup/Exjobb/Josefin/PRS/gwas_sum/UKBB_snp.txt
```

```
module load bioinfo-tools R_packages
R

library(dplyr)
setwd("/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/gwas_sum")

#BREAST CANCER
BC_gwas_sum<-read.table("breast_cancer/BC_gwas_summary.txt",header=T, stringsAsFactors=F, sep=",")
UKB_snps<-read.table("UKBB_snp.txt",header=F)
BC_snps<-which(BC_gwas_sum$var_name %in% UKB_snps$V1)
BC_gwas<-BC_gwas_sum[BC_snps,]
write.table(BC_gwas, file="breast_cancer/BC_gwas_summary_UKB.txt", col.names=T, row.names=F, quote=F, sep=",")

#OVARIAN CANCER
OC_gwas_sum<-read.table("ovarian_cancer/OC_gwas_summary.txt",header=T, stringsAsFactors=F, sep=",")
OC_snps<-which(OC_gwas_sum$var_name %in% UKB_snps$V1)
OC_gwas<-OC_gwas_sum[OC_snps,]
OC_final<-subset(OC_gwas, maf>=0.01)
write.table(OC_final, file="ovarian_cancer/OC_gwas_summary_maf001.txt", col.names=T, row.names=F, quote=F, sep=",")
```

".fam" files for the UKB participants to include in this study was created that contained case/control status for BC and OC. Disease status was extracted from the UKB phenotype files using ICD10, ICD9 and self-reported illness codes.

Breast cancer: ICD10 (C50 - Malignant neoplasm of breast), ICD9 (174 - Malignant neoplasm of female breast), self-reported illness (1002 - Breast cancer).

Ovarian cancer: ICD10 (C56 - Malignant neoplasm of ovary), ICD9 (183 - Malignant neoplasm of ovary and other uterine adnexa), self-reported illness (1039 - Ovarian cancer).

```
module load bioinfo-tools R_packages
R

library(dplyr)

#Load fam file for 41143 project and phenotype information from UKB
fam <- read.table("/proj/sens2017538/nobackup/UKBB_IMP_DOSAGE_V3_bim_bed/ukb41143.fam", header=F, stringsAsFactors=F)
load("/proj/sens2017538/nobackup/UKBB_41143_Data/All.Phenos.20210212.RData")
id<-phenos$f.eid

#BREAST CANCER
source("/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/breast_cancer/GetDiseaseUKB.plinkformat.R")
ICD10<-c("C50") 
ICD9<-c("174") 
SelfRep<-c(1002) 
BC<-GetDiseasesUKB(phenos, ICD10, ICD9, SelfRep)
bc<-data.frame(id,BC)
full_bc<- left_join(fam,bc, by=c("V1"="id"))
full_out<-subset(full_bc,select=-c(V6))
write.table(full_out,file="/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/breast_cancer/ukb41143.BC.tmp.fam",col.names=F,row.names=F)

#OVARIAN CANCER
source("/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/breast_cancer/GetDiseaseUKB.plinkformat.R")
ICD10<-c("C56") 
ICD9<-c("183")
SelfRep<-c(1039)
OC<-GetDiseasesUKB(phenos, ICD10, ICD9, SelfRep)
oc<-data.frame(id,OC)
full_oc<- left_join(fam,oc, by=c("V1"="id"))
oc_out<-subset(full_oc,select=-c(V6))
write.table(oc_out,file="/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/ovarian_cancer/ukb41143.OC.tmp.fam",col.names=F,row.names=F)

```

Individuals missing with missing disease information is coded as NA after joining the two tables together. Changed NA values to -9 to fit plink fam file format.

```
cat ukb41143.BC.tmp.fam | sed 's/NA/-9/´g > ukb41143.BC.fam

cat ukb41143.OC.tmp.fam | sed 's/NA/-9/´g > ukb41143.OC.fam

```

# Running LDpred2
Plink files from UKB is separated by chromosome and therefore LDpred2 was computed for each chromosome separately. The combined results were then 


