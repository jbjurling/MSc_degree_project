# PRS for breast and ovarian cancer
Aim for the project was to construct polygenic risk scores constituted of polygenic factors identified in GWAS for breast cancer (BC) and ovarian cancer (OC) separately. The PRS was created in the UK Biobank (UKB) cohort. To estimate the weights for the PRS, GWAS summary statistics from the Breast Cancer Association Consortium (BCAC) and Ovarian Cancer Association Consortium (OCAC) were used. 

PRS was computed using LDpred2:

Privé F, Arbel J, Vilhjálmsson BJ. LDpred2: better, faster, stronger. Bioinformatics. 2020;36.

# Prepare GWAS summary statistics files and fam files
Downloaded from:

Breast cancer: BCAC - https://bcac.ccge.medschl.cam.ac.uk/

Ovarian cancer: OCAC - https://www.ebi.ac.uk/gwas/studies/GCST004417

The downloaded summary statistics files were rearranged to be used as input for LDpred2. Column order: temporary variant name (var_name), chromosome number (chr), genetic position (pos), non-effect allele (a0), effect allele (a1), effect size estimate (beta), standard error of beta (beta_se), p-value for effect estimate (p), and effect allele frequency (here named maf). The temporary variant name (var_name) is the location of the SNP, so which chromosome it is present on and at which position.


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
The summary files were filtered to only include variants that can also be found in the UKB cohort. In addition, variants with a maf <0.01 were filtered out for the OC GWAS summary. (The downloaded BC GWAS summary statistics was already filtered for maf <0.01).

Created a file containing rsID from the UKB SNPs, also added a var_name column to match agains the GWAS summary statistics. The "dedup.bim" files contain only unique rsID from UKB. 

```
cat *dedup.bim | cut -f 1-2,4 | sort -k1,1V -k3,3n | awk '{print $1"_"$3,$2}' > /proj/sens2017538/nobackup/Exjobb/Josefin/PRS/gwas_sum/UKBB_snp.txt
```

```
module load bioinfo-tools R_packages
R

library(dplyr)
setwd("/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/gwas_sum")

#BREAST CANCER
BC_gwas_sum <- read.table("breast_cancer/BC_gwas_summary.txt", header = T, stringsAsFactors = F, sep = ",")
UKB_snps <- read.table("UKBB_snp.txt", header = F)
BC_snps <- which(BC_gwas_sum$var_name %in% UKB_snps$V1)
BC_gwas <- BC_gwas_sum[BC_snps,]
write.table(BC_gwas, file="breast_cancer/BC_gwas_summary_UKB.txt", col.names = T, row.names = F, quote = F, sep = ",")

#OVARIAN CANCER
OC_gwas_sum <- read.table("ovarian_cancer/OC_gwas_summary.txt", header = T, stringsAsFactors = F, sep = ",")
OC_snps <- which(OC_gwas_sum$var_name %in% UKB_snps$V1)
OC_gwas <- OC_gwas_sum[OC_snps,]
OC_final <- subset(OC_gwas, maf >= 0.01)
write.table(OC_final, file="ovarian_cancer/OC_gwas_summary_maf001.txt", col.names = T, row.names = F, quote = F, sep = ",")
```

".fam" files for the UKB participants to include in this study was created that contained case/control status for BC and OC. Disease status was extracted from the UKB phenotype file using ICD10, ICD9 and self-reported illness codes.

Breast cancer: ICD10 (C50 - Malignant neoplasm of breast), ICD9 (174 - Malignant neoplasm of female breast), self-reported illness (1002 - Breast cancer).

Ovarian cancer: ICD10 (C56 - Malignant neoplasm of ovary), ICD9 (183 - Malignant neoplasm of ovary and other uterine adnexa), self-reported illness (1039 - Ovarian cancer).

```
module load bioinfo-tools R_packages
R

library(dplyr)

#Load fam file for 41143 project and phenotype information from UKB
fam <- read.table("/proj/sens2017538/nobackup/UKBB_IMP_DOSAGE_V3_bim_bed/ukb41143.fam", header = F, stringsAsFactors = F)
load("/proj/sens2017538/nobackup/UKBB_41143_Data/All.Phenos.20210212.RData")
id <- phenos$f.eid

#BREAST CANCER
source("/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/breast_cancer/GetDiseaseUKB.plinkformat.R")
ICD10 <- c("C50") 
ICD9 <- c("174") 
SelfRep <- c(1002) 
BC <- GetDiseasesUKB(phenos, ICD10, ICD9, SelfRep)
bc_df <- data.frame(id,BC)
full_bc <- left_join(fam,bc_df, by = c("V1"="id"))
full_out <- subset(full_bc, select = -c(V6))
write.table(full_out,file="/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/breast_cancer/ukb41143.BC.tmp.fam",col.names=F,row.names=F)

#OVARIAN CANCER
source("/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/breast_cancer/GetDiseaseUKB.plinkformat.R")
ICD10 <- c("C56") 
ICD9 <- c("183")
SelfRep <- c(1039)
OC <- GetDiseasesUKB(phenos, ICD10, ICD9, SelfRep)
oc_df <- data.frame(id,OC)
full_oc <- left_join(fam,oc_df, by = c("V1"="id"))
oc_out <- subset(full_oc, select = -c(V6))
write.table(oc_out,file="/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/ovarian_cancer/ukb41143.OC.tmp.fam", col.names = F, row.names = F)

```

Individuals with missing disease information is coded as NA after joining the two tables together. Changed NA values to -9 to fit plink fam file format.

```
cat ukb41143.BC.tmp.fam | sed 's/NA/-9/´g > ukb41143.BC.fam

cat ukb41143.OC.tmp.fam | sed 's/NA/-9/´g > ukb41143.OC.fam

```

# Running LDpred2
Plink files from UKB is separated by chromosome and therefore LDpred2 was computed for each chromosome separately. The combined results were then combined before evaluating the best-fit PRS model. 

The dataset, consisting of 263,313 females from the UKB, were split into three cohorts. The validation and testing cohort is used in this analysis. The validation cohort contain approximately 10% of the individuals from the dataset. An additional cohort, intended to be used as LD reference, was also created but due to the time-limit of this project LD references from HapMap3 were used instead.

``` 
module load bioinfo-tools R_packages
R

library(bigsnpr)

# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach("/proj/sens2017538/nobackup/UKBB_IMP_DOSAGE_V3_bim_bed/chr1.rds")
# Get aliases for useful slots
G <- obj.bigSNP$genotypes

#Split cohort (for correlation and ld calculations want the Caucasian non related individuals since we have Caucasian GWAS data.)
setwd("/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/")
load("/proj/sens2017538/nobackup/UKBB_41143_Data/All.Phenos.20210212.IDs_to_indlude_NonRelCauc.RData")
ind.females <- which(obj.bigSNP$fam$sample.ID %in% temp$f.eid[temp$f.31.0.0==0])
ind.sample <- sample(ind.females, 1000) 
cohort <- setdiff(rows_along(G), ind.sample)
ind.val <- sample(cohort, 26430) 
ind.test <- setdiff(cohort, ind.val)

save(ind.sample, file = "ldscore_cohort.RData")
save(ind.val, file = "validation_cohort.RData")
save(ind.test, file = "test_cohort.RData")
```

Plink files containing genotype information from UKB was loaded to R and saved as ".rds" files. Filtering was made so the ".rds" files only contained genotypes from the participants for this study (487,409 individuals) and SNPs present in the GWAS summary statistics files from BCAC and OCAC. The script makeRDS.R was used for this step. 

HapMap3 containing 1,054,330 variants based on 362,320 European individuals of the UK biobank was used to get SNP and LD matrices (see run.ldref.R). Ldpred2-grid was then run using the script beta_grid.R.

The genetic load for each individual was calculated using the --score function in plink2. This was done for all PRS models produced by LDpred2-grid. The output from beta_grid.R was concatenated with rsID for each variant and the alternate allele for the variant to use this as input for the allelic scoring.

```
module load bioinfo-tools R_packages
R

library(bigsnpr)
library(dplyr)

setwd("/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/ovarian_cancer/") #/ovarian_cancer/ for OC 
load("working.environment.ldref.RData")
load("beta_grid.RData")

out <- data.frame(df_beta$rsid, df_beta$a1, beta_grid)
write.table(out, file = "beta_grid.score", col.names = T, row.names = F , quote = F, sep = " ")
```
The allelic scoring with plink was run using the script plink.score.sh

Since it took between 20-40h to run each script, they were run in parallell by looping over the columns including the models.

```
for score in {3..170}
do

sbatch --job-name=$score.plink --output=$score.out --export=score=$score plink.score.sh

done
```

The output files from allelic scoring with plink was merged to get the sum of scores across all chromosomes for each PRS model. Log files from plink2 were moved to a separate folder before merging the scores. Function "get_folder" from [@Schmytzi](https://github.com/Schmytzi/ldpred-prs-manual/blob/master/merge_scores.R)

```
module load bioinfo-tools R_packages
R

setwd("/castor/project/proj_nobackup/Exjobb/Josefin/PRS/LDpred2/breast_cancer") #/ovarian_cancer for OC

library(magrittr)
library(data.table)
library(dplyr)

get_folder <- function(folder) {
	data <- list.files(path=folder, pattern=paste0("model_", number), full.names=T) %>%
        lapply(fread) %>%
        rbindlist()
	scores <- data[, .(score = sum (SCORE1_SUM)), by=IID]
	scores
}

#Loop over all models. 
	#In same step standardise the score from each model and concatenate them all into one table.
	
output<-{}

for (number in 3:170){
score <- get_folder("/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/breast_cancer/PRS_score") #/ovarian_cancer/PRS_score for OC
if(number==3){colnames(score) <- c("IID", paste0("model_", number))
score_st <- score %>% mutate_at(c(paste0("model_", number)), ~(scale(.) %>% as.vector)) #Standardize score
output<-cbind(output,score_st)}
else{score_sub <- subset(score, select=-c(IID))
colnames(score_sub) <- c(paste0("model_", number))
score_st <- score_sub %>% mutate_at(c(paste0("model_", number)), ~(scale(.) %>% as.vector))
output <- cbind(output,score_st)}
}

save(output, file="PRS_score_table.RData")
  
```


# Validation of PRS models
To find the best-fit PRS model used Z-score from logistic regression, AUC and how many cases the 10th decile contained. All analysis were carried out using the validation cohort. 

The script glm_prs.R was used to look at the Z-score for all PRS models. The results were then plotted using ggplot2 to look at the differences between using different: 
- h2 estimates
- proportions of variants that are believed to have an effect

```
module load bioinfo-tools R_packages
R

library(ggplot2)
setwd("/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/breast_cancer")
load("params.RData")

plot <- ggplot(params, aes(x = p, y = score, color = as.factor(h2))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
  facet_wrap(~ sparse, labeller = label_both) +
  labs(y = "GLM Z-Score", color = "h2") +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))

```

AUC calculations:

```
library(pROC)

AUC<-{}

#Change to OC for OC analysis
for (number in 1:168){
 model<-glm(BC[ind.val,] ~ PRS_score[,{number}][ind.val], family="binomial")
 predicted <- predict(model, PRS_score[ind.val,], type="response")
 auc<-auc(BC[ind.val,],predicted)
 df<-data.frame(auc)
 rownames(df)<-c(paste0("model_", {number}))
 AUC<-rbind(AUC,df)
}

AUC_cov<-{}

#With covariates
for (number in 1:168){
 model<-glm(BC[ind.val,] ~ PRS_score[,{number}][ind.val]+yob[ind.val,]+age[ind.val,]+BMI[ind.val,], family="binomial")
 predicted <- predict(model, PRS_score[ind.val,], type="response")
 auc_cov<-auc(BC[ind.val,],predicted)
 df<-data.frame(auc_cov)
 rownames(df)<-c(paste0("model_", {number}))
 AUC_cov<-rbind(AUC_cov,df)
}

AUC_PRS<-cbind(AUC,AUC_cov)

save(AUC_PRS, file="AUC_PRS_models.RData")
```

Proportion of cases in the 10th decentile of the validation set was calculated. To see how variables were loaded check glm_prs.R

```
tmp <- inner_join(output, fam, by = c("IID" = "V2"))
data_dec <- subset(tmp, select = -c(IID, V3, V4, V5))

data_dec$V6[data_dec$V6 == "1"] <- 0
data_dec$V6[data_dec$V6 == "2"] <- 1
data_dec$V6[data_dec$V6 == "-9"] <- NA

top10 <- {}

for (model in 1:168){
  scoring <- data.frame(data_dec[,169], data_dec[,170], data_dec[,..model])
  colnames(scoring) <- c("id", "BC", paste0("model_", model))
  sub <- scoring[ind.val,]
  sub$decentile <- dplyr::ntile(sub[,3], 10)
  top.dec <- subset(sub, decentile == 10)
  BC_case <- NROW(which((top.dec$BC == 1) == TRUE)) / NROW(top.dec$BC)
  PRS = c(paste0("model_", model))
  out <- data.frame(PRS, BC_case)
  top10 <- rbind(top10, out)
}
```

# Evaluation of best-fit PRS
The best PRS model was evaluated using the test cohort from UKB. 

```
#load test dataset
load("/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/test_cohort.RData")

# run logistic regression
model_eval <- glm(y[ind.test] ~ PRS_score[,135][ind.test], family = "binomial")
model_eval_cov <- glm(y[ind.test] ~ PRS_score[,135][ind.test] + yob[ind.test,] + age[ind.test,] + BMI[ind.test,], family = "binomial")

# Check AUC for model
predicted <- predict(model_eval, PRS_score[ind.test,], type = "response")
predicted_cov <- predict(model_eval_cov, PRS_score[ind.test,], type = "response")

auc <- auc(BC[ind.test,], predicted)
auc_cov <- auc(BC[ind.test,], predicted_cov)

#OR 
OR <- odds.ratio(model_eval)
OR_cov <- odds.ratio(model_eval_cov)

#Proportion cases in 10th decentile
final_prs <- data.frame(data_dec[,169], data_dec[,170], data_dec[,135])
colnames(final_prs) <- c("id", "BC", "PRS_score")
subset.df <- final_prs[ind.test,]
subset.df$decentile <- dplyr::ntile(subset.df[,3], 10)
top.decentile <- subset(subset.df, decentile == 10)
case_proportion <- NROW(which((top.decentile$BC==1)==TRUE)) / NROW(top.decentile$BC)

#OR for 10th decentile 
ind.dec <- which(subset.df$decentile == 10)
model_top_dec <- glm(y[ind.dec] ~ PRS_score[,135][ind.dec] + yob[ind.dec,] + age[ind.dec,] + BMI[ind.dec,], family = "binomial")
OR_top10 <- odds.ratio(model_top_dec)

```
