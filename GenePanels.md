# Investigating gene panels for breast and ovarian cancer
To investigate potential disease causing variants in genes associated with breast and ovarian cancer, LoF variants were identified in the UK biobank (UKB) whole exome sequencing (WES) data. Variants within 32 genes for BC (BC_genes.txt) and 25 genes for OC (OC_genes.txt) were annotated. The number of LoF variants were then compared between individuals with BC/OC and controls to to check how good the current genetic testing is in detecting females with disease.

For annotation the Ensembl Variant Effect Predictor version 99 (VEP/99) was used:

McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GRS, Thormann A, et al. The Ensembl Variant Effect Predictor. Genome Biol. 2016 Dec;17(1):122.


# Annotation of WES data
The plink files containing WES data for 200K individuals from UKB were converted to VCF files to be used as input for VEP annotation. Genotype information from the WES data was filtered out in the same step (done to save space in the project directory, not necessary to do). The script plink2vcf.sh was used for this. The variants were then annotated using VEP/99 (see VEP.sh).

The VEP output files were filtered to only include variants found in genes that have known association with breast and ovarian cancer (filterVEP.sh). The filtered files were then used as input for filtering the WES plink files from UKB to only contain genotypes for these variants (filterPLINK.sh).

To get one combined file with all variants, all chromosome files were concatenated. 

```
#BREAST CANCER
ls BC* | grep -v "afreq" | grep -v "log" | tr "." "\t" | cut -f 1 | uniq > BC_file.txt
plink --merge-list BC_file.txt --freq --export A --out BC_variants

#OVARIAN CANCER
ls OC* | grep -v "afreq" | grep -v "log" | tr "." "\t" | cut -f 1 | uniq > OC_file.txt
plink --merge-list OC_file.txt --freq --export A --out OC_variants
```

# Extracting phenotype information from UKB
To determine the case/control status (to use for statistical testing), disease status was extracted from phenotype files from UKB. This was done using case_control_status.R 

Scripts to extract disease information (GetDiseasesUKB.R) and age of onset (GetDiseasesUKB_ages.R) was provided by [@AasaJohanssonUU](https://github.com/AasaJohanssonUU) 

# Statistical testing
To explore if there are significant differences between cases and controls ran two-sample t-test using the t.test function in R/4.1.1. This was done for 1) all variants in the gene panel, 2) variants in gene panel from the Swedish National Board of Healt and Welfare, 3) variants per gene from the full gene panel and 4) each variants separately. 

1) T-test for full set of variants:
```
module load bioinfo-tools R_packages
R

setwd("/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/R_analysis")

#Load genotype data
BC_full <- read.table("BC_variants_status.RData") #OC_variants_status.RData for OC

#Added column with case/control status and removed first 11 columns to get table to use for two sample t-test.
BC_full$status <- rep("control", nrow(BC_full))
BC_cases <- rownames(BC_full[BC_full$BC == 1,])
BC_full[BC_cases,]$status="case"
BC_ttest <- BC_full[-c(1:11)]
 
#Get sum of variants per individual
BC_case_table <- subset(BC_ttest, status == "case")
BC_case_table$ind_sum <- rowSums(BC_case_table == "2" | BC_case_table == "1", na.rm = T)
BC_control_table <- subset(BC_ttest, status=="control")
BC_control_table$ind_sum <- rowSums(BC_case_table == "2" | BC_case_table == "1", na.rm = T)

#t.test between cases and controls (sum of variants per individual):
case <- BC_case_table$ind_sum
control <- BC_control_table$ind_sum
ttest <- t.test(case, control)
```

2) T-test for gene panel from the National Board of Health and Welfare:
```
module load bioinfo-tools R_packages
R

setwd("/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/R_analysis/ttest")
list <- read.table("/proj/sens2017538/nobackup/Exjobb/Josefin/Scripts/BC_SOCgenes.txt",header=F)
BC_table <- read.table("/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/R_analysis/BC_variants_status.RData", header = T, stringsAsFactors = F, sep = "\t")
BC_status <- data.frame(BC_table$id, BC_table$status)

for (name in list[,1]){
gene <- read.table(/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/plink_files/genes/paste0(name,"_BC.raw"), header = T, stringsAsFactors = F, sep="")
gene_table <- merge(gene, BC_status, by.x = "FID", by.y = "BC_table.id")
gene_for_test <- gene_table[-c(1:6)]
gene_case <- subset(gene_for_test, BC_table.status == "case")
gene_control <- subset(gene_for_test, BC_table.status == "control")
gene_case$sum_variants <- rowSums(gene_case == "2" | gene_case == "1", na.rm=T)
gene_control$sum_variants <- rowSums(gene_control == "2" | gene_control == "1", na.rm=T)
case <- gene_case$sum_variants
control <- gene_control$sum_variants
ttest <- t.test(case, control)
sink(paste0(name,"_SOC_ttest.log"))
print(ttest)
sink()
}

q()

cat *SOC_ttest.log > socialstyrelsen_BC_genes.txt
rm *SOC_ttest.log
```

Extracted genotype data per gene using plink_genes.sh and then ran t-test per gene to see if specific genes in the gene panels gave the significant differences between cases and controls.

3) T-test for each gene separately.

```
module load bioinfo-tools R_packages
R

setwd("/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/R_analysis/ttest")
list <- read.table("/proj/sens2017538/nobackup/Exjobb/Josefin/Scripts/BC_genes.txt", header = F) #OC_genes for OC analysis
BC_table <- read.table("/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/R_analysis/BC_variants_status.RData", header = T, stringsAsFactors = F, sep = "\t") #OC_variants_status.RData for OC analysis
BC_status <- data.frame(BC_table$id, BC_table$status)

for (name in list[,1]){
gene <- read.table(/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/plink_files/genes/paste0(name,"_BC.raw"), header = T, stringsAsFactors = F, sep = "")
gene_table <- merge(gene, BC_status, by.x = "FID", by.y = "BC_table.id")
gene_for_test <- gene_table[-c(1:6)]
gene_case <- subset(gene_for_test, BC_table.status == "case")
gene_control <- subset(gene_for_test, BC_table.status == "control")
gene_case$sum_variants <- rowSums(gene_case == "2" | gene_case == "1", na.rm = T)
gene_control$sum_variants <- rowSums(gene_control == "2" | gene_control == "1", na.rm = T)
case<-gene_case$sum_variants
control <- gene_control$sum_variants
ttest <- t.test(case, control)
sink(paste0(name,"_BC_ttest.log"))
print(ttest)
sink()
}
```

4) T-test per variant. 

In this step created a temporary directory with output files from the loop. These were then concatenated and moved to the same directory as the rest of the t-test statistics. The temporary files were removed at the end of the test.

```
cd /proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/R_analysis/
mkdir tmp

module load bioinfo-tools R_packages
R

setwd("/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/R_analysis/tmp")
list <- read.table("/proj/sens2017538/nobackup/Exjobb/Josefin/Scripts/BC_variants.txt", header = F)

for (name in list[,1]){
BC_variant <- data.frame(BC_table$FID, BC_table[,name])
variant_table <- merge(BC_variant, BC_status, by.x = "BC_table.FID", by.y = "BC_table.id")
variant_case <- subset(variant_table, BC_table.status == "case")
variant_control <- subset(variant_table, BC_table.status == "control")
variant_case$sum_variants <- rowSums(variant_case == "2" | variant_case == "1", na.rm = T)
variant_control$sum_variants <- rowSums(variant_control == "2" | variant_control == "1", na.rm = T)
case <- variant_case$sum_variants
control <- variant_control$sum_variants
name <- c(name)
ttest <- data.frame(name, t.test(case, control)$p.value)
sink(paste0(name,"_BC.txt"))
print(ttest)
sink()
}

q()

cat *_BC.txt > /proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/R_analysis/ttest/BC_variants_pvalue.txt
rm -r tmp
```

