# Investigating gene panels for breast and ovarian cancer
To investigate potential disease causing variants in genes associated with breast and ovarian cancer, LoF variants were identified in the UKB WES data. Variants within 32 genes for BC (genepanel_BC.txt) and 25 genes for OC (genepanel_OC.txt) were annotated. The number of LoF variants were then compared between individuals with BC/OC and controls to to check the efficacy of the current genetical testing approaches. 

For annotation the Ensembl Variant Effect Predictor version 99 (VEP/99) was used:

McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GRS, Thormann A, et al. The Ensembl Variant Effect Predictor. Genome Biol. 2016 Dec;17(1):122.


# Annotation of WES data
The plink files containing WES data for 200K individuals from UKB were converted to VCF files to be used as input for VEP annotation. To save space in the project directory genotype information from the WES data was filtered out in the same step. The script plink2vcf.sh was used for this. The variants were then annotated using VEP/99 (see VEP.sh).

The VEP output files were filtered to only include variants found in genes that have known association with breast and ovarian cancer (filterVEP.sh). The filtered files were then used as input for filtering the WES plink files from UKB to only contain genotypes for these variants (filterPLINK.sh).

# Extracting phenotype information from UKB
To determine the case/control status for statistical testing, disease status was extracted from phenotype files from UKB. This was done using case_control_status.R 

Scripts to extract disease information and age of onset was provided by Åsa Johansson. 

# Statistical testing
To explore if there are significant differences between cases and controls ran two-sample t-test using the t.test function in R/4.1.1. This was done for all variants in the gene panel, variants per gene from the gene panel and each variants separately. 

T-test for full set of variants:
```
module load bioinfo-tools R_packages
R

setwd("/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/R_analysis")

#Load genotype data
BC_full<-read.table("BC_variants_status.RData") #OC_variants_status.RData for OC

#Added column with case/control status and removed first 11 columns to get table to use for two sample t-test.
BC_full$status <- rep("control",nrow(BC_full))
BC_cases<-rownames(BC_full[BC_full$BC==1,])
BC_full[BC_cases,]$status="case"
BC_ttest<-BC_full[-c(1:11)]
 
#Get sum of variants per individual
BC_case_table<-subset(BC_ttest, status=="case")
BC_case_table$ind_sum<-rowSums(BC_case_table=="2"|BC_case_table=="1",na.rm=T)
BC_control_table<-subset(BC_ttest, status=="control")
BC_control_table$ind_sum<-rowSums(BC_case_table=="2"|BC_case_table=="1",na.rm=T)

#t.test between cases and controls (sum of variants per individual):
case<-BC_case_table$ind_sum
control<-BC_control_table$ind_sum
ttest<-t.test(case,control)
```
Extractted genotype data per gene using plink_genes.sh and then ran t-test per gene to see if specific genes in the gene panels gave the significant differences between cases and controls.

```
module load bioinfo-tools R_packages
R

setwd("/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/R_analysis/ttest")
list<-read.table("/proj/sens2017538/nobackup/Exjobb/Josefin/Scripts/BC_genes.txt",header=F) #OC_genes for OC analysis
BC_table <- read.table("/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/R_analysis/BC_variants_status.RData", header=T,stringsAsFactors=F,sep="\t") #OC_variants_status.RData for OC analysis
BC_status<-data.frame(BC_table$id,BC_table$status)

for (name in list[,1]){
gene <- read.table(/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/plink_files/genes/paste0(name,"_BC.raw"), header = T, stringsAsFactors = F, sep="")
gene_table<-merge(gene,BC_status,by.x="FID",by.y="BC_table.id")
gene_for_test<-gene_table[-c(1:6)]
gene_case<-subset(gene_for_test,BC_table.status=="case")
gene_control<-subset(gene_for_test,BC_table.status=="control")
gene_case$sum_variants<-rowSums(gene_case=="2"|gene_case=="1",na.rm=T)
gene_control$sum_variants<-rowSums(gene_control=="2"|gene_control=="1",na.rm=T)
case<-gene_case$sum_variants
control<-gene_control$sum_variants
ttest<-t.test(case,control)
sink(paste0(name,"_BC_ttest.log"))
print(ttest)
sink()
}
```

T-test per variant:
```

```

