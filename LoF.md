# Investigating LoF variants in WES UKB data

To estimate the amount of variants in well-known cancer genes that may contribute to disease the WES dataset was annotated using the the Ensembl Variant Effect Predictor version 99 (VEP/99) 
To identify potential disease causing variants in genes associated with breast and ovarian cancer, whole exome sequencing data from UKB was annotated. Variants were defined as disease 
Annotation of potential disease causing variants were made using the whole exome sequencing data from 100K individuals in UKB.


To identify potential disease causing variants in the UKB dataset, genes associated with BC and OC was identified. The 32 genes for BC and 25 genes for OC that are part of current genetic testing panels were selected for inclusion in the study. LoF variants were defined as potentially disease causing variants.
Variants within these genes were annotated using whole exome sequencing data from 200K individuals in UK Biobank.

The frecuency of these variants in the cases and controls were determined.

For annotation the Ensembl Variant Effect Predictor version 99 (VEP/99) was used:

McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GRS, Thormann A, et al. The Ensembl Variant Effect Predictor. Genome Biol. 2016 Dec;17(1):122.


# Annotation of WES data
The plink files containing WES data for 200K individuals from UKB were converted to VCF files to be used as input for VEP annotation. To save space in the project directory genotype information from the WES data was filtered out in the same step. The script plink2vcf.sh was used for this. The variants were then annotated using VEP/99 (see VEP.sh).

The vep output files were filtered to only include variants found in genes that have known association with breast and ovarian cancer (filterVEP.sh). The filtered files were then used as input for filtering the WES plink files from UKB to only contain genotypes for these variants (filterPLINK.sh).


