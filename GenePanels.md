# Investigating gene panels for breast and ovarian cancer
To investigate potential disease causing variants in genes associated with breast and ovarian cancer, LoF variants were identified in the UKB WES data. Variants within 32 genes for BC (genepanel_BC.txt) and 25 genes for OC (genepanel_OC.txt) were annotated. The number of LoF variants were then compared between individuals with BC/OC and controls to to check the efficacy of the current genetical testing approaches. 

For annotation the Ensembl Variant Effect Predictor version 99 (VEP/99) was used:

McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GRS, Thormann A, et al. The Ensembl Variant Effect Predictor. Genome Biol. 2016 Dec;17(1):122.


# Annotation of WES data
The plink files containing WES data for 200K individuals from UKB were converted to VCF files to be used as input for VEP annotation. To save space in the project directory genotype information from the WES data was filtered out in the same step. The script plink2vcf.sh was used for this. The variants were then annotated using VEP/99 (see VEP.sh).

The VEP output files were filtered to only include variants found in genes that have known association with breast and ovarian cancer (filterVEP.sh). The filtered files were then used as input for filtering the WES plink files from UKB to only contain genotypes for these variants (filterPLINK.sh).

# Extracting phenotype information from UKB
To determine the case/control status for statistical testing, disease status was extracted from phenotype files from UKB. Script to extract disease information was provided by Åsa Johansson. 


To get genotype information for 

# Statistical testing
To explore if there are significant differences between cases and controls ran two-sample t-test using the t.test function in R/4.1.1. This was done for all variants in the gene panel, variants per gene from the gene panel and each variants separately. 



```

```


