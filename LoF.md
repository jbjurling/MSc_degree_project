# Investigating gene panels for breast and ovarian cancer


For annotation the Ensembl Variant Effect Predictor version 99 (VEP/99) was used:

McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GRS, Thormann A, et al. The Ensembl Variant Effect Predictor. Genome Biol. 2016 Dec;17(1):122.


# Annotation of WES data
The plink files containing WES data for 200K individuals from UKB were converted to VCF files to be used as input for VEP annotation. To save space in the project directory genotype information from the WES data was filtered out in the same step. The script plink2vcf.sh was used for this. The variants were then annotated using VEP/99 (see VEP.sh).

The vep output files were filtered to only include variants found in genes that have known association with breast and ovarian cancer (filterVEP.sh). The filtered files were then used as input for filtering the WES plink files from UKB to only contain genotypes for these variants (filterPLINK.sh).


