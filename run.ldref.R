library(bigsnpr)
library(bigreadr)
setwd("/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/breast_cancer") #remove /breast_cancer for OC analysis

## Information for the variants provided in the LD reference
map_ldref <- readRDS("LD_ref/map.rds") #add breast_cancer/ for OC analysis

## Breast cancer summary statistics
sumstats <- bigreadr::fread2("/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/gwas_sum/breast_cancer/BC_gwas_summary_filtered_UKB.txt") #ovarian_cancer/OC_gwas_summary_filtered_UKB.txt for OC
sumstats$n_eff <- 4 / (1 / 122977 + 1 / 105974) #sumstats$n_eff <- 4 / (1 / 22406 + 1 / 40941) for OC

#Matching variants between GWAS summary statistics and genotype data
info_snp <- snp_match(sumstats, map_ldref)
(info_snp <- tidyr::drop_na(tibble::as_tibble(info_snp)))

sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))
sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2))

is_bad <-
  sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05

df_beta <- info_snp[!is_bad, ]

#Create SNP correlation matrix
tmp <- tempfile(tmpdir = "tmp-data")

for (chr in 1:22) {
  
  cat(chr, ".. ", sep = "")
  
  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'map_ldref'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  ## indices in 'corr_chr'
  ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))
  
  corr_chr <- readRDS(paste0("LD_ref/LD_with_blocks_chr", chr, ".rds"))[ind.chr3, ind.chr3] #breast_cancer/LD_ref/LD_with_blocks_chr for OC
  
  if (chr == 1) {
    corr <- as_SFBM(corr_chr, tmp)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}

#Create LD matrix
for (chr in 1:22) {

  if (chr == 1) {
    ld <- Matrix::colSums(corr_chr^2)
    } else {
      ld <- c(ld, Matrix::colSums(corr_chr^2))
    }
}

#setwd("/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/ovarian_cancer/") for OC

save.image("working.environment.ldref.RData")
