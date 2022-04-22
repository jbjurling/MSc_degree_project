setwd("/proj/sens2017538/nobackup/UKBB_IMP_DOSAGE_V3_bim_bed")
library(bigsnpr)
library(bigassertr)

##snp_readBed to read plink bed file. Created a function to read plink bed files based on the snp_readBed function. 
  #The corresponding ".bim" and ".fam" needs to be in the same directory.
 
#Create variable for individual IDs to include (only females)
famfile <- read.table("ukb41143.fam", header=F)
colnames(famfile) <- c("FID", "IID", "p", "m", "sex", "pheno")
ID <- famfile$FID[famfile$sex == 2] 

#Create variable for SNPs to include (used GWAS summary files only filtered for maf>0.01)
snps <- read.table("/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/gwas_sum/breast_cancer/BC_gwas_summary.txt",header=T,sep=",") #../ovarian_cancer/OC_gwas_summary_filtered_maf001.txt for OC
SNPnames <- snps$var_name

##Function to read plink files to a bigSNP array
JB_snp_readBed <- function(bedfile) {
  
  #Load bed file
  obj.bed <- bed(bedfile)
  
  # Check if backingfile already exists 
  backingfile <- sub_bed(bedfile) #add ".OC" for OC
  assert_noexist(paste0(backingfile, ".bk"))
  
  # Get other files
  fam <- obj.bed$fam; rownames(fam) <- rows_along(fam)
  bim <- obj.bed$map; rownames(bim) <- rows_along(bim)
  
  # Create variable in bim to match against SNPnames
  bim$var_name <- paste(bim$chromosome, bim$physical.pos, sep="_")
  
  # Create index for individuals and SNPs to load (SNPs for all chr in same variable)
  ind.row <- which(fam$sample.ID %in% ID)
  ind.col <- which(bim$var_name %in% SNPnames)
  
  #Include only individuals and SNPs in fam/bim files
  fam <- fam[ind.row, ]
  bim <- bim[ind.col, ]
  
  # Prepare Filebacked Big Matrix
  bigGeno <- FBM.code256(
    nrow = nrow(fam),
    ncol = nrow(bim),
    code = CODE_012,
    backingfile = backingfile,
    init = NULL,
    create_bk = TRUE
  )
  
  # Fill the FBM from bedfile
  reach.eof <- readBin(bedfile, what = 1L) ##changed from readbina according to solution from privefl (https://github.com/privefl/bigsnpr/issues/101)
  # what=1L - return the value as a long (integer type), not an int,
  if (!reach.eof) warning("EOF of bedfile has not been reached.")
  
  # Create the bigSNP object
  snp.list <- structure(list(genotypes = bigGeno, fam = fam, map = bim),
                        class = "bigSNP")
  
  # save it and return the path of the saved object
  rds <- sub_bk(bigGeno$backingfile, ".rds")
  saveRDS(snp.list, rds)
  rds
}

##Command to run function for all chromosomes
list <- read.table("/proj/sens2017538/nobackup/Exjobb/Josefin/Scripts/chr_v2.txt")

for (chr in list[,1]){
  JB_snp_readBed(paste0(chr, ".bed"))
}
