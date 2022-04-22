#Load phenotype data
setwd("/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/R_analysis/")
load("/proj/sens2017538/nobackup/UKBB_41143_Data/All.Phenos.20210212.RData")

#Creating variables of interest:
id <- phenos$f.eid
sex <- phenos$f.31.0.0
age <- phenos$f.21003.0.0 #age at recruitment
yob <- phenos$f.34.0.0   #year of birth

#Get BC/OC case status and age of diagnosis:

#Load R scripts
source("/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/R_analysis/GetDiseasesUKB.R")
source("/proj/sens2017538/nobackup/Exjobb/Josefin/Annotation/R_analysis/GetDiseasesUKB_ages.R")

# Define the codes for the disease you are interested in, here BC: 
ICD10 <- c("C50") #C56 for OC
ICD9 <- c("174") #183 for OC
SelfRep <- c(1002) #1039 for OC

#Run function from loaded scripts (OC instead of BC for ovarian cancer analysis)
BC <- GetDiseasesUKB(phenos, ICD10, ICD9, SelfRep)
BC_age <- GetDiseasesUKB_ages(temp, ICD10, ICD9, SelfRep, yob, BC)

BC_table <- data.frame(id, sex, age, yob, BC, BC_age)

#Genotypes for each individual at each variant position:
BC_variants <- read.table("BC_variants.raw", header=T, stringsAsFactors=F, sep="")
BC_full <- merge(BC_table, BC_variants, by.x="id", by.y="IID")
write.table(BC_full, "BC_variants_status.RData", col.names=T, row.names=F, sep="\t")
