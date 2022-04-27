library(bigsnpr)
library(bigreadr)

#Make IID column into rownames so only the models will be looped over
library(tidyverse)
setwd("/castor/project/proj_nobackup/Exjobb/Josefin/PRS/LDpred2/breast_cancer") #/ovarian_cancer for OC
load("PRS_score_table.RData")
PRS_score <- output %>% remove_rownames %>% column_to_rownames(var="IID")

#Load validation cohort and get case/control status variable (y)
load("/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/validation_cohort.RData")
fam <- read.table("ukb41143.BC.fam", header=F) #ukb41143.OC.fam for OC
y <- fam[,6]
#Change case/control status from plink format to R format
y[y=="1"]<-0
y[y=="2"]<-1
y[y=="-9"]<-NA

#Load covariates (BMI, year of birth, age at recruitment)
cov<- read.table("BC_covariates.txt", header=T, sep=" ",stringsAsFactors=F) #OC_covariates.txt for OC

#OC instead of BC for ovarian cancer analysis
cov$BC[cov$BC=="1"]<-0
cov$BC[cov$BC=="2"]<-1
cov$BC[cov$BC=="-9"]<-NA
BC<-cov[c("BC")]
yob<-cov[c("yob")]
age<-cov[c("age")]
BMI<-cov[c("BMI")]

#See how params variable is loaded in beta_grid.R

params$score <- apply(PRS_score[ind.val, ], 2, function(x) {
  if (all(is.na(x))) return(NA)
  summary(glm(y[ind.val] ~ x, family = "binomial"))$coef["x", 3]
})

#Test again adding covariates to the logisitic regression model.
params$score_cov <- apply(PRS_score[ind.val, ], 2, function(x) {
  if (all(is.na(x))) return(NA)
  summary(glm(BC[ind.val,] ~ x+yob[ind.val,]+age[ind.val,]+BMI[ind.val,], family = "binomial"))$coef["x", 3]
})
