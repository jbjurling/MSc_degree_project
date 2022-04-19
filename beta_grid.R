library(bigsnpr)
library(bigreadr)
setwd("/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/breast_cancer/") #/ovarian_cancer/ for OC

load("working.environment.ldref.RData")

NCORES<-nb_cores()
(ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
                                chi2 = (beta / beta_se)^2,
                                sample_size = n_eff,
                                ncores = NCORES)))
 
h2_est <- ldsc[["h2"]]

#LDpred2-grid 
load("/proj/sens2017538/nobackup/Exjobb/Josefin/PRS/LDpred2/validation_cohort.RData") # load ind.val variable
(h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4))
(p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2))
(params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))) #sparse=some effect size is exactly 0

beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
save(beta_grid,file="beta_grid.RData")
