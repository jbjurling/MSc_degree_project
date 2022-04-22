GetDiseasesUKB_ages <- function(temp, I10, I9, SR, yob, BC){

ind1 <- which(BC == 1)  ## Only loop over individuals with a diagnose

### Self reported

##Data-Field 20007 Description:    Interpolated Age of participant when cancer first diagnosed
indC <- which(strtrim(names(temp), 7) %in% c("f.20001"))  ### diagnose filed in UKB
indA <- which(strtrim(names(temp), 7) %in% c("f.20007"))  ### Age field UKB
names <- SR
BC_age_SR <- replicate(NROW(BC), NA)
for(i in ind1){   ## Loop over particpants with BC
    ind <- which(temp[i, indC] %in% names)
    if(NROW(ind) > 0){
    BC_age_SR[i] <- min(temp[i, indA[ind]], na.rm=T)
}}
BC_age_SR[which(!(BC_age_SR > 0 & BC_age_SR < 100))] <- NA

print(paste("There was",  NROW(which(!is.na(BC_age_SR))),  "ages identified in the SR"))

#Cancer registry:
#Data-Field 40008, Description:    Age at cancer diagnosis, matches 40006
indC10 <- which(strtrim(names(temp), 7) %in% c("f.40006"))  ### ICD10 Data fields i UK biobank.
indC9 <- which(strtrim(names(temp), 7) %in% c("f.40013"))  ### ICD9 Data fields i UK biobank
indD <- which(strtrim(names(temp), 7) %in% c("f.40005"))  ### Date (ICD9/ICD10) fields i UK biobank. Anv<C3><A4>n
indA <- which(strtrim(names(temp), 7) %in% c("f.40008"))  ### Age (ICD9/ICD10) fields i UK biobank.
names <- c(I9, I10);nc <- nchar(names);nc <- nc[which(!duplicated(nc))]  ## The ICD-9 and 10 names and the number of characters in the names
BC_age_CR <- replicate(NROW(temp), NA)
for(i in ind1){   ## Loop over particpants with registry data
    for(j in nc){    ## There are different diagnoses with different number of letters 
    ind <- which(strtrim(as.matrix(temp[i,indC9]), j) %in% names)     ## Which ICD9 match the diganose
   ind <- c(ind, which(strtrim(as.matrix(temp[i, indC10]), j) %in% names)) ## Which ICD10 match the diganose
   if(NROW(ind) > 0 ){BC_age_CR[i] <- min(temp[i, indA[ind]], na.rm=T)}   ## Choose the lowest age matching the diagnose
}}


print(paste("There was",  NROW(which(!is.na(BC_age_CR))),  "ages identified in the Cancer registry"))

## Cause and age at death  -  Same as cancer but with a different variable for age
indC10 <- which(strtrim(names(temp), 7) %in% c("f.40001", "f.40002" ) )
indA <- which(strtrim(names(temp), 7) %in% c("f.40007") )  ### age of death
BC_age_DR <- replicate(NROW(temp), NA)
for(i in ind1){   ## Loop over particpants with registrty data
      for(j in nc){
            ind <- which(strtrim(as.matrix(temp[i, indC10]), j) %in% names)
            if(NROW(ind) > 0 ){ 
        #            print(c(i,NROW(ind), temp[i, indA]) )
                     BC_age_DR[i] <- min(temp[i, indA], na.rm=T)
            }
}}

print(paste("There was",  NROW(which(!is.na(BC_age_DR))),  "ages identified in the Death registry"))

### 41202/41204/dates of first diagnose (main and secondary are merged)
indC10 <- which(strtrim(names(temp), 7) %in% c("f.41270"))  ### ICD10 Data fields i UK biobank.
indD10 <- which(strtrim(names(temp), 7) %in% c("f.41280"))  ### Date first diagnose
indC9 <- which(strtrim(names(temp), 7) %in% c("f.41271"))  ## ICD9 diagnose
indD9 <- which(strtrim(names(temp), 7) %in% c("f.41281"))   ## Date first dianose

BC_age_IH <- replicate(NROW(temp), NA)
for(i in ind1){   ## Loop over particpants with registrty data
    for(k in nc){
        ind <- indD9[which(strtrim(as.matrix(temp[i, indC9]), k) %in% names)]
        ind <- c(ind, indD10[which(strtrim(as.matrix(temp[i, indC10]), k) %in% names)])
    #   year_temp <- {}
        date_temp <- {}
        if(NROW(ind) > 0 ){
                for(j in ind){
                              date_temp <- c(date_temp, (temp[i, j]-as.Date(paste(yob[i], "-01-01", sep=""))) / 365.25) ### This gives us the year of diagnose
                 }
        
        date_temp <- min(date_temp, na.rm=T)
        BC_age_IH[i] <- date_temp
    }
   
}}
BC_age_IH[which(BC_age_IH > 100)] <- NA

print(paste("There was",  NROW(which(!is.na(BC_age_IH))),  "ages identified in the in-hospital registry"))

BC_age <- BC_age_IH
#print(NROW( which(BC_age_DR<BC_age)))
BC_age[which(BC_age_DR<BC_age | is.na(BC_age))] <- BC_age_DR[which(BC_age_DR < BC_age | is.na(BC_age))]
#print(NROW( which(BC_age_CR<BC_age)))
BC_age[which(BC_age_CR < BC_age | is.na(BC_age))] <- BC_age_CR[which(BC_age_CR < BC_age | is.na(BC_age))]
#print(NROW( which(BC_age_SR<BC_age)))
BC_age[which(BC_age_SR < BC_age | is.na(BC_age))] <- BC_age_SR[which(BC_age_SR < BC_age | is.na(BC_age))]

return(BC_age)

}
