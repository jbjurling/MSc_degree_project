GetDiseasesUKB <- function(temp, I10, I9, SR){

### ICD10 OBS!! might need to add variables before running script.
ind<-which(strtrim(names(temp),7) %in% c("f.40006","f.41202","f.41204","f.40001", "f.40002","f.41270"))   ### ICD10 Data fields i UK biobank.
names<-I10;nc<-nchar(names);nc<-nc[which(!duplicated(nc))]
ind1<-{}
for(i in ind){
     for(j in nc){
         ind1<-c(ind1,which(strtrim((temp[,i]),j) %in% names))
 }
}
ind1<-ind1[which(!duplicated(ind1))]
print(paste("There was",  NROW(ind1), "cases identified in the ICD10-coding"))
BC<-replicate(NROW(temp),1)## Controls.
BC[ind1]<-2;

##ICD9

## ICD9
ind<-which(strtrim(names(temp),7) %in% c("f.41203", "f.41205","f.40013","f.41271"))   ### ICD9 Data fields i UK biobank.
names<-I9;nc<-nchar(names);nc<-nc[which(!duplicated(nc))]
ind1<-{}
for(i in ind){
 for(j in nc){
         ind1<-c(ind1,which(strtrim((temp[,i]),j) %in% names))
 }
}


print(paste("There was",  NROW(ind1), "cases identified in the ICD9-coding"))
BC[ind1]<-2;

### Selfreported
ind<-which(strtrim(names(temp),7)=="f.20001")
names<-SR
ind1<-{}
for(i in ind){
    ind1<-c(ind1,which(temp[,i] %in% names))
}

print(paste("There was",  NROW(ind1), "cases identified in the Self-Reported"))
BC[ind1]<-2;


return(BC)

}
