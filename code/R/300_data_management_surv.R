# Preparing data for surival analysis

library(data.table)
library(tictoc)
library(survival)

tic()

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfAc_excess_mortality/processed"
filepath_write <- "C:/ISPM/Data/HIV-mental disorders/AfAc_excess_mortality/processed" 

min_age <- 15
max_age <- 85

load(file=paste(filepath_read,"/long_table.RData",sep=''))

data_surv <- copy(BAS_MHD_new)

print(paste0("Starting number of patients: ",uniqueN(data_surv,"patient")))

# switching origin of the regressions to start of time at risk (i.e. maximum between enrolment date and 2011-01-01)
#data_surv[,enrol_d:=start]

# counting process format, one line per patient, on 'age' scale, in days
data_surv[,`:=`(start=as.numeric(start-birth_d),
                stop=as.numeric(end-birth_d))]

# left-truncataion/right-censoring (age)
data_surv <- data.table(survSplit(Surv(start,stop,event)~.,data=data_surv,cut=c(min_age*365.25,max_age*365.25),episode="i"))
data_surv <- data_surv[i==2]
data_surv[,i:=NULL]

print(paste0("*After excluding patients with no follow-up between ages ",min_age," and ", max_age,": ",uniqueN(data_surv,"patient")))

# switching to 'calendar year' scale = the number of days since 1970-01-01
data_surv[,`:=`(start=as.numeric(birth_d+start),
                stop=as.numeric(birth_d+stop))]
calyear_cutoffs <- sort(as.numeric(as.Date(c(paste0(c(2011,2014,2017,2021),"-01-01"),"2020-03-15"))))   # cutoff dates (note: there is no follow-up past June 2020)
data_surv <- data.table(survSplit(Surv(start,stop,event)~.,data=data_surv,cut=calyear_cutoffs,episode="calyear_cat"))

# switching to 'days since enrolment' scale for the analysis
data_surv[,`:=`(start=(start-as.numeric(enrol_d)),
                stop=(stop-as.numeric(enrol_d)))]

# cleaning up cause of death AFTER censoring and splitting
data_surv[event==0,cod2:="Alive"]
data_surv[,cod2:=factor(cod2,levels=c("Alive","Natural","Unnatural","Unknown"))]
data_surv[,table(cod2,useNA="a")]

# rescaling so everyone starts at time 0 (no left-truncation)
data_surv[,start_min:=min(start),by="patient"]
data_surv[,enrol_d:= enrol_d + start_min]
data_surv[,`:=`(start=start-start_min,stop=stop-start_min)]
data_surv[,start_min:=NULL]

# formatting risk factors
data_surv[,`:=`(sex=factor(sex,levels=c("Male","Female")),
                calyear_cat=factor(calyear_cat),
                age_enrol_cat=cut(as.numeric(enrol_d-birth_d)/365.25,breaks=c(0,25,40,55,75,Inf),right=FALSE))]

save(data_surv,file=file.path(filepath_write,"data_surv.RData"))

toc()
