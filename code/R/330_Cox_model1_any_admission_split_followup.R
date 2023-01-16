# Cox Model 1, all sources (HOS,OPD,MED): combining men and women, no adjustment for psychiatric comorbidities
# All-cause mortality
# Regressors included apart from the exposure:                                      
#  sex, age group, calendar year (time-updated)
# Follow-up divided into intervals, and the analysis is done over each interval separately

library(data.table)
library(tictoc)
library(survival)
library(writexl)

tic()

ubelix <- FALSE

if(ubelix)
{
  filepath_read <- "~/MHD/R/input/AfAc"
  filepath_tables <- "~/MHD/R/output/tables"
  filepath_write <- "~/MHD/R/output/regression objects"
  source(file="~/MHD/R/Code/Utils/timeSplit_DT.R")
} else {
  filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfAc_excess_mortality/processed"
  filepath_tables <- "C:/ISPM/HomeDir/HIV-mental disorders/AfAc_excess_mortality/output/tables/model output"
  filepath_write <- "C:/ISPM/Data/HIV-mental disorders/AfAc_excess_mortality/regression objects" 
  source(file="C:/ISPM/HomeDir/HIV-mental disorders/AfAc_excess_mortality/LYL/code/R/utils/timeSplit_DT.R") # found in Utils folder of github repository
}

load(file=paste(filepath_read,"/data_surv.RData",sep=''))

# exposures to analyse
which_mhd <- c("first_any_admission","first_any_admission_organic","first_any_admission_substance_use_disorder","first_any_admission_psychotic",
               "first_any_admission_mood_disorder","first_any_admission_anxiety","first_any_admission_asc_developmental_disorders","first_any_admission_asc_personality_disorder")
period_cutoffs <- c(1,2)                 # cutoffs for splitting of follow-up

df_aHRs <- NULL

for(v in which_mhd)
{
  DT <- copy(data_surv)
  setnames(DT,v,"diag_d")
  
  # number of days between enrolment and first exposure
  DT[,diag_time:=as.numeric(diag_d-enrol_d)]
  DT <- DT[,.(patient,sex,age_enrol_cat,calyear_cat,popgrp,start,stop,event,diag_time,cod2)]
  
  # splitting at time of first exposure
  DT_ac <- timeSplit_DT(X=DT,vars_date="diag_time",vars_tu="exp",event="event",start_date="start",stop_date="stop",id_var="patient",print_out=FALSE)
  
  print(paste0("Number exposed: ",uniqueN(DT_ac[exp==1],by="patient")))
  
  cformula <- as.formula("Surv(start,stop,event) ~ exp + strata(sex) + strata(age_enrol_cat) + strata(calyear_cat)")

  # splitting follow-up into intervals
  DT_ac_split <- data.table(survSplit(Surv(start,stop,event)~.,cut=365.25*period_cutoffs,episode="fup_period",data=DT_ac))
  vect_aHRs <- NULL
  
  for(i in 1:(length(period_cutoffs)+1))
  {
    print("********")
    print(paste(v,i))
    creg_ac_fup <- coxph(cformula,data=DT_ac_split[fup_period==i])
    print(cox.zph(creg_ac_fup))
    
    Z_ac <- summary(creg_ac_fup)
    X_ac <- Z_ac$conf.int
    Y_ac <- Z_ac$coefficients
    
    vect_aHRs <- c(vect_aHRs,paste0(format(round(X_ac[row.names(X_ac)=="exp",1],2),nsmall=2)," [",
                                    format(round(X_ac[row.names(X_ac)=="exp",3],2),nsmall=2),"; ",
                                    format(round(X_ac[row.names(X_ac)=="exp",4],2),nsmall=2),"]"))

  }
  df_aHRs <- rbind(df_aHRs,vect_aHRs)
}

df_aHRs <- cbind(which_mhd,df_aHRs)
df_aHRs <- data.frame(df_aHRs)
colnames(df_aHRs) <- c("exp",0,period_cutoffs)

write_xlsx(df_aHRs,path=file.path(filepath_tables,"aHRs_any_admission_table_splitfup.xlsx"))

toc()