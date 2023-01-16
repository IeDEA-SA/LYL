# Cox Model 3, all sources (HOS,OPD,MED): combining men and women, adjusting for psychiatric comorbidities
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

# which disorders to analyse per interval
which_mhd_noph <- c("first_any_admission_organic","first_any_admission_substance_use_disorder","first_any_admission_psychotic","first_any_admission_mood_disorder",
                    "first_any_admission_anxiety","first_any_admission_asc_developmental_disorders","first_any_admission_asc_personality_disorder")

# exposures to control for
which_mhd <- c("first_any_admission_organic","first_any_admission_substance_use_disorder","first_any_admission_psychotic","first_any_admission_mood_disorder",
              "first_any_admission_anxiety","first_any_admission_asc_developmental_disorders","first_any_admission_asc_behavioural_disorders",
              "first_any_admission_asc_physical_factors","first_any_admission_asc_personality_disorder")


load(file=paste(filepath_read,"/data_surv.RData",sep=''))

period_cutoffs <- c(1,2)      # cutoffs for splitting of follow-up

print("Model:")
cformula <- as.formula(paste0("Surv(start,stop,event) ~ strata(sex) + strata(age_enrol_cat) + strata(calyear_cat) + ",paste0(paste0("exp_",which_mhd),collapse=" + ")))
print(cformula)

DT <- copy(data_surv)

# splitting over all exposures in model
for(v in which_mhd)
{
  DT[,diag_time:=as.numeric(get(v)-enrol_d)]
  DT <- timeSplit_DT(X=DT,vars_date="diag_time",vars_tu=paste0("exp_",v),event="event",start_date="start",stop_date="stop",id_var="patient",print_out=FALSE)
  DT[,diag_time:=NULL]
}

# splitting follow-up into intervals
DT_split <- data.table(survSplit(Surv(start,stop,event)~.,cut=365.25*period_cutoffs,episode="fup_period",data=DT))

df_aHRs <- NULL

for(i in 1:(length(period_cutoffs)+1))
{
  vect_aHRs <- NULL
  print("********")
  print(i)
  creg_fup <- coxph(cformula,data=DT_split[fup_period==i])
  print(cox.zph(creg_fup))
  
  Z <- summary(creg_fup)
  X <- Z$conf.int
  
  for(v in which_mhd_noph)
  {
    y <- X[rownames(X)==paste0("exp_",v)]
    vect_aHRs <- c(vect_aHRs,paste0(format(round(y[1],2),nsmall=2)," [",
                                    format(round(y[3],2),nsmall=2),"; ",
                                    format(round(y[4],2),nsmall=2),"]"))
    rm(y)
  }
  df_aHRs <- cbind(df_aHRs,vect_aHRs)
  rm(vect_aHRs)
}

df_aHRs <- data.frame(df_aHRs)
df_aHRs <- cbind(which_mhd_noph,df_aHRs)
colnames(df_aHRs) <- c("exposure",0,period_cutoffs)


write_xlsx(df_aHRs,path=file.path(filepath_tables,"aHRs_any_admission_table_splitfup_adjusted.xlsx"))

toc()