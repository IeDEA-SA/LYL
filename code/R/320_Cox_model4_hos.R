# Cox Model 4, hospital admissions (HOS): for men and women and separately, including psychiatric comorbidity
# pschiatric morbidity = comorbidity from ANY source (OPD, HOS, MED), not just hospital admissions
# Separate analyses for all-cause mortality, natural death, and unnatural death 
# Regressors included apart from the exposures:                                      
#  age group, calendar year (time-updated)
# Run time: 20 minutes

library(data.table)
library(tictoc)
library(survival)
library(writexl)

tic()

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfAc_excess_mortality/processed"
filepath_tables <- "C:/ISPM/HomeDir/HIV-mental disorders/AfAc_excess_mortality/output/tables/model output"
filepath_write <- "C:/ISPM/Data/HIV-mental disorders/AfAc_excess_mortality/regression objects" 
source(file="C:/ISPM/HomeDir/HIV-mental disorders/AfAc_excess_mortality/LYL/code/R/utils/timeSplit_DT.R") # found in Utils folder of github repository

# exposures to analyse
which_mhd <- c("first_hos_organic","first_hos_substance_use_disorder","first_hos_psychotic","first_hos_mood_disorder",
              "first_hos_anxiety","first_hos_asc_developmental_disorders","first_hos_asc_behavioural_disorders",
              "first_hos_asc_physical_factors","first_hos_asc_personality_disorder")
which_comorb <- gsub("first_hos_","first_any_admission_",which_mhd)

load(file=paste(filepath_read,"/data_surv.RData",sep=''))

df_aHRs <- data.frame(NULL)
creg_list_ac <- NULL
creg_list_nd <- NULL
creg_list_ud <- NULL

DT_ac <- copy(data_surv)
DT_nd <- copy(DT_ac)
DT_nd[cod2%in%c("Unnatural","Unknown"),event:=0]
DT_ud <- copy(DT_ac)
DT_ud[cod2%in%c("Natural","Unknown"),event:=0]

# splitting over all exposures + comorbidities in one go (~2 minutes)
for(v in c(which_mhd,which_comorb))
{
  DT_ac[,diag_time:=as.numeric(get(v)-enrol_d)]
  DT_ac <- timeSplit_DT(X=DT_ac,vars_date="diag_time",vars_tu=paste0("exp_",v),event="event",start_date="start",stop_date="stop",id_var="patient",print_out=FALSE)
  DT_ac[,diag_time:=NULL]
  
  DT_nd[,diag_time:=as.numeric(get(v)-enrol_d)]
  DT_nd <- timeSplit_DT(X=DT_nd,vars_date="diag_time",vars_tu=paste0("exp_",v),event="event",start_date="start",stop_date="stop",id_var="patient",print_out=FALSE)
  DT_nd[,diag_time:=NULL]
  
  DT_ud[,diag_time:=as.numeric(get(v)-enrol_d)]
  DT_ud <- timeSplit_DT(X=DT_ud,vars_date="diag_time",vars_tu=paste0("exp_",v),event="event",start_date="start",stop_date="stop",id_var="patient",print_out=FALSE)
  DT_ud[,diag_time:=NULL]
}

for(ggg in c("Male","Female"))
{
  for(v in which_mhd)
  {
    # adjusting for mental illnesses of different types than the one being analyzed, from any setting (not just hospital)
    w <- gsub("first_hos_","",v)
    comorb_vect <- paste0("exp_",which_comorb[-grep(w,which_comorb)])
    rm(w)
    
    print("Model:")
    cformula <- as.formula(paste0("Surv(start,stop,event) ~ strata(age_enrol_cat) + strata(calyear_cat) + ",paste0("exp_",v," + "),paste0(comorb_vect,collapse=" + ")))
    print(cformula)
    print("*****")
    
    # Cox regression
    creg_ac <- coxph(cformula,data=DT_ac[sex==ggg])
    creg_nd <- coxph(cformula,data=DT_nd[sex==ggg])
    creg_ud <- coxph(cformula,data=DT_ud[sex==ggg])

    # testing proportional hazards
    zph_ac <- cox.zph(creg_ac)$table
    zph_nd <- cox.zph(creg_nd)$table
    zph_ud <- cox.zph(creg_ud)$table
    
    # storing the regression output
    Z_ac <- summary(creg_ac)
    Z_nd <- summary(creg_nd)
    Z_ud <- summary(creg_ud)
  
    creg_list_ac[[paste(ggg,v,sep="_")]] <- Z_ac
    creg_list_nd[[paste(ggg,v,sep="_")]] <- Z_nd
    creg_list_ud[[paste(ggg,v,sep="_")]] <- Z_ud
  
    X_ac <- Z_ac$conf.int
    X_ac <- X_ac[grepl("exp_first_hos",row.names(X_ac)),]
    X_nd <- Z_nd$conf.int
    X_nd <- X_nd[grepl("exp_first_hos",row.names(X_nd)),]
    X_ud <- Z_ud$conf.int
    X_ud <- X_ud[grepl("exp_first_hos",row.names(X_ud)),]
    
    zph_ac <- zph_ac[grepl("exp_first_hos",row.names(zph_ac)),]
    zph_nd <- zph_nd[grepl("exp_first_hos",row.names(zph_nd)),]
    zph_ud <- zph_ud[grepl("exp_first_hos",row.names(zph_ud)),]
  
    df_aHRs <- rbind(df_aHRs,
                     data.frame(disorder=v,sex=ggg,cause="all",estimate=X_ac[1],
                                lcl=X_ac[3],ucl=X_ac[4],
                                p_ph_exp=zph_ac[3]),
                     data.frame(disorder=v,sex=ggg,cause="Natural",estimate=X_nd[1],
                                lcl=X_nd[3],ucl=X_nd[4],
                                p_ph_exp=zph_nd[3]),
                     data.frame(disorder=v,sex=ggg,cause="Unnatural",estimate=X_ud[1],
                                lcl=X_ud[3],ucl=X_ud[4],
                                p_ph_exp=zph_ud[3]))
  }
}
df_aHRs <- data.table(df_aHRs)
df_aHRs[,disorder:=factor(disorder,levels=which_mhd)]
df_aHRs[,cause:=factor(cause,levels=c("all","Natural","Unnatural"))]
setorder(df_aHRs,"disorder","cause","sex")

write_xlsx(df_aHRs,path=file.path(filepath_tables,paste0("aHRs_model4_hos.xlsx")))
save(creg_list_ac,file=file.path(filepath_write,"obj_model4_ac_hos.RData"))
save(creg_list_nd,file=file.path(filepath_write,"obj_model4_nd_hos.RData"))
save(creg_list_ud,file=file.path(filepath_write,"obj_model4_ud_hos.RData"))

toc()
