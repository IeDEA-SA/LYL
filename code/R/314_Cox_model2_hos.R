# Cox Model 2, hospital admissions (HOS): separate modesl for men and women, no adjustment for psychiatric comorbidities
# Separate analyses for all-cause mortality, natural death, and unnatural death 
# Regressors included apart from the exposure:                                      
#  age group, calendar year (time-updated)
# Run time when analyzing all exposures: 10+ minutes   

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
which_mhd<- c("first_hos_any","first_hos_organic","first_hos_substance_use_disorder","first_hos_psychotic",
              "first_hos_mood_disorder","first_hos_anxiety","first_hos_asc_developmental_disorders","first_hos_asc_personality_disorder",
              "first_hos_alcohol_use_disorder","first_hos_drug_use_disorder","first_hos_bipolar","first_hos_depression",
              "first_hos_generalised_anxiety_disorder","first_hos_PTSD","first_hos_eating_disorder")

load(file=paste(filepath_read,"/data_surv.RData",sep=''))

df_aHRs <- data.frame(NULL)
creg_list_ac <- NULL
creg_list_nd <- NULL
creg_list_ud <- NULL
print("--------------------")
print("Currently analyzing:")

for(v in which_mhd)
{
  print(paste0("*",v))
  for(ggg in c("Male","Female"))
  {
    DT <- data_surv[sex==ggg]
    
    setnames(DT,v,"diag_d")
    
    # number of days between enrolment and first exposure
    DT[,diag_time:=as.numeric(diag_d-enrol_d)]
    DT <- DT[,.(patient,sex,age_enrol_cat,calyear_cat,start,stop,event,diag_time,cod2)]
    DT_nd <- copy(DT)
    DT_nd[cod2%in%c("Unnatural","Unknown"),event:=0]
    DT_ud <- copy(DT)
    DT_ud[cod2%in%c("Natural","Unknown"),event:=0]
    
    # splitting at time of first exposure
    DT_ac <- timeSplit_DT(X=DT,vars_date="diag_time",vars_tu="exp",event="event",start_date="start",stop_date="stop",id_var="patient",print_out=FALSE)
    DT_nd <- timeSplit_DT(X=DT_nd,vars_date="diag_time",vars_tu="exp",event="event",start_date="start",stop_date="stop",id_var="patient",print_out=FALSE)
    DT_ud <- timeSplit_DT(X=DT_ud,vars_date="diag_time",vars_tu="exp",event="event",start_date="start",stop_date="stop",id_var="patient",print_out=FALSE)
    
    # Cox regression
    cformula <- as.formula("Surv(start,stop,event) ~ exp + strata(age_enrol_cat) + strata(calyear_cat)")
    creg_ac <- coxph(cformula,data=DT_ac)
    creg_nd <- coxph(cformula,data=DT_nd)
    creg_ud <- coxph(cformula,data=DT_ud)
    
    # testing proportional hazards
    zph_ac <- cox.zph(creg_ac)$table
    zph_nd <- cox.zph(creg_nd)$table
    zph_ud <- cox.zph(creg_ud)$table
    
    # storing the regression output
    Z_ac <- summary(creg_ac)
    Z_nd <- summary(creg_nd)
    Z_ud <- summary(creg_ud)
    
    creg_list_ac[[paste0(v,"_",ggg)]] <- Z_ac
    creg_list_nd[[paste0(v,"_",ggg)]] <- Z_nd
    creg_list_ud[[paste0(v,"_",ggg)]] <- Z_ud
    
    X_ac <- Z_ac$conf.int
    X_nd <- Z_nd$conf.int
    X_ud <- Z_ud$conf.int
    
    df_aHRs <- rbind(df_aHRs,
                     data.frame(disorder=v,sex=ggg,cause="all",estimate=X_ac[row.names(X_ac)=="exp",1],
                                lcl=X_ac[row.names(X_ac)=="exp",3],ucl=X_ac[row.names(X_ac)=="exp",4],
                                p_ph_exp=zph_ac[row.names(zph_ac)=="exp",3]),
                     data.frame(disorder=v,sex=ggg,cause="Natural",estimate=X_nd[row.names(X_nd)=="exp",1],
                                lcl=X_nd[row.names(X_nd)=="exp",3],ucl=X_nd[row.names(X_nd)=="exp",4],
                                p_ph_exp=zph_nd[row.names(zph_nd)=="exp",3]),
                     data.frame(disorder=v,sex=ggg,cause="Unnatural",estimate=X_ud[row.names(X_ud)=="exp",1],
                                lcl=X_ud[row.names(X_ud)=="exp",3],ucl=X_ud[row.names(X_ud)=="exp",4],
                                p_ph_exp=zph_ud[row.names(zph_ud)=="exp",3]))
  }
}

df_aHRs <- data.table(df_aHRs)
df_aHRs[,disorder:=factor(disorder,levels=which_mhd)]
df_aHRs[,cause:=factor(cause,levels=c("all","Natural","Unnatural"))]
setorder(df_aHRs,"disorder","cause","sex")

write_xlsx(df_aHRs,path=file.path(filepath_tables,"aHRs_model2_hos.xlsx"))
save(creg_list_ac,file=file.path(filepath_write,"obj_model2_ac_hos.RData"))
save(creg_list_nd,file=file.path(filepath_write,"obj_model2_nd_hos.RData"))
save(creg_list_ud,file=file.path(filepath_write,"obj_model2_ud_hos.RData"))

toc()
