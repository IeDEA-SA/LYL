# Cox Model 3, medication (MED_ATC_N): adjusting for sex, including psychiatric comorbidity
# Separate analyses for all-cause mortality, natural death, and unnatural death 
# Regressors included apart from the exposures:                                      
#  sex, age group, calendar year (time-updated)
# Run time: 2 minutes

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
which_mhd<- which_mhd<-  c("first_substance_use_med","first_antipsychotic","first_antidepressant","first_anxiolytic")

load(file=paste(filepath_read,"/data_surv.RData",sep=''))

print("Model:")
cformula <- as.formula(paste0("Surv(start,stop,event) ~ strata(sex) + strata(age_enrol_cat) + strata(calyear_cat) + ",paste0(paste0("exp_",which_mhd),collapse=" + ")))
print(cformula)

df_aHRs <- data.frame(NULL)
creg_list_ac <- NULL
creg_list_nd <- NULL
creg_list_ud <- NULL

DT_ac <- copy(data_surv)

DT_nd <- copy(DT_ac)
DT_nd[cod2%in%c("Unnatural","Unknown"),event:=0]
DT_ud <- copy(DT_ac)
DT_ud[cod2%in%c("Natural","Unknown"),event:=0]

# splitting over all exposures in model
for(v in which_mhd)
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

# Cox regression
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

creg_list_ac <- Z_ac
creg_list_nd <- Z_nd
creg_list_ud <- Z_ud

X_ac <- Z_ac$conf.int
X_ac <- X_ac[grepl("exp",row.names(X_ac)),]
X_nd <- Z_nd$conf.int
X_nd <- X_nd[grepl("exp",row.names(X_nd)),]
X_ud <- Z_ud$conf.int
X_ud <- X_ud[grepl("exp",row.names(X_ud)),]

zph_ac <- zph_ac[grepl("exp",row.names(zph_ac)),]
zph_nd <- zph_nd[grepl("exp",row.names(zph_nd)),]
zph_ud <- zph_ud[grepl("exp",row.names(zph_ud)),]

df_aHRs <- rbind(df_aHRs,
                 data.frame(disorder=which_mhd,sex="Both",cause="all",estimate=X_ac[,1],
                            lcl=X_ac[,3],ucl=X_ac[,4],
                            p_ph_exp=zph_ac[,3]),
                 data.frame(disorder=which_mhd,sex="Both",cause="Natural",estimate=X_nd[,1],
                            lcl=X_nd[,3],ucl=X_nd[,4],
                            p_ph_exp=zph_nd[,3]),
                 data.frame(disorder=which_mhd,sex="Both",cause="Unnatural",estimate=X_ud[,1],
                            lcl=X_ud[,3],ucl=X_ud[,4],
                            p_ph_exp=zph_ud[,3]))

df_aHRs <- data.table(df_aHRs)
df_aHRs[,disorder:=factor(disorder,levels=which_mhd)]
df_aHRs[,cause:=factor(cause,levels=c("all","Natural","Unnatural"))]
setorder(df_aHRs,"disorder","cause","sex")

write_xlsx(df_aHRs,path=file.path(filepath_tables,paste0("aHRs_model3_med.xlsx")))
save(creg_list_ac,file=file.path(filepath_write,"obj_model3_ac_med.RData"))
save(creg_list_nd,file=file.path(filepath_write,"obj_model3_nd_med.RData"))
save(creg_list_ud,file=file.path(filepath_write,"obj_model3_ud_med.RData"))

toc()
