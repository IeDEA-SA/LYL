# point estimates for excess life-years lost, medication
# output: R table "elyl_any_med.RData"
# ~1 hour runtime

library(data.table)
library(survival)
library(ggplot2)
library(lillies)
library(tictoc)

tic()

# filepath_read <- "~/MHD/R/input/AfAc"
# filepath_write <- "~/MHD/R/output/lyl"
# source("~/MHD/R/Code/Utils/timeSplit_DT.R")

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfAc_excess_mortality/processed"
filepath_write <- "C:/ISPM/Data/HIV-mental disorders/AfAc_excess_mortality/lillies/point_estimates" 
source(file="C:/ISPM/HomeDir/HIV-mental disorders/R/Code/utils/timeSplit_DT.R")   # found in Utils folder of github repository

which_mhd <- c("first_any_med","first_substance_use_med","first_antidepressant","first_anxiolytic","first_antipsychotic")
which_sex <- c("Men","Women")

min_age <- 15
max_age <- 85

save_filename <- "elyl_med.RData"

load(file=paste(filepath_read,"/long_table.RData",sep=''))
data_surv <- copy(BAS_MHD_new)

# removing patients with unknown sex
data_surv <- data_surv[data_surv$sex=="Male"|data_surv$sex=="Female"]

#removing patients with NA in birth_d
data_surv<-data_surv[!is.na(birth_d)]

# counting process format, one line per patient, on 'age' scale
data_surv[,`:=`(start=as.numeric(start-birth_d)/365.25,
                stop=as.numeric(end-birth_d)/365.25)]

# left-truncation/right-censoring (age)
data_surv <- data.table(survSplit(Surv(start,stop,event)~.,data=data_surv,cut=c(min_age,max_age),episode="i"))
data_surv <- data_surv[i==2]
data_surv[,i:=NULL]
print(paste0("After excluding patients with no follow-up between ages 15 and 85: ",uniqueN(data_surv,"patient")))

# cleaning up AFTER censoring
data_surv[event==0,cod2:="Alive"]

data_surv[,cod2:=factor(cod2,levels=c("Alive","Natural","Unnatural","Unknown"))]

data_surv[,table(cod2,useNA="a")]

# looping over gender & types of mental disorders
LYL_diff_df <- data.frame(NULL)
for(ggg in which_sex)
{
  for(v in which_mhd)
  {
    tic(paste0(ggg,", ",v))
    DT <- copy(data_surv)
    if(ggg=="Men")
      DT <- DT[sex=="Male"]
    if(ggg=="Women")
      DT <- DT[sex=="Female"]
    
    # age at first psychiatric medication of any kind
    DT[,age_any_treat:=as.numeric(first_any_med-birth_d)/365.25]
    
    setnames(DT,v,"treat_d")
    
    # age at mhd treatment
    DT[,age_treat:=as.numeric(treat_d-birth_d)/365.25]
    DT[age_treat>=stop,age_treat:=NA]                              # excluding disorders occuring post-censoring  
    DT[!is.na(age_treat) & age_treat<start,age_treat:=start]       # moving forward prevalent MHD treatments forward to cohort entry
    age_treat <- DT[!is.na(age_treat),age_treat]
    
    # splitting time-at-risk into exposed/unexposed
    DT_exp <- DT[!is.na(age_treat)]
    DT_unexp <- timeSplit_DT(X=DT,vars_date="age_any_treat",vars_tu="cens",event="event",
                             start_date="start",stop_date="stop",id_var="patient",print_out=FALSE)
    DT_unexp <- DT_unexp[cens==0]
    DT_unexp[,cens:=NULL]
    DT_unexp[event==0,cod2:="Alive"]
    
    # computing excess LYL, exposed vs. unexposed
    LYL_exp <- lyl_range(data=DT_exp,t0=age_treat,t=stop,status=cod2,age_begin=min_age,age_end=max_age-1,tau=max_age)
    LYL_unexp <- lyl_range(data=DT_unexp,t0=start,t=stop,status=cod2,age_begin=min_age,age_end=max_age-1,tau=max_age)
    LYL_diff <- lyl_diff(LYL_exp,LYL_unexp,weights=age_treat)
    
    #Create DF for results
    
    LYL_diff_df <- rbind(LYL_diff_df,data.frame("type"=v,"sex"=ggg,LYL_diff,"nb_patients"=DT[,.N],"nb_exposed"=length(age_treat),
                                                "nb_nat_deaths"=DT[,sum(cod2=="Natural")],"nb_unnat_deaths"=DT[,sum(cod2=="Unnatural")],
                                                "nb_uk_deaths"=DT[,sum(cod2=="Unknown")],"nb_deaths"=DT[,sum(cod2!="Alive")],"py"=DT[,sum(stop-start)]))
    rm(DT_exp,DT_unexp,age_treat)
    toc()
  }
}
toc()

save(LYL_diff_df,file=file.path(filepath_write,save_filename))


