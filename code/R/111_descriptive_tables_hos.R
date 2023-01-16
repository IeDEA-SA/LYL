# descriptive tables, hospitalizations
# one table with overall characterstics, stratified by MHD status at end of follow-up (any/none)
# second table with breakdown of different MHDs per sex

library(data.table)
library(survival)
library(tictoc)
library(tableone)
library(writexl)

tic()

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfAc_excess_mortality/processed"
filepath_write <- "C:/ISPM/HomeDir/HIV-mental disorders/AfAc_excess_mortality/output/tables/descriptive"
source(file="C:/ISPM/HomeDir/HIV-mental disorders/AfAc_excess_mortality/LYL/code/R/utils/timeSplit_DT.R") # found in Utils folder of github repository

which_mhd <- c("first_hos_any","first_hos_organic","first_hos_substance_use_disorder","first_hos_alcohol_use_disorder","first_hos_drug_use_disorder",
               "first_hos_psychotic","first_hos_mood_disorder","first_hos_bipolar","first_hos_depression","first_hos_anxiety","first_hos_generalised_anxiety_disorder",
               "first_hos_PTSD","first_hos_asc_developmental_disorders","first_hos_eating_disorder","first_hos_asc_personality_disorder")
which_sex <- c("Men","Women","Both")

min_age <- 15
max_age <- 85

load(file=paste(filepath_read,"/long_table.RData",sep=''))
data_surv <- copy(BAS_MHD_new)

# counting process format, one line per patient, on 'age' scale
data_surv[,`:=`(start=as.numeric(start-birth_d)/365.25,
                stop=as.numeric(end-birth_d)/365.25)]

# left-truncataion/right-censoring (age)
data_surv <- data.table(survSplit(Surv(start,stop,event)~.,data=data_surv,cut=c(min_age,max_age),episode="i"))
data_surv <- data_surv[i==2]
data_surv[,i:=NULL]
print(paste0("After excluding patients with no follow-up between ages 15 and 85: ",uniqueN(data_surv,"patient")))

# left-truncataion(age)
# data_surv <- data.table(survSplit(Surv(start,stop,event)~.,data=data_surv,cut=min_age,episode="i"))
# data_surv <- data_surv[i==2]
# data_surv[,i:=NULL]
# print(paste0("*After excluding patients with no follow-up aged ",min_age," or above: ",uniqueN(data_surv,"patient")))

# cleaning up AFTER censoring
data_surv[event==0,cod2:="Alive"]

# removing patients with unknown cause of death
#data_surv <- data_surv[cod2%in%c("Alive","Natural","Unnatural")]
#print(paste0("After excluding patients with deaths of unknown cause: ",uniqueN(data_surv,"patient")))

# formatting
data_surv[,`:=`(age_base=start,
                age_base_cat=cut(start,breaks=c(15,20,30,40,50,60,70,85,Inf),right=FALSE),
                calyear_base_cat=cut(year(birth_d+start*365.25),breaks=c(2010,2014,2017,2020,Inf),right=FALSE),
                sex=factor(sex,levels=c("Male","Female")),
                cod2=factor(cod2,levels=c("Natural","Unnatural","Unknown","Alive")))]

# looping over types of mental disorders
MHD_count_df <- data.frame(NULL)

N_total <- data_surv[,.N]
N_men <- data_surv[sex=="Male",.N]
N_women <- data_surv[sex=="Female",.N]

for(v in which_mhd)
{
  DT <- copy(data_surv)
  
  # age at mhd treatment
  DT[,`:=`(age_treat=as.numeric(get(v)-birth_d)/365.25,py=stop-start)]
  DT[age_treat>=stop,age_treat:=NA]                              # excluding disorders occuring post-censoring  
  DT[,exp:=ifelse(is.na(age_treat),"no MHD","MHD")]
  
  if(v=="first_any_admission")   # creating single table with/without mental illness
  {
    overall_df <- CreateTableOne(vars = c("py","age_base_cat","age_base","sex","calyear_base_cat"),strata="exp",data=DT,test=FALSE,addOverall=TRUE)
    overall_df <- print(overall_df,nonnormal=c("age_base","py"),showAllLevels=TRUE)
    overall_df <- data.table(data.frame(cbind(row.names(overall_df),overall_df)))
    deaths_df <- CreateTableOne(vars = "cod2",strata="exp",data=DT[cod2!="Alive"],test=FALSE,addOverall=TRUE)
    deaths_df <- print(deaths_df,showAllLevels=TRUE)
    deaths_df <- deaths_df[-1,]
    deaths_df <- data.table(data.frame(cbind(row.names(deaths_df),deaths_df)))
    overall_df <- rbind(overall_df,deaths_df)
    overall_df[,`:=`(Overall=gsub(")","%)",Overall),
                     MHD=gsub(")","%)",MHD),
                     no.MHD=gsub(")","%)",no.MHD))]
    overall_df[,`:=`(Overall=gsub("\\( ","\\(",Overall),
                     MHD=gsub("\\( ","\\(",MHD),
                     no.MHD=gsub("\\( ","\\(",no.MHD))]
    overall_df[,`:=`(Overall=gsub("\\.","·",Overall),
                     MHD=gsub("\\.","·",MHD),
                     no.MHD=gsub("\\.","·",no.MHD))]
    overall_df <- overall_df[,c(1,2,4,5,3)]
    write_xlsx(overall_df,path=file.path(filepath_write,"descriptive_overall_hos.xlsx"))
  }
  N_exp_total <- DT[exp=="MHD",.N]
  N_exp_men <- DT[exp=="MHD" & sex=="Male",.N]
  N_exp_women <- DT[exp=="MHD" & sex=="Female",.N]
  MHD_count_df <- rbind(MHD_count_df,data.frame(exp=v,
                                                men=paste0(N_exp_men," (",format(round(100*N_exp_men/N_men,digits=1),nsmall=1),"%)"),
                                                women=paste0(N_exp_women," (",format(round(100*N_exp_women/N_women,digits=1),nsmall=1),"%)"),
                                                total=paste0(N_exp_total," (",format(round(100*N_exp_total/N_total,digits=1),nsmall=1),"%)")))
  rm(DT,N_exp_total,N_exp_men,N_exp_women)
}        

MHD_count_df <- data.table(MHD_count_df)
MHD_count_df[,`:=`(men=gsub("\\.","·",men),
                   women=gsub("\\.","·",women),
                   total=gsub("\\.","·",total))]
write_xlsx(MHD_count_df,path=file.path(filepath_write,"descriptive_by_MHD_and_sex_hos.xlsx"))

toc()