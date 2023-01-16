########################################################################
#Data management for LYL and survival analysis                         #              
#Adds variables for number of each type of diagnosis                   #
#Any admission source                                                  #
#Date 21.06.2022                                                       #
#Anja Wettstein and Yann Ruffieux                                      #
#Runtime: ~5 minutes                                                   #
########################################################################

library(foreign)
library(Hmisc)
library(haven)
library(tictoc)
library(plyr)
library(dplyr)
library(data.table)
library(tidyverse)
library(stringr)

library(WriteXLS)

library(survival)

# read in data
filepath_read <- "C:/ISPM/Data/AfA-AfAc/Stata/AfAc_20211020"
filepath_write <- "C:/ISPM/Data/HIV-mental disorders/AfAc_excess_mortality/processed"

which_mhd<- c("organic","substance_use_disorder","alcohol_use_disorder",
              "drug_use_disorder","psychotic","mood_disorder","bipolar",
              "depression","anxiety","generalised_anxiety_disorder","PTSD",
              "asc_developmental_disorders","asc_behavioural_disorders","asc_physical_factors",
              "eating_disorder","asc_personality_disorder")

tic()
  
ICD10_F <-read_dta(file.path(filepath_read,"ICD10_F.dta"), encoding=NULL)
ICD10_R <-read_dta(file.path(filepath_read,"ICD10_R.dta"), encoding=NULL)
BAS<-read_dta(file.path(filepath_read,"BAS.dta"), encoding=NULL)
FUP <- read_dta(file.path(filepath_read,"FUP.dta"), encoding=NULL)
FUPwide<-read_dta(file.path(filepath_read,"FUPwide.dta"), encoding=NULL)
VITAL<-read_dta(file.path(filepath_read,"VITAL.dta"), encoding=NULL)

#To import categorical variables from STATA correcctly

(ICD10_F <- mutate_if(ICD10_F,is.labelled,as_factor))
(ICD10_R <- mutate_if(ICD10_R,is.labelled,as_factor))
(BAS <- mutate_if(BAS,is.labelled,as_factor))
(FUP <- mutate_if(FUP,is.labelled,as_factor))
(FUPwide <- mutate_if(FUPwide,is.labelled,as_factor))
(VITAL <- mutate_if(VITAL,is.labelled,as_factor))

#Set path for saving data frames

#Clean data

#1) Clean BAS

BAS<- select(BAS,c("patient", "birth_d","birth_d_a","sex","popgrp"))

print(paste0("Starting number of patients: ",uniqueN(BAS,"patient")))

# 2) Clean ICD10

#Join ICD10_f and ICD10_R DF

ICD10<-bind_rows(ICD10_F, ICD10_R[substr(ICD10_R$icd10_code,1,5)%in% paste0("R44.",0:3),], .id=NULL)


#Subset ICD10

ICD10_ss<-subset(ICD10,select= c(patient,icd10_date,icd10_code,source, discharge_date, icd_1, icd_23))

#Create disorder categories
ICD10_ss$disorder_class_any<-1
ICD10_ss$disorder_class_organic<-as.numeric(substr(ICD10_ss$icd10_code,1,2)=="F0")                                                                         #Organic mental disorder
#ICD10_ss$disorder_class_substance_use_disorder<-as.numeric(substr(ICD10_ss$icd10_code,1,2)=="F1")                                                          #Substance use disorder                                        
ICD10_ss$disorder_class_substance_use_disorder<-as.numeric(substr(ICD10_ss$icd10_code,1,3)%in%c("F10","F11","F12","F13","F14","F15","F16","F18","F19"))     #Substance use disorder                                        
ICD10_ss$disorder_class_alcohol_use_disorder<-as.numeric(substr(ICD10_ss$icd10_code,1,3)=="F10")                                                            #Alcohol use disorder
#ICD10_ss$disorder_class_drug_use_disorder<-as.numeric(substr(ICD10_ss$icd10_code,1,3)%in%c("F11","F12","F13", "F14","F15","F16","F18", "F19"))             #Drug use disorder
ICD10_ss$disorder_class_drug_use_disorder<-as.numeric(substr(ICD10_ss$icd10_code,1,3) %in% paste0("F",11:16)|substr(ICD10_ss$icd10_code,1,3)%in% paste0("F",18:19)) #Drug use disorder
ICD10_ss$disorder_class_psychotic<-as.numeric(substr(ICD10_ss$icd10_code,1,2)=="F2"|substr(ICD10_ss$icd10_code,1,5)%in% paste0("R44.",0:3))       #Psychotic disorder
ICD10_ss$disorder_class_mood_disorder<-as.numeric(substr(ICD10_ss$icd10_code,1,2)=="F3")                                                                   #Mood disorder
ICD10_ss$disorder_class_bipolar<-as.numeric(substr(ICD10_ss$icd10_code,1,3)=="F31")                                                                        #Bipolar
ICD10_ss$disorder_class_depression<-as.numeric(substr(ICD10_ss$icd10_code,1,3)%in%c("F32","F33")| substr(ICD10_ss$icd10_code,1,5)=="F34.1")                #Depression and Dysthymia
#ICD10_ss$disorder_class_anxiety<-as.numeric(substr(ICD10_ss$icd10_code,1,3)%in%c("F40","F41","F42","F43","F44","F45","F46","F47","F48"))                   #Anxiety
ICD10_ss$disorder_class_anxiety<-as.numeric(substr(ICD10_ss$icd10_code,1,3) %in% paste0("F",40:49))                                                        #Anxiety
ICD10_ss$disorder_class_generalised_anxiety_disorder<-as.numeric(substr(ICD10_ss$icd10_code,1,5)=="F41.1")                                                 #Generalised anxiety disorder
ICD10_ss$disorder_class_PTSD<-as.numeric(substr(ICD10_ss$icd10_code,1,5)=="F43.1")                                                                         #Post-traumatic stress disorder
#ICD10_ss$disorder_class_asc_physical_factors<-as.numeric(substr(ICD10_ss$icd10_code,1,3)%in%c("F50","F51","F52", "F53","F54","F55","F56","F57","F58","F59"))#Behavioural syndromes associated with physical factors    
ICD10_ss$disorder_class_asc_physical_factors<-as.numeric(substr(ICD10_ss$icd10_code,1,3)%in%paste0("F",50:59))                                              #Behavioural syndromes associated with physical factors    
ICD10_ss$disorder_class_eating_disorder<-as.numeric(substr(ICD10_ss$icd10_code,1,3)=="F50")                                                                 #Eating Disorder 
ICD10_ss$disorder_class_non_organic_sleep_disorder<-as.numeric(substr(ICD10_ss$icd10_code,1,3)=="F51")                                                      #Non-organic sleep disorders         #
ICD10_ss$disorder_class_asc_personality_disorder<-as.numeric(substr(ICD10_ss$icd10_code,1,2)=="F6")                                                         #Personality disorder
ICD10_ss$disorder_class_asc_developmental_disorders<-as.numeric(substr(ICD10_ss$icd10_code,1,2)=="F8")                                                      #Developmental disorders
#ICD10_ss$disorder_class_asc_behavioural_disorders<-as.numeric(substr(ICD10_ss$icd10_code,1,2)=="F9")                                                       #Behavioural disorders
ICD10_ss$disorder_class_asc_behavioural_disorders<-as.numeric(substr(ICD10_ss$icd10_code,1,3) %in% paste0("F",90:98))                                       #Behavioural disorders

#Build new categories for source:

ICD10_ss$source_any<-as.numeric(substr(ICD10_ss$source,1,3)%in%c("HOS","OPD","MED"))  #Any source; includes hospitalisation, OPD, and MED
ICD10_ss$source_hos<-as.numeric(substr(ICD10_ss$source,1,3)=="HOS")             #Only hospitalisation as source

ICD10_ss <- data.table(ICD10_ss)

# couting number of diagnoses for each type of disorder
ICD10_ss[,N_any_admission:=sum(disorder_class_any),by="patient"]
for(v in which_mhd)
{
  in_var_name <- paste0("disorder_class_",v)
  new_var_name <- paste0("N_any_admission_",v)
  ICD10_ss[,eval(new_var_name):=sum(get(in_var_name)),by="patient"]
}

# B009290634 - should have no 'first behavioural ASC' since there's only one diagnosis

#i) Create a new table which will have first occurence of disorders for each patient

setorder(ICD10_ss,"patient","icd10_date")
BAS_MHD <- data.table(BAS)          

# date of first any admission for any mental disroder, from any source
BAS_MHD <- ICD10_ss[N_any_admission>0,.(patient,first_any_admission=icd10_date)][BAS_MHD,on="patient",mult="first"]

# date of first admission, for each type of disorder
for(v in which_mhd)
{
  dis <- paste0("disorder_class_",v)
  N_dis <- paste0("N_any_admission_",v)
  BAS_MHD <- ICD10_ss[get(dis)==1 & get(N_dis)>0,.(patient,dummy_name=icd10_date)][BAS_MHD,on="patient",mult="first"]
  setnames(BAS_MHD,"dummy_name",paste0("first_any_admission_",v))
}

#Merge with FUPwide

BAS_MHD_new <- merge(BAS_MHD,FUPwide,by="patient" )

print(paste0("*After excluding patients with no follow-up between Jan 2011 and Jul 2020: ",uniqueN(BAS_MHD_new,"patient")))

#Subset VITAL data frame: exclude all patients who weren't linked to NPR (linked==0 in VITAL) 

VITAL_linked<-subset(VITAL,VITAL$linked==1)

#Merge the VITAL_linked data frame to BAS_MDH_new but only keep death_y, death_d, cod1, and cod2

VITAL_linked_to_merge<-VITAL_linked[,c(1,9,10,11,12)]
BAS_MHD_new <- merge(BAS_MHD_new,VITAL_linked_to_merge,by="patient" )
print(paste0("*After excluding patients not linked to NPR: ",uniqueN(BAS_MHD_new,"patient")))

# Removing patients with missing sex
BAS_MHD_new <- BAS_MHD_new[sex=="Male" | sex=="Female"]
print(paste0("*After excluding patients with unknown sex: ",uniqueN(BAS_MHD_new,"patient")))

# removing patients with missing date of birth
BAS_MHD_new <- BAS_MHD_new[!is.na(birth_d)]
print(paste0("*After excluding patients with missing date of birth: ",uniqueN(BAS_MHD_new,"patient")))

#Censor deaths that occur after a patient leaves the study

BAS_MHD_new$event <- BAS_MHD_new$death_y

BAS_MHD_new$event[BAS_MHD_new$death_y==1&BAS_MHD_new$death_d>BAS_MHD_new$end] <- 0

# appending enrolment date into the program (from FUP table, not FUPwide!)
FUP <- data.table(FUP)
setorder(FUP,"patient","start")
BAS_MHD_new <- FUP[,.(patient,enrol_d=start)][BAS_MHD_new,on="patient",mult="first"]
# some patients are enrolled before birth -> switching enrolment date to the date of birth
BAS_MHD_new[enrol_d<birth_d,enrol_d:=birth_d]

# appending number of diagnoses of each type
var_names <- c("patient","N_any_admission",paste0("N_any_admission_",which_mhd))
BAS_MHD_new <- ICD10_ss[,..var_names][BAS_MHD_new,on="patient",mult="first"]

# fill in entries for patients not appearing in ICD10 table
BAS_MHD_new[is.na(N_any_admission),N_any_admission:=0]
for(v in which_mhd)
{
  N_dis <- paste0("N_any_admission_",v)
  BAS_MHD_new[is.na(get(N_dis)),eval(N_dis):=0]
}

# Save long table

save(BAS_MHD_new,file=file.path(filepath_write,"long_table_Ndiag_any_admission.RData"))

toc()