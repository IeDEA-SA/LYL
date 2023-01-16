#########################################################################
#Data management for LYL and survival analysis                          #
#Not censoring deaths occurring less than one year after transfer out   #
#Date 21.06.2022                                                        #
#Anja Wettstein and Yann Ruffieux                                       #
#Runtime: ~5 minutes                                                    #
#########################################################################

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

tic()
  
ICD10_F <-read_dta(file.path(filepath_read,"ICD10_F.dta"), encoding=NULL)
ICD10_R <-read_dta(file.path(filepath_read,"ICD10_R.dta"), encoding=NULL)
MEDS<-read_dta(file.path(filepath_read,"MED_ATC_N.dta"), encoding=NULL)
BAS<-read_dta(file.path(filepath_read,"BAS.dta"), encoding=NULL)
FUP <- read_dta(file.path(filepath_read,"FUP.dta"), encoding=NULL)
FUPwide<-read_dta(file.path(filepath_read,"FUPwide.dta"), encoding=NULL)
VITAL<-read_dta(file.path(filepath_read,"VITAL.dta"), encoding=NULL)

#To import categorical variables from STATA correcctly

(ICD10_F <- mutate_if(ICD10_F,is.labelled,as_factor))
(ICD10_R <- mutate_if(ICD10_R,is.labelled,as_factor))
(MEDS <- mutate_if(MEDS,is.labelled,as_factor))
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

#Select relevant columns
#ICD10_new0<-ICD10_ss #Any mental disorder
#ICD10_new1<-ICD10_ss[ICD10_ss$icd10_code %like% "F0[0-9]",] #Organic mental disorders
#ICD10_new2<-ICD10_ss[ICD10_ss$icd10_code %like% "F1[0-9]",] #Substance use disorders
#ICD10_new3<-ICD10_ss[ICD10_ss$icd10_code %like% "F2[0-9]",] #Schizophrenia, schizotypal and delusional disorders
#ICD10_new4<-ICD10_ss[ICD10_ss$icd10_code %like% "F3[0-9]",] #Mood affective disorders
#ICD10_new5<-ICD10_ss[ICD10_ss$icd10_code %like% "F4[0-8]",] #Neurotic, stress-related and somatoform disorders (Anxiety)
#ICD10_new6<-ICD10_ss[ICD10_ss$icd10_code %like% "F5[0-9]",] #Behavioural syndromes associated with physical factors
#ICD10_new7<-ICD10_ss[ICD10_ss$icd10_code %like% "F6[0-9]",] #Personality disorders
#ICD10_new8<-ICD10_ss[ICD10_ss$icd10_code %like% "F8[0-9]",] #Developmental disorders
#ICD10_new9<-ICD10_ss[ICD10_ss$icd10_code %like% "F9[0-8]",] #Behavioural disorders


#Bind rows

#ICD10_new<-bind_rows(ICD10_new0,ICD10_new1,ICD10_new2,ICD10_new3,ICD10_new4,ICD10_new5,ICD10_new6,ICD10_new7,ICD10_new8,ICD10_new9, .id=NULL)

#Save this df

#  save(ICD10_new,file=file.path(filepath_write,"ICD10_new_new_categories.csv"))

#Create disorder categories
ICD10_ss$disorder_class_any<-"any"
ICD10_ss$disorder_class_any<-1
ICD10_ss$disorder_class_organic<-as.numeric(substr(ICD10_ss$icd10_code,1,2)=="F0")                                                                         #Organic mental disorder
#ICD10_ss$disorder_class_substance_use_disorder<-as.numeric(substr(ICD10_ss$icd10_code,1,2)=="F1")                                                          #Substance use disorder                                        
ICD10_ss$disorder_class_substance_use_disorder<-as.numeric(substr(ICD10_ss$icd10_code,1,3)%in%c("F10","F11","F12","F13","F14","F15","F16","F18","F19"))     #Substance use disorder                                        
ICD10_ss$disorder_class_alcohol_use_disorder<-as.numeric(substr(ICD10_ss$icd10_code,1,3)=="F10")                                                            #Alcohol use disorder
#ICD10_ss$disorder_class_drug_use_disorder<-as.numeric(substr(ICD10_ss$icd10_code,1,3)%in%c("F11","F12","F13", "F14","F15","F16","F18", "F19"))             #Drug use disorder
ICD10_ss$disorder_class_drug_use_disorder<-as.numeric(substr(ICD10_ss$icd10_code,1,3) %in% paste0("F",11:16)|substr(ICD10_ss$icd10_code,1,3)%in% paste0("F",18:19)) #Drug use disorder
ICD10_ss$disorder_class_psychotic_disorder<-as.numeric(substr(ICD10_ss$icd10_code,1,2)=="F2"|substr(ICD10_ss$icd10_code,1,5)%in% paste0("R44.",0:3))       #Psychotic disorder
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

#sapply(ICD10_ss,class)

#Build new categories for source:

ICD10_ss$source_any<-as.numeric(substr(ICD10_ss$source,1,3)%in%c("HOS","OPD","MED"))  #Any source; includes hospitalisation, OPD, and MED
ICD10_ss$source_hos<-as.numeric(substr(ICD10_ss$source,1,3)=="HOS")             #Only hospitalisation as source

#i) Create a new table which will have first occurence of disorders (per type/source) for each patient

ICD10_new <- data.table(ICD10_ss)
setorder(ICD10_new,"patient","icd10_date")
BAS_MHD <- data.table(BAS)          

# date of first any admission for any mental
BAS_MHD <- ICD10_new[source_any==1&disorder_class_any==1,.(patient,first_any_admission=icd10_date)][BAS_MHD,on="patient",mult="first"]
# date of first hospitalization for ANY mental disorder
BAS_MHD <- ICD10_new[source_hos==1&disorder_class_any==1,.(patient,first_hos_any=icd10_date)][BAS_MHD,on="patient",mult="first"]

# dates of first any admission for organic disorder
BAS_MHD <- ICD10_new[source_any==1 & disorder_class_organic==1,.(patient,first_any_admission_organic=icd10_date)][BAS_MHD,on="patient",mult="first"]
# date of first hospitalization for organic disorder
BAS_MHD <- ICD10_new[source_hos==1 & disorder_class_organic==1,.(patient,first_hos_organic=icd10_date)][BAS_MHD,on="patient",mult="first"]

# dates of first any admission for substance use disorder
BAS_MHD <- ICD10_new[source_any==1 & disorder_class_substance_use_disorder==1,.(patient,first_any_admission_substance_use_disorder=icd10_date)][BAS_MHD,on="patient",mult="first"]
# date of first hospitalization for substance use disorder
BAS_MHD <- ICD10_new[source_hos==1 & disorder_class_substance_use_disorder==1,.(patient,first_hos_substance_use_disorder=icd10_date)][BAS_MHD,on="patient",mult="first"]


#dates of first any admission for alcohol use disorder
BAS_MHD <- ICD10_new[source_any==1 & disorder_class_alcohol_use_disorder==1,.(patient,first_any_admission_alcohol_use_disorder=icd10_date)][BAS_MHD,on="patient",mult="first"]
# date of first hospitalization for alcohol use disorder
BAS_MHD <- ICD10_new[source_hos==1 & disorder_class_alcohol_use_disorder==1,.(patient,first_hos_alcohol_use_disorder=icd10_date)][BAS_MHD,on="patient",mult="first"]

#dates of first any admission for drug use disorder
BAS_MHD <- ICD10_new[source_any==1 & disorder_class_drug_use_disorder==1,.(patient,first_any_admission_drug_use_disorder=icd10_date)][BAS_MHD,on="patient",mult="first"]
# date of first hospitalization for drug use disorder
BAS_MHD <- ICD10_new[source_hos==1 & disorder_class_drug_use_disorder==1,.(patient,first_hos_drug_use_disorder=icd10_date)][BAS_MHD,on="patient",mult="first"]

#dates of first any admission for psychotic disorder
BAS_MHD <- ICD10_new[source_any==1 & disorder_class_psychotic_disorder==1,.(patient,first_any_admission_psychotic=icd10_date)][BAS_MHD,on="patient",mult="first"]
# date of first hospitalization for psychotic disorder
BAS_MHD <- ICD10_new[source_hos==1 & disorder_class_psychotic_disorder==1,.(patient,first_hos_psychotic=icd10_date)][BAS_MHD,on="patient",mult="first"]

#dates of first any admission for mood disorder
BAS_MHD <- ICD10_new[source_any==1 & disorder_class_mood_disorder==1,.(patient,first_any_admission_mood_disorder=icd10_date)][BAS_MHD,on="patient",mult="first"]
# date of first hospitalization for mood disorder
BAS_MHD <- ICD10_new[source_hos==1 & disorder_class_mood_disorder==1,.(patient,first_hos_mood_disorder=icd10_date)][BAS_MHD,on="patient",mult="first"]

#dates of first any admission for bipolar
BAS_MHD <- ICD10_new[source_any==1 & disorder_class_bipolar==1,.(patient,first_any_admission_bipolar=icd10_date)][BAS_MHD,on="patient",mult="first"]
# date of first hospitalization for bipolar
BAS_MHD <- ICD10_new[source_hos==1 & disorder_class_bipolar==1,.(patient,first_hos_bipolar=icd10_date)][BAS_MHD,on="patient",mult="first"]

#dates of first any admission for depression
BAS_MHD <- ICD10_new[source_any==1 & disorder_class_depression==1,.(patient,first_any_admission_depression=icd10_date)][BAS_MHD,on="patient",mult="first"]
# date of first hospitalization for depression
BAS_MHD <- ICD10_new[source_hos==1 & disorder_class_depression==1,.(patient,first_hos_depression=icd10_date)][BAS_MHD,on="patient",mult="first"]

#dates of first any admission for anxiety
BAS_MHD <- ICD10_new[source_any==1 & disorder_class_anxiety==1,.(patient,first_any_admission_anxiety=icd10_date)][BAS_MHD,on="patient",mult="first"]
# date of first hospitalization for anxiety
BAS_MHD <- ICD10_new[source_hos==1 & disorder_class_anxiety==1,.(patient,first_hos_anxiety=icd10_date)][BAS_MHD,on="patient",mult="first"]

#dates of first any admission for generalised anxiety disorder
BAS_MHD <- ICD10_new[source_any==1 & disorder_class_generalised_anxiety_disorder==1,.(patient,first_any_admission_generalised_anxiety_disorder=icd10_date)][BAS_MHD,on="patient",mult="first"]
# date of first hospitalization for generalised anxiety disorder
BAS_MHD <- ICD10_new[source_hos==1 & disorder_class_generalised_anxiety_disorder==1,.(patient,first_hos_generalised_anxiety_disorder=icd10_date)][BAS_MHD,on="patient",mult="first"]

#dates of first any admission for PTSD
BAS_MHD <- ICD10_new[source_any==1 & disorder_class_PTSD==1,.(patient,first_any_admission_PTSD=icd10_date)][BAS_MHD,on="patient",mult="first"]
# date of first hospitalization for PTSD
BAS_MHD <- ICD10_new[source_hos==1 & disorder_class_PTSD==1,.(patient,first_hos_PTSD=icd10_date)][BAS_MHD,on="patient",mult="first"]

#dates of first any admission for associated physical factors
BAS_MHD <- ICD10_new[source_any==1 & disorder_class_asc_physical_factors==1,.(patient,first_any_admission_asc_physical_factors=icd10_date)][BAS_MHD,on="patient",mult="first"]
# date of first hospitalization for  associated physical factors
BAS_MHD <- ICD10_new[source_hos==1 & disorder_class_asc_physical_factors==1,.(patient,first_hos_asc_physical_factors=icd10_date)][BAS_MHD,on="patient",mult="first"]

#dates of first any admission for eating disorder
BAS_MHD <- ICD10_new[source_any==1 & disorder_class_eating_disorder==1,.(patient,first_any_admission_eating_disorder=icd10_date)][BAS_MHD,on="patient",mult="first"]
# date of first hospitalization for eating isorder
BAS_MHD <- ICD10_new[source_hos==1 & disorder_class_eating_disorder==1,.(patient,first_hos_eating_disorder=icd10_date)][BAS_MHD,on="patient",mult="first"]

#dates of first any admission for non-organic sleep disorders
BAS_MHD <- ICD10_new[source_any==1 & disorder_class_non_organic_sleep_disorder==1,.(patient,first_any_admission_non_organic_sleep_disorder=icd10_date)][BAS_MHD,on="patient",mult="first"]
# date of first hospitalization for eating isorder
BAS_MHD <- ICD10_new[source_hos==1 & disorder_class_non_organic_sleep_disorder==1,.(patient,first_hos_non_organic_sleep_disorder=icd10_date)][BAS_MHD,on="patient",mult="first"]

#dates of first any admission for  associated personality
BAS_MHD <- ICD10_new[source_any==1 & disorder_class_asc_personality_disorder==1,.(patient,first_any_admission_asc_personality_disorder=icd10_date)][BAS_MHD,on="patient",mult="first"]
# date of first hospitalization for  associated asc personality disorder
BAS_MHD <- ICD10_new[source_hos==1 & disorder_class_asc_personality_disorder==1,.(patient,first_hos_asc_personality_disorder=icd10_date)][BAS_MHD,on="patient",mult="first"]

#dates of first any admission for associated developmental disorder
BAS_MHD <- ICD10_new[source_any==1 & disorder_class_asc_behavioural_disorders==1,.(patient,first_any_admission_asc_behavioural_disorders=icd10_date)][BAS_MHD,on="patient",mult="first"]
# date of first hospitalization for associated developmental
BAS_MHD <- ICD10_new[source_hos==1 & disorder_class_asc_behavioural_disorders==1,.(patient,first_hos_asc_behavioural_disorders=icd10_date)][BAS_MHD,on="patient",mult="first"]

#dates of first any admission for associated behavioural disorder
BAS_MHD <- ICD10_new[source_any==1 & disorder_class_asc_developmental_disorders==1,.(patient,first_any_admission_asc_developmental_disorders=icd10_date)][BAS_MHD,on="patient",mult="first"]
# date of first hospitalization for  associated behavioural developmental
BAS_MHD <- ICD10_new[source_hos==1 & disorder_class_asc_developmental_disorders==1,.(patient,first_hos_asc_developmental_disorders=icd10_date)][BAS_MHD,on="patient",mult="first"]


# 3) Clean MEDS

#Subset MEDS

MEDS_ss<-subset(MEDS,select= c(patient,med_sd,med_id,quantity))

#Select relevant columns

MEDS_new1<-MEDS_ss[MEDS_ss$med_id %like% "N0[5-6]",]
MEDS_new2<-MEDS_ss[MEDS_ss$med_id %like% "N07B",]

#Bind rows

MEDS_new<-bind_rows(MEDS_new1,MEDS_new2,.id=NULL)

#Save this df

save(MEDS_new,file=file.path(filepath_write,"MEDS_new.RData"))

#Create disorder categories

MEDS_new$med_class_any<-as.numeric(substr(MEDS_new$med_id,1,4)=="N06A"|substr(MEDS_new$med_id,1,4)%in% c("N05A","N05B"))
MEDS_new$med_class_antidepressant<-as.numeric(substr(MEDS_new$med_id,1,4)=="N06A")
MEDS_new$med_class_antipsychotic<-as.numeric(substr(MEDS_new$med_id,1,4)=="N05A")
MEDS_new$med_class_anxiolytic<-as.numeric(substr(MEDS_new$med_id,1,4)=="N05B")
MEDS_new$med_class_substance_use<-as.numeric(substr(MEDS_new$med_id,1,4)=="N07B")

# Create a new table which will have first adminstration of specific medication per type for each patient
MEDS_new <- data.table(MEDS_new)
setorder(MEDS_new,"patient","med_sd")

# dates of first any med
BAS_MHD <- MEDS_new[med_class_any==1,.(patient,first_any_med=med_sd)][BAS_MHD,on="patient",mult="first"]

# dates of first antidepressant
BAS_MHD <- MEDS_new[med_class_antidepressant==1,.(patient,first_antidepressant=med_sd)][BAS_MHD,on="patient",mult="first"]

# dates of first antipsychotic
BAS_MHD <- MEDS_new[med_class_antipsychotic==1,.(patient,first_antipsychotic=med_sd)][BAS_MHD,on="patient",mult="first"]

# dates of first anxiolytic
BAS_MHD <- MEDS_new[med_class_anxiolytic==1,.(patient,first_anxiolytic=med_sd)][BAS_MHD,on="patient",mult="first"]

# dates of first stubstance use med
BAS_MHD <- MEDS_new[med_class_substance_use==1,.(patient,first_substance_use_med=med_sd)][BAS_MHD,on="patient",mult="first"]


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

#Censor deaths that occur more than one year after a patient leaves the study, or after database closure
BAS_MHD_new$event <- BAS_MHD_new$death_y
BAS_MHD_new[death_y==1 & (death_d>end+365.25 | death_d>=as.Date("2020-07-01")),event:=0]
BAS_MHD_new[event==1,end:=pmax(end,death_d,na.rm=TRUE)]

# appending enrolment date into the program (from FUP table, not FUPwide!)
FUP <- data.table(FUP)
setorder(FUP,"patient","start")
BAS_MHD_new <- FUP[,.(patient,enrol_d=start)][BAS_MHD_new,on="patient",mult="first"]
# some patients are enrolled before birth -> switching enrolment date to the date of birth
BAS_MHD_new[enrol_d<birth_d,enrol_d:=birth_d]

# Save long table

save(BAS_MHD_new,file=file.path(filepath_write,"long_table_alt_censoring.RData"))

toc()