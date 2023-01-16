##########################
#ELYL Table any Admission#
#Anja Wettstein          #
# 10/08/2022             #
##########################

library(data.table)
library(dplyr)
library(writexl)
library(tidyverse)

filepath_lookup <- "C:/ISPM/HomeDir/HIV-mental disorders/AfAc_excess_mortality/LYL/lookup"
filepath_est <- "C:/ISPM/Data/HIV-mental disorders/AfAc_excess_mortality/lillies/point_estimates"
filepath_boot <- "C:/ISPM/Data/HIV-mental disorders/AfAc_excess_mortality/lillies/bootstrap/joined"
filepath_write <- "C:/ISPM/HomeDir/HIV-mental disorders/AfAc_excess_mortality/output/tables/model output"

# ordered exposures, nested disorders at the end
target <- c("Any Mental Health Diagnosis", "Organic Disorders","Substance Use Disorders","Psychotic Disorders","Mood Disorders","Anxiety Disorders",
            "Developmental Disorders","Eating Disorders","Personality Disorders","Alcohol Use Disorders","Drug Use Disorders",
            "Bipolar Disorders","Depressive Disorders","Generalised Anxiety Disorders","PTSD")

load(file.path(filepath_est,"elyl_any_admission.RData"))
estimates <- data.table(LYL_diff_df)
rm(LYL_diff_df)
load(file.path(filepath_boot,"CIs_any_admission.RData"))
CI <- data.table(CIs_any_admission)
rm(CIs_any_admission)

expo <- read.csv(file.path(filepath_lookup,"exposures.csv"), header=TRUE, sep=",")
gender_lookup<- read.csv(file.path(filepath_lookup,"sex.csv"), header=TRUE, sep=",")

gender_lookup <- rbind(gender_lookup,c("Both","Both"))

names(CI)[names(CI)=="exposure"] <- "type"
names(estimates)[names(estimates)=="sex"] <- "gender"

elyl_all_a <- merge(estimates,expo, by="type")

elyl_all_b <- merge(CI,expo, by="type")
elyl_all_b <- merge(elyl_all_b,gender_lookup, by="sex")
elyl_all_b$sex <- NULL

elyl_all <- merge(elyl_all_a,elyl_all_b, by =c("type","gender","d_lables"))

elyl_all <- data.table(elyl_all)
elyl_all <- elyl_all[gender!="Both"]
elyl_all[,`:=`(d_lables=factor(d_lables,levels=target),gender=factor(gender,levels=c("Men","Women")))]
elyl_all <- elyl_all[!is.na(d_lables)]
setorder(elyl_all,d_lables,gender)

mod_cols <- c("TotalLYL","Natural","Unnatural","Unknown","lcl_allcause","ucl_allcause","lcl_natural","ucl_natural",
              "lcl_unnatural","ucl_unnatural","lcl_uk","ucl_uk")
elyl_all[ , (mod_cols) := lapply(.SD,round,digits=2), .SDcols = mod_cols] 

elyl_all[,`:=`(All_cause=paste0(format(TotalLYL,nsmall=2)," [",format(lcl_allcause,nsmall=2),"; ",format(ucl_allcause,nsmall=2),"]"),
               Natural=paste0(format(Natural,nsmall=2)," [",format(lcl_natural,nsmall=2),"; ",format(ucl_natural,nsmall=2),"]"),
               Unnatural=paste0(format(Unnatural,nsmall=2)," [",format(lcl_unnatural,nsmall=2),"; ",format(ucl_unnatural,nsmall=2),"]"),
               Unknown=paste0(format(Unknown,nsmall=2)," [",format(lcl_uk,nsmall=2),"; ",format(ucl_uk,nsmall=2),"]"))]

elyl_all <- elyl_all[,.(d_lables,gender,All_cause,Natural,Unnatural,Unknown)]

write_xlsx(x = elyl_all, path=file.path(filepath_write,"eLYL_any_admission_table.xlsx"), col_names = TRUE)
