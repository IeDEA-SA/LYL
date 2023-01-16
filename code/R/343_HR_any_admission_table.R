# produces tables of hazard ratios, any admission type

library(data.table)
library(writexl)
library(readxl)

filepath_lookup <- "C:/ISPM/HomeDir/HIV-mental disorders/Anja/lookup"
filepath_read <- "C:/ISPM/HomeDir/HIV-mental disorders/Anja/output/tables/model output"
filepath_write <- "C:/ISPM/HomeDir/HIV-mental disorders/Anja/output/tables/model output"

#read in data
expo <- read.csv(file.path(filepath_lookup,"exposures.csv"),header=TRUE, sep=",")
mod1 <- data.table(read_xlsx(file.path(filepath_read,"aHRs_model1_any_admission.xlsx")))
mod2 <- data.table(read_xlsx(file.path(filepath_read,"aHRs_model2_any_admission.xlsx")))
mod3 <- data.table(read_xlsx(file.path(filepath_read,"aHRs_model3_any_admission.xlsx")))
mod4 <- data.table(read_xlsx(file.path(filepath_read,"aHRs_model4_any_admission.xlsx")))

names(expo)[1] <- "disorder"

# ordered exposures, nested disorders at the end
target <- c("Any Mental Health Diagnosis", "Organic Disorders","Substance Use Disorders","Psychotic Disorders","Mood Disorders","Anxiety Disorders",
            "Developmental Disorders","Personality Disorders","Alcohol Use Disorders","Drug Use Disorders",
            "Bipolar Disorders","Depressive Disorders","Generalised Anxiety Disorders","PTSD","Eating Disorders")

###### Models 1 and 3 in one table #########

mod1[cause=="all",cause:="All"]
mod1 <- merge(mod1,expo,by="disorder")
mod1[,ucl:=as.numeric(ucl)]
mod1[,inf_indicator:=0]
mod1[ucl>200,`:=`(inf_indicator=1,ucl=NA)]           # to remove cod/exposure/sex combinations with no event in exposure group
mod1_table <- mod1[,.(Disorder=d_lables,Cause=cause,
                          HR1=paste0(format(round(estimate,2),nsmall=2)," [",format(round(lcl,2),nsmall=2),"; ",format(round(ucl,2),nsmall=2),"]"),
                          p_val=as.character(format(round(p_int_sex_exp,digits=3)),nsmall=2),
                          p_ph_exp,inf_indicator)]
#mod1_table[p_ph_exp<0.05,HR1:=paste0(HR1,"*")]
mod1_table[p_val<0.001,p_val:="<0.001"]
mod1_table[inf_indicator==1,`:=`(HR1=".",p_val=".")]
mod1_table[,`:=`(p_ph_exp=NULL,inf_indicator=NULL)]

mod3[cause=="all",cause:="All"]
mod3 <- merge(mod3,expo,by="disorder")
mod3[,ucl:=as.numeric(ucl)]
mod3[,inf_indicator:=0]
mod3[ucl>200,`:=`(inf_indicator=1,ucl=NA)]           # to remove cod/exposure/sex combinations with no event in exposure group
mod3_table <- mod3[,.(Disorder=d_lables,Cause=cause,
                          HR3=paste0(format(round(estimate,2),nsmall=2)," [",format(round(lcl,2),nsmall=2),"; ",format(round(ucl,2),nsmall=2),"]"),
                          p_ph_exp,inf_indicator)]
#mod3_table[p_ph_exp<0.05,HR3:=paste0(HR3,"*")]
mod3_table[inf_indicator==1,HR3:="."]
mod3_table[,`:=`(p_ph_exp=NULL,inf_indicator=NULL)]

mod13_table <- merge(mod1_table,mod3_table,by=c("Disorder","Cause"),all.x=TRUE)
mod13_table[,`:=`(Disorder=factor(Disorder,levels=target),
                  Cause=factor(Cause,levels=c("All","Natural","Unnatural")))]
mod13_table[is.na(HR3),HR3:="-"]
setorder(mod13_table,"Disorder","Cause")
mod13_table <- mod13_table[,.(Disorder,Cause,HR1,HR3,p_val)]

rm(mod1_table,mod3_table)

###### Models 2 and 4 in one table #########

mod2[cause=="all",cause:="All"]
mod2 <- merge(mod2,expo,by="disorder")
mod2[,ucl:=as.numeric(ucl)]
mod2[,inf_indicator:=0]
mod2[ucl>200,`:=`(inf_indicator=1,ucl=NA)]           # to remove cod/exposure/sex combinations with no event in exposure group
mod2_table <- mod2[,.(Disorder=d_lables,Cause=cause,Sex=sex,
                          HR2=paste0(format(round(estimate,2),nsmall=2)," [",format(round(lcl,2),nsmall=2),"; ",format(round(ucl,2),nsmall=2),"]"),
                          p_ph_exp,inf_indicator)]
#mod2_table[p_ph_exp<0.05,HR2:=paste0(HR2,"*")]
mod2_table[inf_indicator==1,HR2:="."]
mod2_table[,`:=`(p_ph_exp=NULL,inf_indicator=NULL)]

mod4[cause=="all",cause:="All"]
mod4 <- merge(mod4,expo,by="disorder")
mod4[,inf_indicator:=0]
mod4[,ucl:=as.numeric(ucl)]
mod4[ucl>200,`:=`(inf_indicator=1,ucl=NA)]           # to remove cod/exposure/sex combinations with no event in exposure group
mod4_table <- mod4[,.(Disorder=d_lables,Cause=cause,Sex=sex,
                          HR4=paste0(format(round(estimate,2),nsmall=2)," [",format(round(lcl,2),nsmall=2),"; ",format(round(ucl,2),nsmall=2),"]"),
                          p_ph_exp,inf_indicator)]
#mod4_table[p_ph_exp<0.05,HR4:=paste0(HR4,"*")]
mod4_table[inf_indicator==1,HR4:="."]
mod4_table[,`:=`(p_ph_exp=NULL,inf_indicator=NULL)]

mod24_table <- merge(mod2_table,mod4_table,by=c("Disorder","Cause","Sex"),all.x=TRUE)
mod24_table[,`:=`(Disorder=factor(Disorder,levels=target),
                  Cause=factor(Cause,levels=c("All","Natural","Unnatural")),
                  Sex=factor(Sex,levels=c("Female","Male")))]
mod24_table[is.na(HR4),HR4:="-"]
setorder(mod24_table,"Disorder","Cause","Sex")

mod24_table <- mod24_table[,.(Disorder,Sex,Cause,HR2,HR4)]

rm(mod2_table,mod4_table)

###### Everything in one table ##########

setnames(mod13_table,c("HR1","HR3"),c("HRu","HRa"))
setnames(mod24_table,c("HR2","HR4"),c("HRu","HRa"))
mod13_table[,Sex:="Both"]
mod24_table[,`:=`(Sex=as.character(Sex),p_val="-")]

full_table <- rbind(mod24_table,mod13_table,use.names=TRUE)
full_table <- full_table[,.(Disorder,Cause,Sex,HRu,HRa,p_val)]
full_table[,Sex:=factor(Sex,levels=c("Both","Female","Male"))]
full_table <- full_table[!is.na(Disorder)]        # removing exposures that we don't want to appear in table
setorder(full_table,"Disorder","Cause","Sex")
full_table[duplicated(paste(Disorder,Cause)),Cause:=""]
full_table[duplicated(Disorder),Disorder:=""]
full_table[,`:=`(HRu=gsub("\\.","·",HRu),
                 HRa=gsub("\\.","·",HRa),
                 p_val=gsub("\\.","·",p_val))]
write_xlsx(full_table,file.path(filepath_write,"aHRs_any_admission_table.xlsx"))

