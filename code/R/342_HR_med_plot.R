# Forest Plots from Cox Regression, ATC medication
# Adjusted hazard ratios, unadjusted and adjusted for psychiatric morbidity, by sex

library(data.table)
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyverse)
library(cowplot)

filepath_lookup <- "C:/ISPM/HomeDir/HIV-mental disorders/AfAc_excess_mortality/LYL/lookup"   # found in 'lookup' folder of GitHub rep
filepath_read <- "C:/ISPM/HomeDir/HIV-mental disorders/AfAc_excess_mortality/output/tables/model output"
filepath_plot <- "C:/ISPM/HomeDir/HIV-mental disorders/AfAc_excess_mortality/output/plots"

# ordered exposures
target <- c("Any Medication", "Substance Use Medication","Antipsychotic Medication","Antidepressant Medication","Anxiolytic Medication")

#read in data
gender_lookup <- read.csv(file.path(filepath_lookup,"sex.csv"),header=TRUE, sep=",")
expo <- read.csv(file.path(filepath_lookup,"exposures.csv"),header=TRUE, sep=",")
mod1_any <- read_xlsx(file.path(filepath_read,"aHRs_model1_med.xlsx"))
mod2_any <- read_xlsx(file.path(filepath_read,"aHRs_model2_med.xlsx"))
mod3_any <- read_xlsx(file.path(filepath_read,"aHRs_model3_med.xlsx"))
mod4_any <- read_xlsx(file.path(filepath_read,"aHRs_model4_med.xlsx"))

names(expo)[1] <- "disorder"
gender_lookup <- rbind(gender_lookup,c("Both","Both sexes"))

#Row bind 1&2

mod1_any$p_int_sex_exp <- NULL
mod2_any$ucl <- as.numeric(mod2_any$ucl)
modA <- bind_rows(mod1_any,mod2_any)

##Row bind 3&4
mod3_any$p_int_sex_exp <- NULL
mod4_any$ucl <- as.numeric(mod4_any$ucl)   # for some reason the upper limit for this model is a character class
modB <- bind_rows(mod3_any,mod4_any)

# removing cod/exposure/sex combinations with no event in exposure group
modA <- data.table(modA)
modA[ucl>200,`:=`(estimate=NA,lcl=NA,ucl=NA)]
print("Plot A: exposure/sex/outcome combinations with no events:")
print(modA[is.na(estimate),.(disorder,sex,cause)])
modA <- as_tibble(modA)

modB <- data.table(modB)
modB[ucl>200,`:=`(estimate=NA,lcl=NA,ucl=NA)]
print("Plot B: exposure/sex/outcome combinations with no events:")
print(modB[is.na(estimate),.(disorder,sex,cause)])
modB <- as_tibble(modB)

#Add column with new labels
modA<-merge(modA,expo, by="disorder")
modA<-merge(modA,gender_lookup, by="sex")

#Add column with new labels
modB<-merge(modB,expo, by="disorder")
modB<-merge(modB,gender_lookup, by="sex")

#Change cause all to All
modA$cause[modA$cause == "all"] <- "All"
modB$cause[modB$cause == "all"] <- "All" 

#prepare data set for displaying causes of death

modA$gender <- factor(modA$gender,levels=c("Men","Women","Both sexes"))

#prepare data set for displaying causes of death

modB$gender <- factor(modB$gender,levels=c("Men","Women","Both sexes"))

modA$d_lables <- factor(modA$d_lables,levels=rev(target))

modB$d_lables <- factor(modB$d_lables,levels=rev(target))

#add column to address adjustment
attach(modA)
modA$adjusted<-"no"
detach(modA)
attach(modB)
modB$adjusted<-"yes"
detach(modB)

#Create long table

mod_any <- data.table(rbind(modA,modB))

mod_any <- mod_any[!is.na(d_lables) & cause=="All"]
mod_any[,`:=`(adjusted=factor(adjusted,levels=c("yes","no")),
              gender=factor(gender,levels=c("Men","Women","Both sexes")))]

# appending 'empty' estimates for those exposures not in fully adjusted model
combos <- data.table(crossing(target,c("Men","Women","Both sexes"),c("yes","no")))
names(combos) <- c("d_lables","gender","adjusted")
combos[,`:=`(d_lables=factor(d_lables,levels=rev(target)),
             adjusted=factor(adjusted,levels=c("yes","no")),
             gender=factor(gender,levels=c("Men","Women","Both sexes")))]
mod_any <- merge(mod_any,combos,all.y=TRUE)

#Create plot

p <- ggplot(data=mod_any,aes(x=estimate,y=d_lables))+
  geom_point(aes(color=adjusted),position=ggstance::position_dodgev(height=0.6))+
  facet_grid(cols=vars(gender))+  
  geom_errorbarh(aes(xmin=lcl,xmax=ucl,color=adjusted),size=1,height=.2,position=ggstance::position_dodgev(height=0.6)) +
  theme_bw() +
  theme(
    plot.title= element_text(hjust=0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.background = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(),
    axis.ticks.y=element_blank(),
    strip.text.y = element_blank())+
  geom_vline(xintercept = 1)+
  scale_color_manual(values=c(rgb(52,116,180,maxColorValue=255),rgb(205,76,95,maxColorValue=255)),
                     limits=c("no","yes"),labels=c("Model 1","Model 2")) +
  scale_x_continuous(trans="log2",limits=c(0.5,4),breaks = c(0.5,1,2,4),labels=c(0.5,1,2,4))+
  xlab("")+
  ylab("")

ggsave(p,filename=file.path(filepath_plot,"aHRs_med.png"),height=8,width=7,dpi=600)