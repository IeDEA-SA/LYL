# Single plot with two panels: A) eLYL any admission, B) eLYL hospitalizations

library(data.table)
library(ggplot2)

filepath_lookup <- "C:/ISPM/HomeDir/HIV-mental disorders/Anja/lookup"
filepath_est <- "C:/ISPM/Data/HIV-mental disorders/AfAc_excess_mortality/lillies/point_estimates"
filepath_boot <- "C:/ISPM/Data/HIV-mental disorders/AfAc_excess_mortality/lillies/bootstrap/joined"
filepath_plot <- "C:/ISPM/HomeDir/HIV-mental disorders/Anja/output/plots/LYL"

expo <- read.csv(file.path(filepath_lookup,"exposures.csv"),header=TRUE, sep=",")

# ordered exposures, nested disorders at the end
target <- c("Any Mental Health Diagnosis", "Organic Disorders","Substance Use Disorders","Alcohol Use Disorders","Drug Use Disorders",
            "Psychotic Disorders","Mood Disorders","Bipolar Disorders","Depressive Disorders","Anxiety Disorders","Generalised Anxiety Disorders",
            "PTSD","Developmental Disorders","Eating Disorders","Personality Disorders")

##### Any admission #####

load(file.path(filepath_est,"elyl_any_admission.RData"))
DF <- data.table(LYL_diff_df)
rm(LYL_diff_df)
load(file.path(filepath_boot,"CIs_any_admission.RData"))
DF_CI <- data.table(CIs_any_admission)
rm(CIs_any_admission)

#Remove columns containing "BOTH"

DF_m_f <- DF[sex!= "Both"]
DF_CI_m_f<-  DF_CI[sex!= "Both"]

#merge lables dataframe with DF_m_f

DF_new <- merge(DF_m_f,expo, by="type")
DF_new[,type:=as.character(type)]

#Rename exposure to type in DF_CI_m_f to merge later
setnames(DF_CI_m_f,"exposure","type")

# recoding sex variable
DF_CI_m_f[sex=="Male",sex:="Men"]
DF_CI_m_f[sex=="Female",sex:="Women"]

# appending CIs
DF_new_new <- merge(DF_new,DF_CI_m_f, by=c("type","sex"))

#Add column for admission type
DF_new_new[,admission_type:="Any health care setting"]

###### Hospitalizations ##########

load(file.path(filepath_est,"elyl_hos.RData"))
DF <- data.table(LYL_diff_df)
rm(LYL_diff_df)
load(file.path(filepath_boot,"CIs_hos.RData"))
DF_CI <- data.table(CIs_hos)
rm(CIs_hos)

#Remove columns containing "BOTH"

DF_m_f <- DF[sex!= "Both"]
DF_CI_m_f<-  DF_CI[sex!= "Both"]

#merge lables dataframe with DF_m_f

DF_new <- merge(DF_m_f,expo, by="type")
DF_new[,type:=as.character(type)]

#Rename exposure to type in DF_CI_m_f to merge later
setnames(DF_CI_m_f,"exposure","type")

# recoding sex variable
DF_CI_m_f[sex=="Male",sex:="Men"]
DF_CI_m_f[sex=="Female",sex:="Women"]

# appending CIs
DF_new_new_hos <- merge(DF_new,DF_CI_m_f, by=c("type","sex"))

#Add column for admission type
DF_new_new_hos[,admission_type:="Hospital setting"]

#####  Combine the two data sources into one plot ######

DF_combined <- rbind(DF_new_new,DF_new_new_hos)
# negative CI values set to 0
DF_combined[lcl_allcause<0,lcl_allcause:=0]
DF_combined[ucl_allcause<0,ucl_allcause:=0]
# negative point estimate -> remove entirely
DF_combined[TotalLYL<0,`:=`(TotalLYL=NA,lcl_allcause=NA,ucl_allcause=NA)]
# hospitalization for PTSD has negative eLYL in women -> set to missing
DF_combined[,d_lables:=factor(d_lables,levels=target)]
DF_combined <- DF_combined[!is.na(d_lables)]            # removing exposures that we don't want shown in the plot

p1 <- ggplot(data=DF_combined,aes(x=TotalLYL,y= d_lables, fill=sex))+
  facet_grid(cols=vars(admission_type))+ 
  geom_col(data=subset(DF_combined,sex=="Men"), aes(x=TotalLYL*(-1),fill="Men"),width=0.75) +
  geom_errorbarh(data=subset(DF_combined,sex=="Men"),aes(xmin = lcl_allcause*(-1), xmax =ucl_allcause*(-1) ),size=.5, colour = 'black',height = .2) +
  geom_col(data=subset(DF_combined,sex=="Women"),aes(fill="Women"),width=0.75) +
  geom_errorbarh(data=subset(DF_combined,sex=="Women"),aes(xmin = lcl_allcause, xmax =ucl_allcause),size=.5, colour = 'grey57',height = .2) +
  scale_x_continuous(breaks=seq(-15,15,5),labels=abs(seq(-15,15,5))) +  
  scale_y_discrete(limits =rev(target)) +
  scale_fill_manual(labels=c("Men","Women"),values=c(rgb(159,174,213,maxColorValue=255),rgb(226,165,159,maxColorValue=255)))+
  theme_bw() +
  theme(
    plot.title= element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.background = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y=element_blank())+
  ylab(NULL)+
  xlab(NULL)+
  geom_vline(xintercept = 0)
  #geom_hline(yintercept = c(6.5,14.5), linetype="dotted")

ggsave(p1,filename=file.path(filepath_plot,"ELYL_any_admission_and_hos.png"),height=8,width=7,dpi=600)