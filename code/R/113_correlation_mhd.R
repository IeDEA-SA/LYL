# correlation + percentage matrix for the different types of mental illness

library(data.table)
library(survival)
library(tictoc)
library(tableone)
library(writexl)
library(reshape2)
library(ggplot2)

tic()

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfAc_excess_mortality/processed"
filepath_write <- "C:/ISPM/HomeDir/HIV-mental disorders/Anja/output/tables/descriptive"
filepath_plot <- "C:/ISPM/HomeDir/HIV-mental disorders/Anja/output/plots/comorbidity"
source(file="C:/ISPM/HomeDir/HIV-mental disorders/R/Code/utils/timeSplit_DT.R")

which_mhd<- c("first_any_admission_organic","first_any_admission_substance_use_disorder","first_any_admission_psychotic","first_any_admission_mood_disorder",
              "first_any_admission_anxiety","first_any_admission_asc_developmental_disorders","first_any_admission_asc_personality_disorder")
which_mhd_short <- c("Org","SU","Psy","Mood","Anx","Dev","Pers")

which_mhd <- rev(which_mhd)
which_mhd_short <- rev(which_mhd_short)

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

# cleaning up causes of death AFTER censoring
data_surv[event==0,cod2:="Alive"]

# formatting
data_surv[,`:=`(age_base=start,
                age_base_cat=cut(start,breaks=c(15,20,30,40,50,60,70,85),right=FALSE),
                calyear_base_cat=cut(year(birth_d+start*365.25),breaks=c(2010,2014,2017,2020,Inf),right=FALSE),
                sex=factor(sex,levels=c("Male","Female")),
                cod2=factor(cod2,levels=c("Natural","Unnatural","Unknown","Alive")))]

# looping over types of mental disorders
MHD_count_df <- data.frame(NULL)

N_total <- data_surv[,.N]
N_men <- data_surv[sex=="Male",.N]
N_women <- data_surv[sex=="Female",.N]

DT <- copy(data_surv)

for(v in which_mhd)
{
  mhd <- gsub("first_any_admission_","",v)
  
  # age at mhd treatment
  DT[,`:=`(age_treat=as.numeric(get(v)-birth_d)/365.25,py=stop-start)]
  DT[age_treat>=stop,age_treat:=NA]                              # excluding disorders occuring post-censoring  
  DT[,eval(mhd):=as.numeric(!is.na(age_treat))]
  DT[,age_treat:=NULL]
}

mhd_vect <- gsub("first_any_admission_","",which_mhd)
DT <- DT[,..mhd_vect]

corr_matrix_out <- cbind(data.frame(disorder=colnames(DT)),data.frame(cor(DT)))

write_xlsx(corr_matrix_out,path=file.path(filepath_write,"MHD_correlation_matrix.xlsx"))

pct_matrix <- matrix(NA,nrow=ncol(DT),ncol=ncol(DT))
for(i in 1:length(mhd_vect))
 for (j in 1:length(mhd_vect)) 
  pct_matrix[i,j] <- DT[get(mhd_vect[i])==1 & get(mhd_vect[j])==1,.N]/DT[get(mhd_vect[i])==1,.N]

pct_matrix <- data.frame(pct_matrix)
colnames(pct_matrix) <- mhd_vect
pct_matrix_out <- cbind(data.frame(DISORDER=mhd_vect),pct_matrix)
row.names(pct_matrix_out) <- NULL
write_xlsx(pct_matrix_out,path=file.path(filepath_write,"MHD_percentage_matrix.xlsx"))

DT[,N_MHD:=organic+substance_use_disorder+psychotic+mood_disorder+anxiety+asc_developmental_disorders+asc_personality_disorder]
DTmhd <- DT[N_MHD>0]
DTmhd[,prop.table(table(N_MHD>1))]

# plotting percentages (heatmap)
pct_matrix <- round(100*as.matrix(pct_matrix),digits=1)
rownames(pct_matrix) <- which_mhd_short
colnames(pct_matrix) <- which_mhd_short
melted_pct_matrix <- melt(pct_matrix)

melted_pct_matrix$value_txt <- paste0(melted_pct_matrix$value,"%")

melted_pct_matrix$Var1 <- factor(melted_pct_matrix$Var1)

pp <- ggplot(data=melted_pct_matrix,aes(x=Var2,y=Var1,fill=value)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low = "white",high = "red",limit = c(0,100)) +
  geom_text(aes(Var2, Var1, label = value_txt), color = "black", size = 3.5) +
  theme(axis.text=element_text(size=16),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

ggsave(pp,filename=file.path(filepath_plot,"mhd_percent_matrix.png"),height=8,width=7,dpi=600)


