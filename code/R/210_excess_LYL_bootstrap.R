# Bootstrap samples for excess LYLs, for a given exposure and sex
# To be run as job arrays on UBELIX (HPC server at University of Bern)
# ~4 hours for 50 iterations
# for sensitivity analyses, change 'long_table.RData' to appropriate dataset, and change the name of the savefiles, 
# for analysis including only single-diagnosis disorders, reincorporate the appropriate lines

library(data.table)
library(survival)
library(lillies)
library(readxl)
library(tictoc)
library(scriptName)

tic()

filepath_read <- "~/MHD/R/input/AfAc"
filepath_write <- "~/MHD/R/output/lyl"
filepath_lookup <- "~/MHD/R/lookup"            # see'lookup' folder of github repository
source("~/MHD/R/Code/Utils/timeSplit_DT.R")

min_age <- 15
max_age <- 85

# looking up parameters from Excel file, these will be determined by the number in the filename of the script currently running
mhd_lookup <- readxl::read_xlsx(path=file.path(filepath_lookup,"lookup.xlsx"))
x <- unlist(strsplit(current_filename(),"/"))
id_num <- as.numeric(gsub("[^\\d]+", "",x[length(x)], perl=TRUE))
which_exposure <- mhd_lookup$exposure[mhd_lookup$id==id_num]
which_sex <- mhd_lookup$sex[mhd_lookup$id==id_num]
chunk <- mhd_lookup$chunk[mhd_lookup$id==id_num]
nb_boot <- mhd_lookup$nb_boot[mhd_lookup$id==id_num]

set.seed(chunk*16)       # creating RNG

save_filename <- paste0("elyl_boot_",which_exposure,"_",which_sex,"_c",chunk,".RData")
print(save_filename)

load(file=paste(filepath_read,"/long_table.RData",sep=''))
data_surv <- copy(BAS_MHD_new)

# removing patients with unknown sex or missing date of birth
data_surv <- data_surv[(sex=="Male" | sex=="Female") & !is.na(birth_d)]

# counting process format, one line per patient, on 'age' scale
data_surv[,`:=`(start=(start-birth_d)/365.25,
                stop=(end-birth_d)/365.25)]

# left-truncataion/right-censoring (age)
data_surv <- data.table(survSplit(Surv(start,stop,event)~.,data=data_surv,cut=c(min_age,max_age),episode="i"))
print(data_surv[is.na(start) | is.na(stop) | is.na(event) | is.na(i),.(patient,start,stop,event,i)])
data_surv <- data_surv[i==2]
data_surv[,i:=NULL]

# cleaning up AFTER censoring
data_surv[event==0,cod2:="Alive"]

data_surv[,cod2:=factor(cod2,levels=c("Alive","Natural","Unnatural","Unknown"))]

LYL_diff_boot <- data.frame(NULL)
no_deaths_vect <- NULL

# Bootstrap sampling
for(i in 1:nb_boot)
{
 print(i)

 DT_boot <- data_surv[sample(1:.N,replace=TRUE)]    # generating bootstrap sample
 if(which_sex=="men")
  DT_boot <- DT_boot[sex=="Male"]
 if(which_sex=="women")
  DT_boot <- DT_boot[sex=="Female"]
  
 # excluding patients with two or more diagnoses of the disorder (sensitivity analysis only!)
 #N_diag <- gsub("first_","N_",which_exposure)
 #DT_boot <- DT_boot[get(N_diag)<2]
 
 # age at first mhd treatment of any kind and from any source, or at first psychiatric medication
 if(substr(which_exposure,1,19)=="first_any_admission" || substr(which_exposure,1,9)=="first_hos")
 {
   DT_boot[,age_any_treat:=as.numeric(first_any_admission-birth_d)/365.25]
 } else
 {
   DT_boot[,age_any_treat:=as.numeric(first_any_med-birth_d)/365.25]
 }
 
 setnames(DT_boot,which_exposure,"treat_d")
 
 # age at mhd treatment
 DT_boot[,`:=`(age_treat=as.numeric(treat_d-birth_d)/365.25)]
 DT_boot[age_treat>=stop,age_treat:=NA]                             # excluding disorders occuring post-censoring  
 DT_boot[!is.na(age_treat) & age_treat<start,age_treat:=start]      # moving forward prevalent exposures to cohort entry
 age_treat <- DT_boot[!is.na(age_treat),age_treat]
 
 # splitting time-at-risk into exposed/unexposed
 DT_exp <- DT_boot[!is.na(age_treat)]
 DT_unexp <- timeSplit_DT(X=DT_boot,vars_date="age_any_treat",vars_tu="cens",event="event",
                          start_date="start",stop_date="stop",id_var="patient",print_out=FALSE)
 DT_unexp <- DT_unexp[cens==0]
 DT_unexp[,cens:=NULL]
 DT_unexp[event==0,cod2:="Alive"]
 
 if(DT_exp[,any(cod2!="Alive")])   # checking that there is at least one death among the exposed, otherwise ELYL must be determined manually, without lillies
 {
 # computing excess LYL, exposed vs. unexposed
   LYL_exp <- suppressMessages(lyl_range(data=DT_exp,t0=age_treat,t=stop,status=cod2,age_begin=min_age,age_end=max_age-1,tau=max_age))
   LYL_unexp <- suppressMessages(lyl_range(data=DT_unexp,t0=start,t=stop,status=cod2,age_begin=min_age,age_end=max_age-1,tau=max_age))
   LYL_diff <- lyl_diff(LYL_exp,LYL_unexp,weights=age_treat)
   LYL_diff_boot <- rbind(LYL_diff_boot,data.frame(LYL_diff))
 } else        # in case of no deaths in exposed, LYL is zero in that population
 {
   print("no deaths in exposed!")
   no_deaths_vect <- c(no_deaths_vect,i)
   W <- DT_exp[,floor(age_treat)]
   suppressMessages(LYL_unexp <- lyl_range(data=DT_unexp,t0=start,t=stop,status=cod2,age_begin=min_age,age_end=max_age-1,tau=max_age))
   LYL <- LYL_unexp$LYL
   ind <- match(W,LYL$age)
   lyl_nd <- -mean(LYL$Natural[ind])
   lyl_ud <- -mean(LYL$Unnatural[ind])
   lyl_uk <- -mean(LYL$Unknown[ind])
   LYL_diff_boot <- rbind(LYL_diff_boot,data.frame(life_exp=-(lyl_nd+lyl_ud+lyl_uk),TotalLYL=lyl_nd+lyl_ud+lyl_uk,
                                                   Natural=lyl_nd,Unnatural=lyl_ud,Unknown=lyl_uk))
   rm(W,LYL,ind,lyl_nd,lyl_ud,lyl_uk)
 }
 rm(DT_exp,DT_unexp,DT_boot,age_treat)
 gc()
}

print(no_deaths_vect)
write.table(LYL_diff_boot,file=file.path(filepath_write,save_filename))

toc()
