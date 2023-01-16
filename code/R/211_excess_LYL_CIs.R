# produce 95% CIs for excess LYL associated with each exposure, based on bootstrap sampling
# creates a separate file of CIs for each source (any admission, hospital, medication)
# the bootstrapping is done upstream on UBELIX via job arrays

library(data.table)

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfAc_excess_mortality/lillies/bootstrap"
filepath_write <- "C:/ISPM/Data/HIV-mental disorders/AfAc_excess_mortality/lillies/bootstrap/joined"

which_mhd_any <- c("first_any_admission","first_any_admission_organic","first_any_admission_substance_use_disorder","first_any_admission_alcohol_use_disorder",
                   "first_any_admission_drug_use_disorder","first_any_admission_psychotic","first_any_admission_mood_disorder","first_any_admission_bipolar",
                   "first_any_admission_depression","first_any_admission_anxiety","first_any_admission_generalised_anxiety_disorder","first_any_admission_PTSD",
                   "first_any_admission_asc_developmental_disorders","first_any_admission_asc_behavioural_disorders","first_any_admission_asc_physical_factors",
                   "first_any_admission_eating_disorder","first_any_admission_asc_personality_disorder")

which_mhd_hos<- c("first_hos_any","first_hos_organic","first_hos_substance_use_disorder","first_hos_alcohol_use_disorder","first_hos_drug_use_disorder",
                  "first_hos_psychotic","first_hos_mood_disorder","first_hos_bipolar","first_hos_depression","first_hos_anxiety",
                  "first_hos_generalised_anxiety_disorder","first_hos_PTSD","first_hos_asc_developmental_disorders","first_hos_asc_behavioural_disorders",
                  "first_hos_asc_physical_factors","first_hos_eating_disorder","first_hos_asc_personality_disorder")

which_mhd_med <- c("first_any_med","first_substance_use_med","first_antipsychotic","first_antidepressant","first_anxiolytic")

which_sex <- c("men","women")

# merging bootstrap samples per data source

eLYLboot_any_admission <- data.table(NULL)
for(v in which_mhd_any)
  for(w in which_sex)
     for(i in 1:10)
     {
        boot_sample <- data.table(read.table(file.path(filepath_read,"any",paste0("elyl_boot_",v,"_",w,"_c",i,".RData"))))
        boot_sample <- boot_sample[,.(exposure=v,sex=w,AllCause=TotalLYL,Natural,Unnatural,Unknown)]
        eLYLboot_any_admission <- rbind(eLYLboot_any_admission,boot_sample)
        rm(boot_sample)
     }
eLYLboot_any_admission[sex=="men",sex:="Male"]
eLYLboot_any_admission[sex=="women",sex:="Female"]
save(eLYLboot_any_admission,file=file.path(filepath_write,"eLYLboot_any_admission.RData"))

eLYLboot_hos <- data.table(NULL)
for(v in which_mhd_hos)
   for(w in which_sex)
      for(i in 1:10)
      {
         boot_sample <- data.table(read.table(file.path(filepath_read,"hos",paste0("elyl_boot_",v,"_",w,"_c",i,".RData"))))
         boot_sample <- boot_sample[,.(exposure=v,sex=w,AllCause=TotalLYL,Natural,Unnatural,Unknown)]
         eLYLboot_hos <- rbind(eLYLboot_hos,boot_sample)
         rm(boot_sample)
      }
eLYLboot_hos[sex=="men",sex:="Male"]
eLYLboot_hos[sex=="women",sex:="Female"]
save(eLYLboot_hos,file=file.path(filepath_write,"eLYLboot_hos.RData"))

eLYLboot_med <- data.table(NULL)
for(v in which_mhd_med)
  for(w in which_sex)
    for(i in 1:10)
    {
      boot_sample <- data.table(read.table(file.path(filepath_read,"med",paste0("elyl_boot_",v,"_",w,"_c",i,".RData"))))
      boot_sample <- boot_sample[,.(exposure=v,sex=w,AllCause=TotalLYL,Natural,Unnatural,Unknown)]
      eLYLboot_med <- rbind(eLYLboot_med,boot_sample)
      rm(boot_sample)
    }
eLYLboot_med[sex=="men",sex:="Male"]
eLYLboot_med[sex=="women",sex:="Female"]
save(eLYLboot_med,file=file.path(filepath_write,"eLYLboot_med.RData"))

# 95% CIs

CIs_any_admission <- eLYLboot_any_admission[,.(lcl_allcause=quantile(AllCause,0.025),ucl_allcause=quantile(AllCause,0.975),
                                               lcl_natural=quantile(Natural,0.025),ucl_natural=quantile(Natural,0.975),
                                               lcl_unnatural=quantile(Unnatural,0.025),ucl_unnatural=quantile(Unnatural,0.975),
                                               lcl_uk=quantile(Unknown,0.025),ucl_uk=quantile(Unknown,0.975)),
                                               by=.(exposure,sex)]

save(CIs_any_admission,file=file.path(filepath_write,"CIs_any_admission.RData"))

CIs_hos <- eLYLboot_hos[,.(lcl_allcause=quantile(AllCause,0.025),ucl_allcause=quantile(AllCause,0.975),
                           lcl_natural=quantile(Natural,0.025),ucl_natural=quantile(Natural,0.975),
                           lcl_unnatural=quantile(Unnatural,0.025),ucl_unnatural=quantile(Unnatural,0.975),
                           lcl_uk=quantile(Unknown,0.025),ucl_uk=quantile(Unknown,0.975)),
                           by=.(exposure,sex)]

save(CIs_hos,file=file.path(filepath_write,"CIs_hos.RData"))
                    
CIs_med <- eLYLboot_med[,.(lcl_allcause=quantile(AllCause,0.025),ucl_allcause=quantile(AllCause,0.975),
                           lcl_natural=quantile(Natural,0.025),ucl_natural=quantile(Natural,0.975),
                           lcl_unnatural=quantile(Unnatural,0.025),ucl_unnatural=quantile(Unnatural,0.975),
                           lcl_uk=quantile(Unknown,0.025),ucl_uk=quantile(Unknown,0.975)),
                           by=.(exposure,sex)]

save(CIs_med,file=file.path(filepath_write,"CIs_med.RData"))
