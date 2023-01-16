# produce 95% CIs for excess LYL associated with each exposure, based on bootstrap sampling
# exposures = one diagnosis only
# the bootstrapping is done upstream on UBELIX via job arrays

library(data.table)

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfAc_excess_mortality/lillies/bootstrap"
filepath_write <- "C:/ISPM/Data/HIV-mental disorders/AfAc_excess_mortality/lillies/bootstrap/joined"

which_mhd_any <- c("first_any_admission","first_any_admission_organic","first_any_admission_substance_use_disorder","first_any_admission_alcohol_use_disorder",
                   "first_any_admission_drug_use_disorder","first_any_admission_psychotic","first_any_admission_mood_disorder","first_any_admission_bipolar",
                   "first_any_admission_depression","first_any_admission_anxiety","first_any_admission_generalised_anxiety_disorder","first_any_admission_PTSD",
                   "first_any_admission_asc_developmental_disorders","first_any_admission_asc_behavioural_disorders","first_any_admission_asc_physical_factors",
                   "first_any_admission_eating_disorder","first_any_admission_asc_personality_disorder")

which_sex <- c("men","women")

# merging bootstrap samples per data source

eLYLboot_any_admission <- data.table(NULL)
for(v in which_mhd_any)
  for(w in which_sex)
     for(i in 1:10)
     {
       filename <- file.path(filepath_read,paste0("elyl_boot_1diag_",v,"_",w,"_c",i,".RData"))
       if(file.exists(filename))
       {
         boot_sample <- data.table(read.table(filename))
         boot_sample <- boot_sample[,.(exposure=v,sex=w,AllCause=TotalLYL,Natural,Unnatural,Unknown)]
         eLYLboot_any_admission <- rbind(eLYLboot_any_admission,boot_sample)
         rm(boot_sample)
       } else
       {
         print(paste0(filename," not found"))
       }
     }

eLYLboot_any_admission[sex=="men",sex:="Male"]
eLYLboot_any_admission[sex=="women",sex:="Female"]
save(eLYLboot_any_admission,file=file.path(filepath_write,"eLYLboot_any_admission_1diag.RData"))


# 95% CIs

CIs_any_admission_2diag <- eLYLboot_any_admission[,.(lcl_allcause=quantile(AllCause,0.025),ucl_allcause=quantile(AllCause,0.975),
                                               lcl_natural=quantile(Natural,0.025),ucl_natural=quantile(Natural,0.975),
                                               lcl_unnatural=quantile(Unnatural,0.025),ucl_unnatural=quantile(Unnatural,0.975),
                                               lcl_uk=quantile(Unknown,0.025),ucl_uk=quantile(Unknown,0.975)),
                                               by=.(exposure,sex)]

save(CIs_any_admission_2diag,file=file.path(filepath_write,"CIs_any_admission_1diag.RData"))
