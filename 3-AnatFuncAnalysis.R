################ import
library(doMC)
registerDoMC(detectCores())
library(data.table)
library(pracma)
library(coin)
library(magrittr)

############### base information prepare
covars.all <- fread('subjects_71_info.csv')
covars.all[, gender := as.factor(gender)]
covars <- covars.all[,c('age', 'gender', 'group')]
covars.id <- paste0('sub-', covars.all$participant_id)
############### data prepare
## graph prepare
#df.func.g <- fread('Func_abs_binary_unweighted_G.csv') # graph g
#df.func.v <- fread('Func_abs_binary_unweighted_V.csv') # graph v
df.anat.g <- fread('Thickness_resids_binary_unweighted_G.csv') # graph g
df.anat.v <- fread('Thickness_resids_binary_unweighted_V.csv') # graph v

df.anat.g[, .(threshold, AUC=assort.class, group)] %>%
  .[,group:=as.factor(group)] %>%
  .[, .(Wilp=wilcox.test(AUC~group, paired=F)$p.value,
        Ttp=t.test(AUC~group, paired=F)$p.value,
        Permp=pvalue(oneway_test(AUC~group, distribution=approximate(nresample=10000)))[1])]

[1] "threshold"        "Cp"               "Lp"               "E.global"         "mod"              "density"         
[7] "max.comp"         "num.tri"          "diameter"         "transitivity"     "assort"           "assort.lobe"     
[13] "assort.lobe.hemi" "asymm"            "spatial.dist"     "assort.class"     "num.hubs"         "E.local"         
[19] "vulnerability"    "atlas"            "modality"         "group"   

df.anat.g[, .(threshold, AUC=asymm, group)] %>%
  .[,group:=as.factor(group)] %>%
  .[,.(mean(AUC)), by=group]

####################### ####################### ####################### #######################
####################### statistic - sub level
## 1. the structural covariance network indicated us that cortical
## function displayed a degeneration. So, there has a idea
## Hypothesis 1: the more degeneration, the severe disease / long disease course
df.sub.g <- fread('Func_abs_binary_unweighted_G.csv') # graph g
df.sub.v <- fread('Func_abs_binary_unweighted_V.csv') # graph v

### compute indicators in AUC
df.sub.g[, .(threshold, assort, participant_id, group)] %>%
        .[, .(AUC=trapz(threshold, assort), group), by=participant_id] %>%
            unique(., by='participant_id') %>%
               .[,group:=as.factor(group)] %>%
                .[, .(Wilp=wilcox.test(AUC~group, paired=F)$p.value,
                      Ttp=t.test(AUC~group, paired=F)$p.value,
                      Permp=pvalue(oneway_test(AUC~group, distribution=approximate(nresample=10000)))[1])]

[1] "threshold"        "Cp"               "Lp"               "E.global"         "mod"              "density"         
[7] "max.comp"         "num.tri"          "diameter"         "transitivity"     "assort"           "assort.lobe"     
[13] "assort.lobe.hemi" "asymm"            "spatial.dist"     "assort.class"     "num.hubs"         "E.local"         
[19] "vulnerability"    "participant_id"   "atlas"            "modality"         "weighting"        "group"  

### compute GLM comparison
covars_compare <- df.sub.g[, .(threshold, Indicator=vulnerability, participant_id, group)] %>%
                       .[, .(AUC=trapz(threshold, Indicator)), by=participant_id] %>%
                        covars.all[., on=.(participant_id=participant_id)] %>%
                              .[, group:=as.factor(group)]
with(covars_compare,
  summary(aov(AUC~ group + age + gender + BMI))
)
cor.test(covars_compare$AUC, covars_compare$age, method = 'pearson')

[1] "threshold"        "Cp"               "Lp"               "E.global"         "mod"              "density"         
[7] "max.comp"         "num.tri"          "diameter"         "transitivity"     "assort"           "assort.lobe"     
[13] "assort.lobe.hemi" "asymm"            "spatial.dist"     "assort.class"     "num.hubs"         "E.local"         
[19] "vulnerability"    "participant_id"   "atlas"            "modality"         "weighting"        "group" 
####################### ####################### ####################### #######################
####################### statistic - group level
df.func.g <- fread('FuncGroup_abs_binary_unweighted_G.csv') # graph g
df.func.v <- fread('FuncGroup_abs_binary_unweighted_V.csv') # graph v

df.func.g[, .(threshold, AUC=asymm, group)] %>%
  .[,group:=as.factor(group)] %>%
  .[, .(Wilp=wilcox.test(AUC~group, paired=F)$p.value,
        Ttp=t.test(AUC~group, paired=F)$p.value,
        Permp=pvalue(oneway_test(AUC~group, distribution=approximate(nresample=10000)))[1])]

[1] "threshold"        "Cp"               "Lp"               "E.global"         "mod"              "density"         
[7] "max.comp"         "num.tri"          "diameter"         "transitivity"     "assort"1          "assort.lobe"     
[13] "assort.lobe.hemi" "asymm"1            "spatial.dist"     "assort.class"     "num.hubs"1         "E.local"         
[19] "vulnerability"    "participant_id"   "atlas"            "modality"         "weighting"        "group"    

df.func.g[, .(threshold, AUC=assort.class, group)] %>%
  .[,group:=as.factor(group)] %>%
  .[,.(mean(AUC)), by=group]

####################### statistic - sub level
df.sub.g <- fread('Func_abs_binary_unweighted_G.csv') # graph g
df.sub.v <- fread('Func_abs_binary_unweighted_V.csv') # graph v
# indicator
patient_df <-  df.sub.g[group=='patient', .(threshold, assort.class, participant_id, group)] %>%
                              .[,.(AUC=trapz(threshold, assort.class)), by=participant_id]
cor.test(patient_df$AUC, covars.all[group=='patient',]$disease_month, method = 'pearson')
cor.test(patient_df$AUC, covars.all[group=='patient',]$disease_month, method = 'spearman')
cor.test(patient_df$AUC, covars.all[group=='patient',]$`SF-36`, method = 'pearson')
cor.test(patient_df$AUC, covars.all[group=='patient',]$`SF-36`, method = 'spearman')
cor.test(patient_df$AUC, covars.all[group=='patient',]$PRI, method = 'pearson')
cor.test(patient_df$AUC, covars.all[group=='patient',]$PRI, method = 'spearman')
cor.test(patient_df$AUC, covars.all[group=='patient',]$`FS-14`, method = 'pearson')
cor.test(patient_df$AUC, covars.all[group=='patient',]$`FS-14`, method = 'spearman')


                        



[1] "threshold"        "Cp"               "Lp"               "E.global"         "mod"              "density"         
[7] "max.comp"         "num.tri"          "diameter"         "transitivity"     "assort"           "assort.lobe"     
[13] "assort.lobe.hemi" "asymm"            "spatial.dist"     "assort.class"1     "num.hubs"         "E.local"         
[19] "vulnerability"    "participant_id"   "atlas"            "modality"         "weighting"        "group"    
