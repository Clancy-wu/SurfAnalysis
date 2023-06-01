# import
library(brainGraph)
library(doMC)
registerDoMC(detectCores())
library(data.table)
library(pracma)
library(coin)
library(magrittr)
library(ppcor)
library(ggplot2)
# import graphs
options(bg.subject_id='participant_id', bg.group='group')
grps = c('health', 'patient')
densities <- seq(0.01, 0.34, 0.01)
setkey(destrieux, index)
covars.all <- fread('subjects_71_info.csv')
covars.all[, gender := as.factor(gender)]
inds = lapply(grps, function(x) covars.all[group == x, which = TRUE])
atlas = 'destrieux'
anat.gw.group <- readRDS('brainGraphRDS//anat.gw.group.rds')
anat.gb.group <- readRDS('brainGraphRDS//anat.gb.group.rds')
func.gw <- readRDS('brainGraphRDS//func.gw.rds')
func.gb <- readRDS('brainGraphRDS//func.gb.rds')
func.gw.group <- readRDS('brainGraphRDS//func.gw.group.rds')
func.gb.group <- readRDS('brainGraphRDS//func.gb.group.rds')
fiber.gw <- readRDS('brainGraphRDS//fiber.gw.rds')
fiber.gb <- readRDS('brainGraphRDS//fiber.gb.rds')
fiber.gw.group <- readRDS('brainGraphRDS//fiber.gw.group.rds')
fiber.gb.group <- readRDS('brainGraphRDS//fiber.gb.group.rds')

############## multiple modal network
df.gw.g <- rbindlist(lapply(func.gw.group, graph_attr_dt)) 
df.gw.v <- rbindlist(lapply(func.gw.group, vertex_attr_dt)) 

df.gw.g[,.(num.hubs.wt, asymm, group, threshold)] %>% .[,group:=as.factor(group)] %>%
  .[, .(Wp=wilcox.test(num.hubs.wt~group)$p.value, Tp=t.test(asymm~group)$p.value)]
[1] "threshold"        "Cp"               "Lp"               "E.global"         "mod"              "density"         
[7] "max.comp"         "num.tri"          "diameter"         "transitivity"     "assort"           "strength"        
[13] "E.local.wt"       "E.global.wt"      "diameter.wt"      "Lp.wt"            "num.hubs.wt"      "mod.wt"          
[19] "assort.lobe"      "assort.lobe.hemi" "asymm"            "spatial.dist"     "assort.class"     "num.hubs"        
[25] "E.local"          "vulnerability"    "atlas"            "modality"         "group"       

plot_global(func.gw.group, xvar='threshold')
############## 
[1] "Frontal"   "Occipital" "Parietal"  "Limbic"    "Insula"    "Temporal" 

## sub level - L1 lateralzation
df.gw.sub.v <- rbindlist(lapply(func.gw, vertex_attr_dt))
df.gw.sub.v[, 
            .(mean(E.nodal.wt)), by=.(hemi,participant_id, group, threshold)] %>%
            .[, .(Wp=wilcox.test(V1~group)$p.value, Tp=t.test(V1~group)$p.value), by=.(threshold, hemi)]

df <- df.gw.sub.v[, .(mean(E.nodal.wt)), by=.(hemi,participant_id, group, threshold)] 
df.l <- df[ hemi=='L',]
df.r <- df[hemi=='R',]
df.lateral <- df.l[,.(participant_id,group, threshold)]
df.lateral$lateral <- (df.l$V1 - df.r$V1) / (df.l$V1 + df.r$V1)
df.lateral[, .(Wp=wilcox.test(lateral~group)$p.value, Tp=t.test(lateral~group)$p.value),by=threshold]

## sub level - numbers of hubs
df.gw.sub.v <- rbindlist(lapply(func.gw, vertex_attr_dt))
df <- df.gw.sub.v[hubs.wt>2 & hemi=='R', 
                  .(region, lobe, hemi, hubs.wt, hubnum=1, E.nodal.wt, participant_id, group, threshold)]
df[,.(Sum = sum(hubnum)), by=.(participant_id, group, threshold)] %>%
  .[, .(Wp=wilcox.test(Sum~group)$p.value, Tp=t.test(Sum~group)$p.value), by=threshold]

df[,.(Sum = sum(hubnum)), by=.(participant_id, group, threshold)] %>%
  .[, .(AUC=trapz(threshold, Sum)), by= .(participant_id,group)] %>%
  .[, .(Wp=wilcox.test(AUC~group)$p.value, Tp=t.test(AUC~group)$p.value)]
  
## sub level - hubs nodal efficiency  =====  total efficiency
df.gw.sub.v <- rbindlist(lapply(func.gw, vertex_attr_dt))
df <- df.gw.sub.v[hubs.wt>2 & hemi=='R' ,
                  .(region, lobe, hemi, hubs.wt, hubnum=1, E.nodal.wt, participant_id, group, threshold)]
df[, .(HubSum=sum(hubnum), NodSum=sum(E.nodal.wt)),by=.(participant_id, group, threshold)] %>%
    .[, .(HubSum, NodSum, Mean=(NodSum / HubSum)), by=.(participant_id, group, threshold)] %>%
  .[, .(Wp=wilcox.test(Mean~group)$p.value, Tp=t.test(Mean~group)$p.value), by=threshold]
      
df[, .(HubSum=sum(hubnum), NodSum=sum(E.nodal.wt)),by=.(participant_id, group, threshold)] %>%
  .[, .(HubSum, NodSum, Mean=(NodSum / HubSum)), by=.(participant_id, group, threshold)] %>%
  .[,.(AUC=trapz(threshold, NodSum)), by= .(participant_id,group)] %>%
  .[, .(Wp=wilcox.test(AUC~group)$p.value, Tp=t.test(AUC~group)$p.value)]

## sub level - different lobe
[1] "Frontal"   "Occipital" "Parietal"  "Limbic"    "Insula"    "Temporal" 

df.gw.sub.v <- rbindlist(lapply(fiber.gw, vertex_attr_dt))
# total chagned, E.nodal.wt HC>CFS
df.gw.sub.v <- rbindlist(lapply(func.gw, vertex_attr_dt))
df.gw.sub.v[, 
            .(V1=E.nodal.wt), by=.(region,hemi,participant_id, group, threshold)] %>%
  .[, .(Wp=wilcox.test(V1~group)$p.value, Tp=t.test(V1~group)$p.value), by=.(threshold)]

df.gw.sub.v[, 
            .(V1=E.nodal.wt), by=.(region,hemi,participant_id, group, threshold)] %>%
  .[, .(MeanP=mean(V1), MeanH=mean(V1)), by=.(threshold, group)]

## lobe lateral - no result
df.gw.sub.v <- rbindlist(lapply(func.gw, vertex_attr_dt))
df.gw.sub.v[,  .(V1=hubs.wt), by=.(region,hemi,participant_id, group, threshold)] %>%
  .[, .(Wp=wilcox.test(V1~group)$p.value, Tp=t.test(V1~group)$p.value), by=.(threshold)]
## lateral - E.nodal.wt
df <- df.gw.sub.v[ , .(V1=E.nodal.wt), by=.(hemi,participant_id, group, threshold)] 
df.l <- df[ hemi=='L',]
df.r <- df[hemi=='R',]
df.lateral <- df.l[,.(participant_id,group, threshold)]
df.lateral$lateral <- (df.l$V1 - df.r$V1) / (df.l$V1 + df.r$V1)
df.lateral[, .(Wp=wilcox.test(lateral~group)$p.value, Tp=t.test(lateral~group)$p.value),by=threshold]

## number of hubs - only different on Frontal & hubs=0
[1] "Frontal"   "Occipital" "Parietal"  "Limbic"    "Insula"    "Temporal" 
df.gw.sub.v <- rbindlist(lapply(func.gw, vertex_attr_dt))
 # only effect on Frontal lobe
df <- df.gw.sub.v[hubs.wt>1,
                  .(region, lobe, hemi, hubs.wt, hubnum=1, E.nodal.wt, participant_id, group, threshold)]

df.gw.sub.v[hubs.wt>1 , .(lobe, hemi, hubs.wt, hubnum=1, E.nodal.wt, participant_id, group, threshold)] %>%
  .[, .(Wp=wilcox.test(E.nodal.wt~group)$p.value, Tp=t.test(E.nodal.wt~group)$p.value), by=threshold]

df[hubs.wt>1,.(Sum = sum(hubnum)), by=.(participant_id, group, threshold)] %>%
  .[, .(Wp=wilcox.test(Sum~group)$p.value, Tp=t.test(Sum~group)$p.value), by=threshold]

df[,.(Sum = sum(hubnum)), by=.(participant_id, group, threshold)] %>%
  .[, .(AUC=trapz(threshold, Sum)), by= .(participant_id,group)] %>%
  .[, .(Wp=wilcox.test(AUC~group)$p.value, Tp=t.test(AUC~group)$p.value)]  

## sub level - hubs nodal efficiency  =====  total efficiency
df.gw.sub.v <- rbindlist(lapply(func.gw, vertex_attr_dt))
df <- df.gw.sub.v[hubs.wt==0 &  hemi=='L',
                  .(region, lobe, hemi, hubs.wt, hubnum=1, E.nodal.wt, participant_id, group, threshold)]

df[, .(HubSum=sum(hubnum), NodSum=sum(E.nodal.wt)),by=.(participant_id, group, threshold)] %>%
  .[, .(HubSum, NodSum, Mean=(NodSum / HubSum)), by=.(participant_id, group, threshold)] %>%
  .[, .(Wp=wilcox.test(NodSum~group)$p.value, Tp=t.test(NodSum~group)$p.value), by=threshold]

df[, .(HubSum=sum(hubnum), NodSum=sum(E.nodal.wt)),by=.(participant_id, group, threshold)] %>%
  .[, .(HubSum, NodSum, Mean=(NodSum / HubSum)), by=.(participant_id, group, threshold)] %>%
  .[,.(AUC=trapz(threshold, NodSum)), by= .(participant_id,group)] %>%
  .[, .(Wp=wilcox.test(AUC~group)$p.value, Tp=t.test(AUC~group)$p.value)]  

### sub level, vertex nodal efficiency
covars.glm <- covars.all[, .(participant_id, group, age, gender, BMI)]
X <- brainGraph_GLM_design(covars.glm, coding = 'effects',factorize = TRUE, binarize = 'gender')
con.mat <- matrix(c(rep(0, 4), -2), nrow = 1, dimnames = list('Control > Patient'))
summary(with(covars.glm,
             brainGraph_GLM(fiber.gw[[1]], measure='E.nodal.wt', covars=covars.glm, X=X, mean.center=TRUE,
                            contrasts = con.mat, alt='greater')))













## MTPC -- error
mtpcVars <- data.table(level='vertex', outcome='E.nodal.wt', alt='greater')
mtpcVars[level=='vertex', N := 5e3]
library(permute)
mtpcPerms <- list(vertex=shuffleSet(n=nrow(covars.glm), nset=mtpcVars[level=='vertex', unique(N)]))
con.mat <- matrix(c(rep(0, 4), -2), nrow = 1, dimnames = list('Control > Patient'))
mtpc.diffs.list <- sapply(mtpcVars[, unique(level)], function(x) NULL)
for (x in names(mtpc.diffs.list)) { 
  mtpc.diffs.list[[x]] <- sapply(mtpcVars[level==x, unique(outcome)], function(x) NULL) 
  for (y in mtpcVars[level==x, outcome]) { 
    print(paste('Level:', x, '; Outcome:', y, ';', format(Sys.time(), '%H:%M:%S'))) 
    mtpc.diffs.list[[x]][[y]] <- 
      mtpc(func.gw, densities, covars=covars.glm, measure=y, contrasts=con.mat, 
           con.type='t', level=x, N=mtpcVars[level==x&outcome==y, N], perms=mtpcPerms[[x]],
           part.method = 'ridgway', alt=mtpcVars[level==x&outcome==y, alt]) 
  } 
}

# data result display
mtpc.diffs.sig.dt <- 
  rbindlist(lapply(mtpc.diffs.list, function(x) 
    rbindlist(lapply(x, function(y) 
      y$DT[A.mtpc > A.crit, .SD[1], by=region]))))

res.mat <- mtpc(fiber.gw, densities, covars=covars.glm, measure=y, contrasts=con.mat, 
     con.type='t', level=x, N=5e3, perms=mtpcPerms[[x]],
     part.method = 'guttman', alt='greater') 



### compute GLM comparison
covars_compare <- df.gw.sub.v[hemi=='L',
                              .(region, hemi, threshold, participant_id,group, hubs.wt, E.nodal.wt)] %>%
  .[,trapz(E.nodal.wt, threshold), by=.(participant_id, group)] %>%
  covars.all[., on=.(participant_id=participant_id)] %>%
  .[, group:=as.factor(group)]
covars_compare$group <- ifelse(covars_compare$group=='patient',1,0)
with(covars_compare, summary(aov(V1~ group + age )) )

####################### statistic - clinical correlation
# indicator
covars.all[, gender:=as.numeric(as.factor(gender))]
patient_df <- df.gw.sub.v[hubs.wt>2 & group=='patient',
                          .(region, hemi, threshold, participant_id,group, hubs.wt, E.nodal.wt)] %>%
  .[,.(AUC=trapz(E.nodal.wt, threshold)), by=.(participant_id)] %>%
  covars.all[., on=.(participant_id=participant_id)]

cor.test(patient_df$AUC, covars.all[group=='patient',]$disease_month, method = 'pearson')
cor.test(patient_df$AUC, covars.all[group=='patient',]$disease_month, method = 'spearman')
cor.test(patient_df$AUC, covars.all[group=='patient',]$`SF-36`, method = 'pearson')
cor.test(patient_df$AUC, covars.all[group=='patient',]$`SF-36`, method = 'spearman')
cor.test(patient_df$AUC, covars.all[group=='patient',]$PRI, method = 'pearson')
cor.test(patient_df$AUC, covars.all[group=='patient',]$PRI, method = 'spearman')
cor.test(patient_df$AUC, covars.all[group=='patient',]$`FS-14`, method = 'pearson')
cor.test(patient_df$AUC, covars.all[group=='patient',]$`FS-14`, method = 'spearman')

pcor.test(patient_df$AUC, covars.all[group=='patient',]$disease_month, 
          covars.all[group=='patient',.(age, gender, BMI)], method = 'pearson')# 0.05
pcor.test(patient_df$AUC, covars.all[group=='patient',]$disease_month, 
          covars.all[group=='patient',.(age, gender, BMI)], method = 'spearman')
pcor.test(patient_df$AUC, covars.all[group=='patient',]$`SF-36`, 
          covars.all[group=='patient',.(age, gender, BMI)], method = 'pearson')
pcor.test(patient_df$AUC, covars.all[group=='patient',]$`SF-36`, 
          covars.all[group=='patient',.(age, gender, BMI)], method = 'spearman')
pcor.test(patient_df$AUC, covars.all[group=='patient',]$PRI, 
          covars.all[group=='patient',.(age, gender, BMI)], method = 'pearson') 
pcor.test(patient_df$AUC, covars.all[group=='patient',]$PRI, 
          covars.all[group=='patient',.(age, gender, BMI)], method = 'spearman')
pcor.test(patient_df$AUC, covars.all[group=='patient',]$`FS-14`, 
          covars.all[group=='patient',.(age, gender, BMI)], method = 'pearson')
pcor.test(patient_df$AUC, covars.all[group=='patient',]$`FS-14`, 
          covars.all[group=='patient',.(age, gender, BMI)], method = 'spearman')



