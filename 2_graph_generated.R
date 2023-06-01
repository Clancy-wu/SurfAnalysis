# 2. graph generated
## Annotation: gb=group binary, gw=group weighted, 
##             sb=subject binary, sw=subject weighted
library(brainGraph)
library(data.table)
library(doMC)
registerDoMC(detectCores()) # define the number of core
options(bg.subject_id='participant_id', bg.group='group')
grps = c('health', 'patient')
############## multiple modal network
############## graph generated
densities <- seq(0.01, 0.34, 0.01)
setkey(destrieux, index)
covars.all <- fread('subjects_71_info.csv')
covars.all[, gender := as.factor(gender)]
inds = lapply(grps, function(x) covars.all[group == x, which = TRUE])
atlas = 'destrieux'
if(!dir.exists('brainGraphRDS')) dir.create('brainGraphRDS')
######################################################################
## Structural Covariance Network
df.anat.file <- fread('NetAnat/aparc.a2009s_both_thickness.csv')
lhrh <- df.anat.file[,2:149]
library(stringr) # sort name
oldname <- colnames(lhrh)
colnames(lhrh) <- str_replace_all(oldname, c('h_'='', '-'='.'))
covars <- covars.all[,c('age', 'gender', 'group')]
covars.id <- paste0('sub-', covars.all$participant_id)
lhrh$participant_id <- covars$participant_id <- covars.id
myResids <- get.resid(lhrh, covars, atlas = 'destrieux')
source('corr_matrix.R') # self-define, add function weighted = F (default) or T
## weighted , weighted=T in corr.matrix
corrs <- corr.matrix(myResids, densities=densities, weighted=T, what='resids', type='pearson')
anat.gw.group <- lapply(seq_along(densities), function(x) 
              make_brainGraphList(corrs[x], modality='thickness', level='group', 
                                  weighted = TRUE, grpNames = grps))
saveRDS(anat.gw.group, file=file.path('brainGraphRDS', 'anat.gw.group.rds'), compress = 'xz')
## binary , weighted=F in corr.matrix
corrs <- corr.matrix(myResids, densities=densities, weighted=F, what='resids', type='pearson')
anat.gb.group <- lapply(seq_along(densities), function(x) 
  make_brainGraphList(corrs[x], modality='thickness', level='group', 
                      weighted = TRUE, grpNames = grps))
saveRDS(anat.gb.group, file=file.path('brainGraphRDS', 'anat.gb.group.rds'), compress = 'xz')

######################################################################
## Functional Network
matfiles <- list.files('NetFun/pos', pattern='sub-sub', full.names = T)
my.mats <- create_mats(matfiles, modality = 'fmri',threshold.by = 'density',
                       mat.thresh = densities, inds = inds)
## weighted
A.norm.sub <- my.mats$A.norm.sub; A.norm.mean <- my.mats$A.norm.mean
func.gw <- func.gw.group <- vector('list', length(densities))
for (j in seq_along(densities)) {
  func.gw[[j]] <- make_brainGraphList(A.norm.sub[[j]], atlas, level='subject',
                                modality = 'fmri', threshold = densities[j],
                                weighted = TRUE, gnames = covars.all$participant_id,
                                grpNames = covars.all$group )
  func.gw.group[[j]] <- make_brainGraphList(A.norm.mean[[j]], atlas, level='group',
                                 modality = 'fmri',threshold = densities[j],
                                 weighted = TRUE, grpNames = grps )
}
saveRDS(func.gw, file=file.path('brainGraphRDS', 'func.gw.rds'), compress = 'xz')
saveRDS(func.gw.group, file=file.path('brainGraphRDS', 'func.gw.group.rds'), compress = 'xz')

## binary
A.norm.sub <- my.mats$A.norm.sub; A.norm.mean <- my.mats$A.norm.mean
for (i in seq_along(densities)){ A.norm.sub[[i]][A.norm.sub[[i]] > 0 ] = 1 }
for (i in seq_along(densities)){ A.norm.mean[[i]][A.norm.mean[[i]] > 0 ] = 1 }
func.gb <- func.gb.group <- vector('list', length(densities))
for (j in seq_along(densities)) {
  func.gb[[j]] <- make_brainGraphList(A.norm.sub[[j]], atlas, level='subject',
                                modality = 'fmri', threshold = densities[j],
                                weighted = TRUE, gnames = covars.all$participant_id,
                                grpNames = covars.all$group )
  func.gb.group[[j]] <- make_brainGraphList(A.norm.mean[[j]], atlas, level='group',
                                      modality = 'fmri',threshold = densities[j],
                                      weighted = TRUE, grpNames = grps )
}
saveRDS(func.gb, file=file.path('brainGraphRDS', 'func.gb.rds'), compress = 'xz')
saveRDS(func.gb.group, file=file.path('brainGraphRDS', 'func.gb.group.rds'), compress = 'xz')    

######################################################################
## Fiber Network
matfiles <- list.files('NetFiber/', pattern='sub-sub', full.names = T)
my.mats <- create_mats(matfiles, modality = 'dti',threshold.by = 'density',
                       mat.thresh = densities, inds = inds)
## weighted
A.norm.sub <- my.mats$A.norm.sub; A.norm.mean <- my.mats$A.norm.mean
fiber.gw <- fiber.gw.group <- vector('list', length(densities))
for (j in seq_along(densities)) {
  fiber.gw[[j]] <- make_brainGraphList(A.norm.sub[[j]], atlas, level='subject',
                                modality = 'dti', threshold = densities[j],
                                weighted = TRUE, gnames = covars.all$participant_id,
                                grpNames = covars.all$group )
  fiber.gw.group[[j]] <- make_brainGraphList(A.norm.mean[[j]], atlas, level='group',
                                      modality = 'dti',threshold = densities[j],
                                      weighted = TRUE, grpNames = grps )
}
saveRDS(fiber.gw, file=file.path('brainGraphRDS', 'fiber.gw.rds'), compress = 'xz')
saveRDS(fiber.gw.group, file=file.path('brainGraphRDS', 'fiber.gw.group.rds'), compress = 'xz')  

## binary
A.norm.sub <- my.mats$A.norm.sub; A.norm.mean <- my.mats$A.norm.mean
for (i in seq_along(densities)){ A.norm.sub[[i]][A.norm.sub[[i]] > 0 ] = 1 }
for (i in seq_along(densities)){ A.norm.mean[[i]][A.norm.mean[[i]] > 0 ] = 1 }
fiber.gb <- fiber.gb.group <- vector('list', length(densities))
for (j in seq_along(densities)) {
  fiber.gb[[j]] <- make_brainGraphList(A.norm.sub[[j]], atlas, level='subject',
                                modality = 'dti', threshold = densities[j],
                                weighted = TRUE, gnames = covars.all$participant_id,
                                grpNames = covars.all$group )
  fiber.gb.group[[j]] <- make_brainGraphList(A.norm.mean[[j]], atlas, level='group',
                                      modality = 'dti',threshold = densities[j],
                                      weighted = TRUE, grpNames = grps )
}
saveRDS(fiber.gb, file=file.path('brainGraphRDS', 'fiber.gb.rds'), compress = 'xz')
saveRDS(fiber.gb.group, file=file.path('brainGraphRDS', 'fiber.gb.group.rds'), compress = 'xz')  

######################################################################
## END
######################################################################
