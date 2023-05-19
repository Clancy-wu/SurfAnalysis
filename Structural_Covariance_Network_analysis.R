library(brainGraph)
library(data.table)
options(bg.subject_id='participant_id', bg.group='group')
data.l.file <- fread('aparc.a2009s_lh_thickness.csv')
data.r.file <- fread('aparc.a2009s_rh_thickness.csv')
data.l <- data.l.file[,2:75]; data.r <- data.r.file[, 2:75]
lhrh <- cbind(data.l, data.r)
# sort name
library(stringr)
oldname <- colnames(lhrh)
newname <- str_replace_all(oldname, c('h_'='', '&'='_and_', '_thickness'='', '-'='.'))
colnames(lhrh) <- newname
# begin
covars.all <- fread('subjects_71_info.csv')
covars.all[, gender := as.factor(gender)]
covars <- covars.all[,c('age', 'gender', 'group')]
covars.id <- paste0('sub-', covars.all$participant_id)
covars$participant_id <- covars.id
lhrh$participant_id <- covars.id

myResids <- get.resid(lhrh, covars, atlas = 'destrieux')
densities <- seq(0.01, 0.34, 0.01) 
corrs <- corr.matrix(myResids, densities=densities)

g <- lapply(seq_along(densities), function(x) 
  make_brainGraphList(corrs[x], modality='thickness'))
dt.G <- rbindlist(lapply(g, graph_attr_dt))  
dt.V <- rbindlist(lapply(g, vertex_attr_dt))

plot(myResids, region='lS_calcarine', cols = T)[[1]]

setkey(dt.V, density, group)
dt.V[density < 0.1, .SD[which.max(degree), .(region, degree)], 
     by=.(group, density)]

library(corrplot)
corrplot(corrs$R[,,2], method='square', tl.pos="n")

library(coin)
dt.G[,group:=as.factor(group)]
coin::pvalue(oneway_test(dt.G$assort~ dt.G$group, distribution=approximate(nresample=10000)))

library(ggplot2)
ggplot(dt.G, aes(x=density, y=max.comp, colour =group))+
  geom_point()
######################################
library(freesurferformats)
Template.file <- read.nifti1.data('aparc.a2009s+aseg.nii.gz')
Template.header <- read.nifti1.header('aparc.a2009s+aseg.nii.gz')
mydata <- Template.file
#Colormap <- read.table('FreeSurferColorLUT.txt',  comment.char='#')
mydata[mydata<300]=0; mydata[mydata==11142]=0; mydata[mydata==12142]=0;  
############################################################################
library(fsbrain)

template_subject = 'fsaverage'
template_subject_dir = 'C:/Users/Clancy/Desktop/Structure_Dataset_71/freesurfer/'
A <- vis.subject.annot(subjects_dir = template_subject_dir, 
                  subject_id = 'fsaverage', atlas = 'aparc.a2009s', hemi='both',
                  surface='pial', views = c('t4') )
#export(A,output_img='fig1.tiff')

morph_data_both = runif(327684, min=-1, max=1)
B=vis.data.on.fsaverage(subjects_dir = template_subject_dir, 
                      vis_subject_id = 'fsaverage', 
                      morph_data_both = morph_data_both, 
                      draw_colorbar = F,
                      surface = 'pial',
                      makecmap_options=list('colFn'=squash::jet))
export(B,output_img='B_test.tiff', colorbar_legend = 'B TEST')

############################################################################
avpath <- 'C:/Users/Clancy/Desktop/Structure_Dataset_71/freesurfer/'
lh_annot <- subject.annot(avpath, "fsaverage", "lh", "aparc.a2009s")
rh_annot <- subject.annot(avpath, "fsaverage", "rh", "aparc.a2009s")
lh_inf <- subject.surface(avpath, "fsaverage", "inflated", "lh")
rh_inf <- subject.surface(avpath, "fsaverage", "inflated", "rh")
aparc_outline <- list(
  lh = fsbrain::annot.outline(annotdata = tempalte.l.file,
                              surface_mesh = lh_inf,
                              outline_color = NULL),
  rh = fsbrain::annot.outline(annotdata = tempalte.r.file,
                              surface_mesh = rh_inf,
                              outline_color = NULL)
)
############################################################################
library(stringr)
corrs$R[,,1][,1] #health
corrs$R[,,2][,1] #patient

template.l.file <- read.fs.annot('C:/Users/Clancy/Desktop/Structure_Dataset_71/freesurfer/fsaverage/label/lh.aparc.a2009s.annot')
template.r.file <- read.fs.annot('C:/Users/Clancy/Desktop/Structure_Dataset_71/freesurfer/fsaverage/label/rh.aparc.a2009s.annot')
data_input.l <- data_input.r <- rep(0, length(tempalte.l.file$vertices))
for (i in unique(tempalte.l.file$label_names)){
  if (i == 'Unknown' | i == 'Medial_wall') {  
    Mask.l <- tempalte.l.file$label_names == i
    Mask.r <- tempalte.r.file$label_names == i
    data_input.l[Mask.l] <- Value.l <- NA
    data_input.r[Mask.r] <- Value.r <- NA
  }else{
  # i is the name of freesurfer
  Mask.l <- tempalte.l.file$label_names == i
  ifelse(str_detect(i,'-'), 
         bG.name.l <- paste0('l', str_replace_all(i, c('-'='.'))),
         bG.name.l <- paste0('l',i))
  Value.l <- corrs$R[,,2][,1][[bG.name.l]]
  data_input.l[Mask.l] <- Value.l
  # right
  Mask.r <- tempalte.r.file$label_names == i
  ifelse(str_detect(i,'-'), 
         bG.name.r <- paste0('r', str_replace_all(i, c('-'='.'))),
         bG.name.r <- paste0('r',i))
  Value.r <- corrs$R[,,2][,1][[bG.name.r]]
  data_input.r[Mask.r] <- Value.r
  }
}
# fix left bug
data_input.l[102162:102163] <- NA

fsbrain.set.default.figsize(1200, 1200);
#rgla = list('trans_fun'=limit_fun_na(-1,1));
RdBu = colorRampPalette(rev(RColorBrewer::brewer.pal(5, name='RdBu')));
makecmap_options = list('colFn'=RdBu, 'symm'=TRUE);
vis.data.on.fsaverage(subjects_dir = template_subject_dir, 
                      vis_subject_id = 'fsaverage', 
                      morph_data_lh = data_input.l,
                      morph_data_rh = data_input.r,
                      draw_colorbar = T,
                      surface = 'inflated',
                      bg='curv',
                      makecmap_options=makecmap_options ) 
export(output_img='patient_test.png', colorbar_legend = 'Patient group')
############################################################################





