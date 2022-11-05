###############################################################################################
## Prepare
library(freesurferformats)
library(tidyverse)
library(gifti)
library(coin)

subject_info = read.csv('subject_info.csv', header = T)
subject_info$group = as.factor(subject_info$group)

template_file = 'fsaverage5_template/label/' # template path

# LH indicators
lh.area = 'surf_result/AnatSurfLH/Area/fsaverage5/' # area
lh.curv = 'surf_result/AnatSurfLH/Curv/fsaverage5/' # curv
lh.sulc = 'surf_result/AnatSurfLH/Sulc/fsaverage5/' # sulc
lh.thickness = 'surf_result/AnatSurfLH/Thickness/fsaverage5/' # thickness
lh.volume = 'surf_result/AnatSurfLH/Volume/fsaverage5/' # volume
lh.alff = 'surf_result/FunSurfLH/ALFF_FunSurfWIglobalC/' # alff
lh.falff = 'surf_result/FunSurfLH/fALFF_FunSurfWIglobalC/' # falff
lh.reho = 'surf_result/FunSurfLH/ReHo_FunSurfWIglobalCF/' # reho

# RH indicators
rh.area = 'surf_result/AnatSurfRH/Area/fsaverage5/' # area
rh.curv = 'surf_result/AnatSurfRH/Curv/fsaverage5/' # curv
rh.sulc = 'surf_result/AnatSurfRH/Sulc/fsaverage5/' # sulc
rh.thickness = 'surf_result/AnatSurfRH/Thickness/fsaverage5/' # thickness
rh.volume = 'surf_result/AnatSurfRH/Volume/fsaverage5/' # volume
rh.alff = 'surf_result/FunSurfRH/ALFF_FunSurfWIglobalC/' # alff
rh.falff = 'surf_result/FunSurfRH/fALFF_FunSurfWIglobalC/' # falff
rh.reho = 'surf_result/FunSurfRH/ReHo_FunSurfWIglobalCF/' # reho

# function
roi_indicator = function(sub_indicator_path, roi_number, roi_index){
  # sub_indicator_path: file path that includes .gii files.
  # roi_number: Int, unique value. The order of roi in label table
  # roi_index: Int, 10242 values per vertex (in fsaverage5).
  
  Indicator = list.files(sub_indicator_path, full.names = T)
  nrow = length(Indicator)
  ncol = length(roi_number)
  Matrix = matrix(nrow = nrow, ncol = ncol)
  subject_name = c()
  
  for (i in 1:length(Indicator)){
    subject_name[i] = str_extract(Indicator[i], "sub+\\d++") # subject name
    for (k in 1:length(roi_number)){
      subject_data = readgii(Indicator[i])[[1]][[1]][ ,1]
      subject_value = mean(subject_data[roi_index==roi_number[k]])
      Matrix[i, k] = subject_value
    }
  }
  rownames(Matrix) = subject_name
  colnames(Matrix) = paste0('ROI', roi_number)
  # output matrix
  return(Matrix)
}

roi_comparison = function(Matrix, subject_group, cor_value, roi_name){
  # Matrix: matrix data, output from function roi_indicator.
  # subject_group: factor.
  # fs14: special comparison, if unnecessary, delete it via last part.
  # col_value: for correlation test with ROI
  # roi_name: roi name in label table
  
  # subject_group should be factor
  if(is.factor(subject_group)==FALSE){subject_group = as.factor(subject_group)}
  
  roi_index_2 = c() # similar with roi_index, the short name is in order to escape from repeat.
  roi_name_2 = c() # similar with roi_index, the short name is in order to escape from repeat.
  t_p = c() # test p value
  w_p = c() # wilcox p value
  perm_p = c() # perm p value
  cor_p = c() # cor p value with fs-14
  
  Numbers = colnames(Matrix) # roi_index
  compare_time = dim(Matrix)[2]
  for(i in seq(compare_time)){
    if(sum(Matrix[ ,i]==0)>1 | sum(Matrix[ ,i]=='NaN')>1){
      roi_index_2[i] = Numbers[i]
      roi_name_2[i] = roi_name[i]
      t_p[i] = 'NaN'
      w_p[i] = 'NaN'
      perm_p[i] = 'NaN'
      cor_p[i] = 'NaN'
    }else{
      roi_index_2[i] = Numbers[i]
      roi_name_2[i] = roi_name[i]
      t_p[i] = t.test( Matrix[ ,i]~subject_group, paired=FALSE, var=TRUE)$p.value
      w_p[i] = wilcox.test(Matrix[ ,i]~subject_group, paired=FALSE)$p.value
      perm_p[i] = pvalue(oneway_test(Matrix[ ,i]~subject_group, distribution=approximate(nresample=1000)))[1]
      cor_p[i] = cor.test(Matrix[ ,i], cor_value, method='pearson')$p.value
    }
  }
  
  result_df = data.frame(roi_index=roi_index_2, roi_name=roi_name_2, t_p=t_p, w_p=w_p, perm_p=perm_p, cor_p=cor_p)
  return(result_df)
}
################################################################################################################



###################################  HCP Begin         ###########################################
#####################  HCP Information

### LH info
hcp_template.l = read.fs.annot.gii(paste0(template_file, 'fsaverage5_lh_HCP-MMP1.label.gii'))
hcp_roi_l.index = hcp_template.l$label_codes # 10242 roi index
hcp_roi_l.name =  hcp_template.l$colortable_df$struct_name # 181 roi label name
hcp_roi_l.number = hcp_template.l$colortable_df$code # 181 roi ID

### RH info
hcp_template.r = read.fs.annot.gii(paste0(template_file, 'fsaverage5_rh_HCP-MMP1.label.gii'))
hcp_roi_r.index = hcp_template.r$label_codes # 10242 roi index
hcp_roi_r.name = hcp_template.r$colortable_df$struct_name # 181 roi label name
hcp_roi_r.number = hcp_template.r$colortable_df$code # 181 roi ID

#####################  HCP Comparison

## LH
### Area
lh.area %>%
  roi_indicator(roi_number=hcp_roi_l.number, roi_index=hcp_roi_l.index) %>%
    roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=hcp_roi_l.name) %>%
      write.csv('HCP_LH_Area_Result.csv')

### Curv
lh.curv %>%
  roi_indicator(roi_number=hcp_roi_l.number, roi_index=hcp_roi_l.index) %>%
     roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=hcp_roi_l.name) %>%
        write.csv('HCP_LH_Curv_Result.csv')

### sulc
lh.sulc %>%
  roi_indicator(roi_number=hcp_roi_l.number, roi_index=hcp_roi_l.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=hcp_roi_l.name) %>%
  write.csv('HCP_LH_Sulc_Result.csv')

### thickness
lh.thickness %>%
  roi_indicator(roi_number=hcp_roi_l.number, roi_index=hcp_roi_l.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=hcp_roi_l.name) %>%
  write.csv('HCP_LH_Thickness_Result.csv')

### volume
lh.volume %>%
  roi_indicator(roi_number=hcp_roi_l.number, roi_index=hcp_roi_l.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=hcp_roi_l.name) %>%
  write.csv('HCP_LH_Volume_Result.csv')

### alff
lh.alff %>%
  roi_indicator(roi_number=hcp_roi_l.number, roi_index=hcp_roi_l.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=hcp_roi_l.name) %>%
  write.csv('HCP_LH_Alff_Result.csv')

### falff
lh.falff %>%
  roi_indicator(roi_number=hcp_roi_l.number, roi_index=hcp_roi_l.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=hcp_roi_l.name) %>%
  write.csv('HCP_LH_FAlff_Result.csv')

### reho
lh.reho %>%
  roi_indicator(roi_number=hcp_roi_l.number, roi_index=hcp_roi_l.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=hcp_roi_l.name) %>%
  write.csv('HCP_LH_Reho_Result.csv')


## RH
### Area
rh.area %>%
  roi_indicator(roi_number=hcp_roi_r.number, roi_index=hcp_roi_r.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=hcp_roi_r.name) %>%
  write.csv('HCP_RH_Area_Result.csv')

### Curv
rh.curv %>%
  roi_indicator(roi_number=hcp_roi_r.number, roi_index=hcp_roi_r.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=hcp_roi_r.name) %>%
  write.csv('HCP_RH_Curv_Result.csv')

### sulc
rh.sulc %>%
  roi_indicator(roi_number=hcp_roi_r.number, roi_index=hcp_roi_r.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=hcp_roi_r.name) %>%
  write.csv('HCP_RH_Sulc_Result.csv')

### thickness
rh.thickness %>%
  roi_indicator(roi_number=hcp_roi_r.number, roi_index=hcp_roi_r.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=hcp_roi_r.name) %>%
  write.csv('HCP_RH_Thickness_Result.csv')

### volume
rh.volume %>%
  roi_indicator(roi_number=hcp_roi_r.number, roi_index=hcp_roi_r.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=hcp_roi_r.name) %>%
  write.csv('HCP_RH_Volume_Result.csv')

### alff
rh.alff %>%
  roi_indicator(roi_number=hcp_roi_r.number, roi_index=hcp_roi_r.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=hcp_roi_r.name) %>%
  write.csv('HCP_RH_Alff_Result.csv')

### falff
rh.falff %>%
  roi_indicator(roi_number=hcp_roi_r.number, roi_index=hcp_roi_r.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=hcp_roi_r.name) %>%
  write.csv('HCP_RH_FAlff_Result.csv')

### reho
rh.reho %>%
  roi_indicator(roi_number=hcp_roi_r.number, roi_index=hcp_roi_r.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=hcp_roi_r.name) %>%
  write.csv('HCP_RH_Reho_Result.csv')


###################################  HCP finished         ###########################################


####################################    Des Begin     ##########################################
#####################  Destrieux 2009 atlas Information

### LH info
des_template.l = read.fs.annot(paste0(template_file, 'lh.aparc.a2009s.annot'))
des_roi_l.index = des_template.l$label_codes # 10242 roi index
des_roi_l.name =  des_template.l$colortable_df$struct_name # 76 roi label name
des_roi_l.number = des_template.l$colortable_df$code # 76 roi ID

### RH info
des_template.r = read.fs.annot(paste0(template_file, 'rh.aparc.a2009s.annot'))
des_roi_r.index = des_template.r$label_codes # 10242 roi index
des_roi_r.name = des_template.r$colortable_df$struct_name # 76 roi label name
des_roi_r.number = des_template.r$colortable_df$code # 76 roi ID

#####################  Destrieux 2009 atlas Comparison

## LH
### Area
lh.area %>%
  roi_indicator(roi_number=des_roi_l.number, roi_index=des_roi_l.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=des_roi_l.name) %>%
  write.csv('des_LH_Area_Result.csv')

### Curv
lh.curv %>%
  roi_indicator(roi_number=des_roi_l.number, roi_index=des_roi_l.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=des_roi_l.name) %>%
  write.csv('des_LH_Curv_Result.csv')

### sulc
lh.sulc %>%
  roi_indicator(roi_number=des_roi_l.number, roi_index=des_roi_l.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=des_roi_l.name) %>%
  write.csv('des_LH_Sulc_Result.csv')

### thickness
lh.thickness %>%
  roi_indicator(roi_number=des_roi_l.number, roi_index=des_roi_l.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=des_roi_l.name) %>%
  write.csv('des_LH_Thickness_Result.csv')

### volume
lh.volume %>%
  roi_indicator(roi_number=des_roi_l.number, roi_index=des_roi_l.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=des_roi_l.name) %>%
  write.csv('des_LH_Volume_Result.csv')

### alff
lh.alff %>%
  roi_indicator(roi_number=des_roi_l.number, roi_index=des_roi_l.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=des_roi_l.name) %>%
  write.csv('des_LH_Alff_Result.csv')
### falff
lh.falff %>%
  roi_indicator(roi_number=des_roi_l.number, roi_index=des_roi_l.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=des_roi_l.name) %>%
  write.csv('des_LH_FAlff_Result.csv')
### reho
lh.reho %>%
  roi_indicator(roi_number=des_roi_l.number, roi_index=des_roi_l.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=des_roi_l.name) %>%
  write.csv('des_LH_Reho_Result.csv')

## RH
### Area
rh.area %>%
  roi_indicator(roi_number=des_roi_r.number, roi_index=des_roi_r.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=des_roi_r.name) %>%
  write.csv('des_RH_Area_Result.csv')

### Curv
rh.curv %>%
  roi_indicator(roi_number=des_roi_r.number, roi_index=des_roi_r.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=des_roi_r.name) %>%
  write.csv('des_RH_Curv_Result.csv')

### sulc
rh.sulc %>%
  roi_indicator(roi_number=des_roi_r.number, roi_index=des_roi_r.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=des_roi_r.name) %>%
  write.csv('des_RH_Sulc_Result.csv')

### thickness
rh.thickness %>%
  roi_indicator(roi_number=des_roi_r.number, roi_index=des_roi_r.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=des_roi_r.name) %>%
  write.csv('des_RH_Thickness_Result.csv')

### volume
rh.volume %>%
  roi_indicator(roi_number=des_roi_r.number, roi_index=des_roi_r.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=des_roi_r.name) %>%
  write.csv('des_RH_Volume_Result.csv')

### alff
rh.alff %>%
  roi_indicator(roi_number=des_roi_r.number, roi_index=des_roi_r.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=des_roi_r.name) %>%
  write.csv('des_RH_Alff_Result.csv')

### falff
rh.falff %>%
  roi_indicator(roi_number=des_roi_r.number, roi_index=des_roi_r.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=des_roi_r.name) %>%
  write.csv('des_RH_FAlff_Result.csv')

### reho
rh.reho %>%
  roi_indicator(roi_number=des_roi_r.number, roi_index=des_roi_r.index) %>%
  roi_comparison(subject_group=subject_info$group, cor_value=subject_info$FS14, roi_name=des_roi_r.name) %>%
  write.csv('des_RH_Reho_Result.csv')

###################################  Des finished         ###########################################

