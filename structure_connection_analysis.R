library(readr)
library(tidyverse)
library(coin)

###############################################################################
structure_connect_file = 'structure_connection/'

subject_connect = list.files(structure_connect_file, full.names = T)
subject_info = read.csv('subject_info.csv', header = T)
subject_info$group = as.factor(subject_info$group) # group

Destrieux_cortex.l = seq(1,74)
cerebellum_cortex.l = c(75)
Destrieux_subcortex.l = seq(76,82)
Destrieux_subcortex.r = seq(83,89)
Destrieux_cortex.r = seq(90,163)
cerebellum_cortex.r = c(164)

###############################################################################
## function defination
all_roi_extract = function(subject_file, roi_list){
  # subject_file_path: the subject path you want to extract.
  # roi_list: Series, the interested roi index you want to extact.
  
  # create extract_matrix to save data
  nrow = length(subject_file)
  ncol = length(roi_list) ** 2
  extract_matrix = matrix(nrow = nrow, ncol = ncol)
  subject_name = c()
  
  # i is the subject
  for(i in 1:length(subject_file)){
    subject_name[i] = str_extract(subject_file[i], "sub+\\d++") # subject name
    subject.org = read_csv(subject_file[i], col_names = F) # read data
    subject.interest = subject.org[roi_list, roi_list]
    #with t, col by col; no t, row by row
    extract_matrix[i, ]  = as.vector(t(subject.interest))
  }
  
  # output extract_matrix
  rownames(extract_matrix) = subject_name
  colnames(extract_matrix) = paste0('ROI', seq(ncol))
  return(extract_matrix)
}

roi_extract_compare = function(Matrix, subject_group, cor_value){
  # Matrix: matrix format file
  # subject_group: subject group
  # col_value: for correlation test with ROI
  
  # subject_group should be factor
  if(is.factor(subject_group)==FALSE){subject_group = as.factor(subject_group)}
  
  # result contain
  roi_name = c()
  t_p = c()
  w_p = c()
  perm_p = c()
  cor_p = c()
  
  # compute
  Subjects = length(subject_group)
  Names = colnames(Matrix)
  compare_time = dim(Matrix)[2]
  for(i in seq(compare_time)){
    # if the roi equals 0 in all subjects, then put 0
    if(sum(Matrix[ ,i]==0)==Subjects | sum(Matrix[ ,i]=="NaN")==Subjects){
      roi_name[i] = Names[i]
      t_p[i] = "NaN"
      w_p[i] = "NaN"
      perm_p[i] = "NaN"
      cor_p[i] = "NaN"
    }else{
      roi_name[i] = Names[i]
      t_p[i] = t.test(Matrix[ ,i]~subject_group, paired=F, var=T)$p.value
      w_p[i] = wilcox.test(Matrix[ ,i]~subject_group, paired=F)$p.value
      perm_p[i] = pvalue(oneway_test(Matrix[ ,i]~subject_group, distribution=approximate(nresample=1000)))[1]
      cor_p[i] = cor.test(Matrix[ ,i], cor_value, method='pearson')$p.value
    }
  }
  
  ## output
  roi_length = sqrt(compare_time)
  roi_each = paste0("ROI", seq(roi_length))
  roi1 = rep(roi_each, each=roi_length)
  roi2 = rep(roi_each, roi_length)
  result_df = data.frame(roi_name=roi_name, roi1=roi1, roi2=roi2, t_p=t_p, w_p=w_p, perm_p=perm_p, cor_p=cor_p)
  return(result_df)
}

###############################################################################
## begin
interest_roi = c(Destrieux_cortex.l, Destrieux_cortex.r)
subject_connect %>%
  all_roi_extract(roi_list=interest_roi) %>%
  roi_extract_compare(subject_group=subject_info$group,cor_value=subject_info$FS14) %>%
  write.csv('Des_128roi_compare_result.csv')
