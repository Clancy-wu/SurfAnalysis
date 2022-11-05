##########################################################################
## Des 2009 模板ROI
library(tidyverse)
library(freesurferformats)
template_file = 'fsaverage5_template/label/'
Des_template.l = read.fs.annot(paste0(template_file, 'lh.aparc.a2009s.annot'))
Des_template.r = read.fs.annot(paste0(template_file, 'rh.aparc.a2009s.annot'))

Des_roi.l = Des_template.l$label_codes #10242个roi
roi_number = unique(Des_roi.l) #76个

Des_roi.r = Des_template.r$label_codes #10242个roi

## 数据框准备
mydata = data.frame()
file_name = c()
## 被试原始数据
library(gifti)
######################### 
subject_l = list.files('Area_dataset/fsaverage5/', pattern = "(hemi-L)",full.names=T)
subject_r = list.files('Area_dataset/fsaverage5/', pattern = "(hemi-R)",full.names=T)
  
# 第i行第k列，左侧
for(i in 1:length(subject_l)){
  #文件名
  file_name[i] = str_extract(subject_l[i], 'sub+\\d++')
  
  for(k in 1:length(roi_number)){
    subject_data = readgii(subject_l[i])[[1]][[1]][ ,1]
    value = mean(subject_data[Des_roi.l==roi_number[k]])
    mydata[i,k] = value
  }
}
# 第i行第j列，右侧
for(i in 1:length(subject_r)){
  for(k in 1:length(roi_number)){
    subject_data = readgii(subject_r[i])[[1]][[1]][ ,1]
    value = mean(subject_data[Des_roi.r==roi_number[k]])
    j = k+length(roi_number)
    mydata[i,j] = value
  }
}
######################### 数据库输出
col_name = c(paste0('L', roi_number), paste0('R', roi_number))
colnames(mydata) = col_name
mydata$file_name = file_name
write.csv(mydata, 'Des_Area_machine_learning.csv',row.names=FALSE)
###################################  END  #######################################
