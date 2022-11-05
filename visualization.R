library(freesurferformats)

##########################################################################
## HCP模板ROI
hcp_template_L = read.fs.annot.gii('fsaverage5_template/label/fsaverage5_lh_HCP-MMP1.label.gii')
hcp_roi_L = hcp_template_L$label_codes #10242个roi

hcp_template_R = read.fs.annot.gii('fsaverage5_template/label/fsaverage5_rh_HCP-MMP1.label.gii')
hcp_roi_R = hcp_template_R$label_codes # 10242个roi

##########################################################################
## data input 模板
input_template_L = read.fs.mgh('data_input_L.mgz') # 左侧模板，包含各种格式
data_input_L = input_template_L[ ,1,1,1] #提取数值, 10242个数值
data_input_L = data_input_L * 0  # 左侧模板初始化

input_template_R = read.fs.mgh('data_input_R.mgz') #右侧模板
data_input_R = input_template_R[ ,1,1,1] # 10242个数值
data_input_R = data_input_R * 0  # 右侧模板初始化

##########################################################################
## 数据写入
library(readxl)
library(tidyverse)
result_data = read_excel('hcp_withz_roi_p.xlsX', sheet = 'target')
for (i in 1:length(result_data$roi_name)){
  result_data.hemi = str_extract(result_data$roi_name[i], '[:upper:]') #提取半球
  result_data.roi = str_extract(result_data$roi_name[i], '[:digit:]')  #提取roi
  if(result_data.hemi == 'L'){
    data_input_L[hcp_roi_L == result_data.roi] = result_data$w_p[i] # 只输入wilcox P
  }else if(result_data.hemi == 'R'){
    data_input_R[hcp_roi_R == result_data.roi] = result_data$w_p[i]
  }else{next}  # 用L和R为了增加一层保险，保证结果可靠
}

#########################  模板输出可视化
write.fs.mgh('hcp_sulc_L.mgz', data_input_L)
write.fs.mgh('hcp_sulc_R.mgz', data_input_R)
###################################  END  #######################################



library(ggseg)
