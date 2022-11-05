library(freesurferformats)
hcp.l = read.fs.annot.gii('fsaverage5_template/label/fsaverage5_lh_HCP-MMP1.label.gii')
code = hcp.l$colortable_df$code
name = hcp.l$colortable_df$struct_name
mydata = data.frame(code=code, name=name)
write.csv(mydata, 'HCP_L.csv')


hcp.r = read.fs.annot.gii('fsaverage5_template/label/fsaverage5_rh_HCP-MMP1.label.gii')
code = hcp.r$colortable_df$code
name = hcp.r$colortable_df$struct_name
mydata = data.frame(code=code, name=name)
write.csv(mydata, 'HCP_R.csv')

des.l = read.fs.annot('fsaverage5_template/label/lh.aparc.a2009s.annot')
code = des.l$colortable_df$code
name = des.l$colortable_df$struct_name
mydata = data.frame(code=code, name=name)
write.csv(mydata, 'Des_L.csv')


des.r = read.fs.annot('fsaverage5_template/label/rh.aparc.a2009s.annot')
code = des.r$colortable_df$code
name = des.r$colortable_df$struct_name
mydata = data.frame(code=code, name=name)
write.csv(mydata, 'Des_R.csv')