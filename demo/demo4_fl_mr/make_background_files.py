import os
import glob
import pandas as pd

folder='out_data_gencode/out_dmr_vs_deg/out4mr_not_in_tss_enhancer/'
file_prefix='data/*raw*.*'
out_file_str='background_samples_list.tsv'
out_folder='out_data_gencode/out_dmr_vs_deg/'

record_files=[]
record_df=[]
for i in range(1, 25):
  if i==23:
     i='X'
  elif i==24:
     i='Y'
  else:
     i=str(i)
  #
  tmp_files=glob.glob(os.path.join(folder,'chr'+i, file_prefix))
  tmp_out_file=os.path.join(out_folder, out_file_str+'_chr'+i+ '.tsv')
  tmp_df=pd.DataFrame(data=tmp_files,columns=['files'])
  #print(tmp_out_file) 
  record_files.append(tmp_out_file)
  #tmp_df.to_csv(tmp_out_file,sep='\t',index=False, header=None)
  record_df.append(tmp_df.copy())

out_file= os.path.join(out_folder, out_file_str)
print(out_file)
out_df=pd.concat(record_df,axis=0).copy()
out_df.to_csv(out_file,sep='\t',header=None, index=False)




