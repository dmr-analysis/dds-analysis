#this script is used to add chr before the chromsome number
import pandas as pd

#file1='Table1_common_boundaries_merged_sorted.bed'
file1='Table4_TAD_annotations_sorted.bed'
in_df=pd.read_csv(file1,sep='\t',header=None)
in_df[0]=in_df[0].apply(lambda x: 'chr'+x)

out_file=file1.replace('.bed', '_chr.bed')
in_df.to_csv(out_file,sep='\t', index=False, header=None)




