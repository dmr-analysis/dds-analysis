#this script is used to find DEG genes in TSS and 5distnace regions, then export all of them with a bed format such
#as chr, start, end, mr_ID, Pvalues et al.
import pandas as pd
import numpy as np
import os
from scipy.stats import zscore
#exec(open("find_DEG_in_tss_5dist_regions_hap1.py").read())

###-- INPUT parameters which may need be changed in differen runs-- #####
#input file folder for MRs overlapping to TSS and 5Distance regions that was done in dmr_analysis
in_folder='../../data/hap1_cell/in_data/final_demo_data/hap1_cell/out_data/DMR_CpG_context/out_map2genome/'
in_data_str='_hap1'

#file name of TSS regions that overlapping with DMRs that located in in_folder
tss_file='3_chroms_all_mr_data_range_dmrRanking_TSS_Up5000_Down1000_removedShort_overlap1e-09.bed'
in_tss_file=os.path.join( in_folder , tss_file)
print('Read DMR file in TSS regions: ')
print(in_tss_file,'\n')

#file name of 5'distance regions that overlapping with DMRs that located in in_folder
dist_file='3_chroms_all_mr_data_range_dmrRanking_noGenes_5dist_Up1000000_Up5000removedShort_overlap1e-09.bed'
in_5dist_file=os.path.join( in_folder , dist_file)
print('Read DMR file in 5distance regions: ')
print(in_5dist_file,'\n')



#DEG: differetial expression gene file in tab delimited format exported by bpb3
in_deg_file='../../data/hap1_cell/in_data/final_demo_data/hap1_cell/in_data/DEG/HAP1_P1_vs_HAP1_KO1_differentially_expressed_genes_min1.1Fd_min1RPKM.txt'
print('Read DEG file: ')
print(in_deg_file,'\n')

###---OUTPUT parameters -- ######
#output file path ?
out_result_folder='../../data/hap1_cell/out_data/'

###---START to Run --- #######
#read all data such as DMR-in-TSS, DMR-in-5dist, and DEG files
tss_df=pd.read_csv(in_tss_file,sep='\t',header=None)
dist_df=pd.read_csv(in_5dist_file,sep='\t',header=None)
deg_df=pd.read_csv(in_deg_file,sep='\t')
#deg_df=out_df.copy()

selected_col_genename='#gene'

#combine DMR-in-TSS, DMR-in-5distance and DEG data to a singl data frame
combined_df=pd.concat([tss_df,dist_df]).copy()
combined_df.reset_index(inplace=True,drop=True)
combined_df.columns=['chrom','start_pos','end_pos','gene_type','mr_chrom','mr_start_pos','mr_end_pos','mr_id','pval']

#find overlappping between DMR-in-TSS or DMR-in-5istance and DEG , here 'D' means predicted DMR
#or user can use other criteria to select DMRs from full list.
dmr_combined_df=combined_df[combined_df.mr_id.apply(lambda x: ':D' in x)].copy()
dmr_combined_df.reset_index(inplace=True, drop=True)

#combine DMR with DEG
deg_genes=set(deg_df[selected_col_genename].to_list())
dmr_combined_df['gene_name']=dmr_combined_df.gene_type.apply(lambda x: x.split('||')[2].split(':')[0] )

deg_dmr_combined_df=dmr_combined_df[dmr_combined_df.gene_name.apply(lambda x : x in deg_genes)].copy()
deg_dmr_combined_df.reset_index(inplace=True,drop=True)

#export DMR regions are overlapping with either TSS or 5dist of DEGs
out_df=deg_dmr_combined_df[['mr_chrom','mr_start_pos','mr_end_pos','mr_id','pval','gene_type','gene_name']]
out_file=os.path.join(out_result_folder, 'dmr_regions_in_deg_tss_5dist'+in_data_str+'.bed')
out_df.to_csv(out_file,sep='\t',header=None,index=False)
print('Export combined DMR-in-TSS or DMR-in-5distnace with DEG to file:')
print('\t ', out_file)


