#this script is used to generate data for DMR vs. expressiondata
import os
import pandas as pd
import glob 

#this script is used to prepare for files before running dds_Analysis dTarget_methy_vs_express 

#13.
#1. Path for results of dmr_analysis
mr_data_folder='../../data/rat_data/in_data/final_demo_data/rat_data/out_data/DMR_CpG_context/'
#file name of all DMRs in the result folder of dmr_analysis
full_mr_file=os.path.join( mr_data_folder, '5_chroms_all_mr_data_range_dmrRanking.bed')
print('Read ranked MRs from results of dmr_analysis :')
print(full_mr_file,'\n')

#2. output path for exported gene lists files that contain DMR overlapping with TSS or enhancers regions
out_result_folder='../../data/rat_data/out_data/'

#3. genome file path which contains refFlat file and enhancer files
in_genome_folder='../../data/rat_data/in_data//final_demo_data/genome'
in_genome_file='rn6/rn6.enhancers_all_rn5_merged_rn6liftOvered_4dmr.bed'

#read enhancer regions file from genome folder
enhancer_file=os.path.join( in_genome_folder , in_genome_file)
#read TSS regions file from results of dmr_analysis
tss_file=os.path.join(mr_data_folder , 'data/TSS_Up5000_Down1000_removedShort.bed')
print('Read TSS region file: ')
print(tss_file,'\n')

#4. to find mr not located in either enhacner or tss regions that will be used to build a list of background MRs
# To build a list of MRs for randomly selected background MRs , we first need to remove MRs that are located in either TSS or enhancers regins
out_file=os.path.join(out_result_folder, 'mr_regions_not_in_enhancers.bed')
os.system('bedtools intersect -a ' + full_mr_file + ' -b ' + enhancer_file + \
                  ' -v > ' + out_file)
out_file2mr_not_in_enhancer_tss=out_file.replace('.bed','_tss.bed')
os.system('bedtools intersect -a ' + out_file + ' -b ' + tss_file + \
          ' -v > ' + out_file2mr_not_in_enhancer_tss )
os.system('rm -f ' + out_file)
print('Export MR regions ddo not locate in TSS or Enhancer regions:')
print('\t', out_file2mr_not_in_enhancer_tss,'\n')

#5. input DMR regions that overlapping with TSS and 5distance regions of DEGs that export by find_DEG_in_tss_5dist_regions_test.py"
in_file=os.path.join(out_result_folder, 'dmr_regions_in_deg_tss_5dist_rat.bed')
in_df=pd.read_csv(in_file,sep='\t',header=None)
print('Read DMRs that are overlapping with either TSS or 5distance regions of DEGs : ')
print(in_file,'\n')

#14a.
#6. to find unique dmr located in tss and 5dist regions of DEG
#here assue columns 0,1,2,3 are chr, star, end, and name, and 5 is the genome type
out_df=in_df.groupby([0,1,2,3])[5].apply(';'.join).reset_index().copy()
out_file=in_file.replace('dmr_','uqdmr_')
out_df.to_csv(out_file,sep='\t',header=False,index=False)
print('Export unique DMR regions that overlapping to TSS or 5distance regions of DEG: ')
print('\t',out_file,'\n')

#7. to continue find unique gene dmrs in tss and 5dist regions of DEG
#column 3 is DMR id, column 5 is the region/gene type such as TSS, 5dist et al
in_df['mr_site']=in_df[3].apply(lambda x: ':'.join(x.split(':')[0:2]))
in_df['gene_type']=in_df[5].apply(lambda x: x.split('||')[1])

out_df2=in_df.groupby([6,'gene_type'])['mr_site'].apply('~'.join).reset_index().copy()
out_df2.columns=['gene_name','gene_type','new_mr_sites']
out_file2=in_file.replace('dmr_','uqGeneDmr_')
out_df2.to_csv(out_file2,sep='\t',index=False)
print('Export unique genes in DMRs: ')

#8. to find tss overlapping with DMRs
out_file3=out_file2.replace('_5dist_','_')
cmd='grep TSS ' + out_file2 + '  > ' + out_file3
os.system(cmd)
print('Export TSS overlapping with DMRs: ')
#print('\t',out_file3)

#9. to add new column name
out_df=pd.read_csv(out_file3,sep='\t',header=None)
out_columns=['gene_name','gene_type','new_mr_sites']
out_df.columns=out_columns
out_df.to_csv(out_file3, sep='\t',index=False)

#10. to find 5dist regions overlapping with enhancers
out_file3=in_file.replace('_tss_','_')
cmd='grep 5dist ' + in_file + '  > ' + out_file3
print('Export 5distance region overlapping with enhancers')
os.system(cmd)

out_file4=out_file3.replace('_rat','_rat_overlap_enhancer')
cmd2='bedtools intersect -a '+ out_file3 + \
       ' -b ' + os.path.join(in_genome_folder , in_genome_file) + ' -wa -wb ' + \
       '  > ' + out_file4
os.system(cmd2)
print('Export 5dist overlapping with enhancers: ')
#print(out_file4)

#11. reformato file  out_file4 for output
in_file2=out_file4
in_df2=pd.read_csv(in_file2,sep='\t',header=None)
in_df2['mr_site']=in_df2[3].apply(lambda x: ':'.join(x.split(':')[0:2]))
in_df2['gene_type']=in_df2[5].apply(lambda x: x.split('||')[1])
out_df2=in_df2.groupby([6,'gene_type'])['mr_site'].apply('~'.join).reset_index().copy()
out_df2.columns=out_columns
out_file5=in_file2.replace('dmr_','uqGeneDmr_')
out_df2.to_csv(out_file5,sep='\t',index=False)
#here new_mr_sites was changed manually to block_id for doing chromSegment_test4blocks 
print('Convert file format to unique gene linked DMRs ..')
print(out_file5,'\n')

print('Remove temporary files: ')
print('\n', out_file2,'\n', out_file3,'\n', out_file4,'\n')
os.system('rm -f ' + out_file2)
os.system('rm -f ' + out_file3)
os.system('rm -f ' + out_file4)
print('\n')






