#this script is used to combine genome and enhancer annotation for block and find its relevant differentially expressed genes
#exec(open('find_geneExp4block.py').read())
import pandas as pd
import glob
import os
import argparse
#from bpb3.script.script_high.other.common import check_folder, is_non_zero_file
from .script_high.functions_from_bpb3 import check_folder, is_non_zero_file


def my_parser(parser):
  required= parser.add_argument_group('Required')
  required.add_argument('-in_BGfolder','--in_blocks_genome_folder',required=True, help='Path of blocks mapped to genomic region files, which are exported by map_block2genome ')
  required.add_argument('-in_SBfile','--in_sortedBlock_file',required=True, help='Bpb3 exported mutation block summary position file in BED format and sorted by chromosome and position')
  required.add_argument('-in_DGfile','--in_deg_file', required=True,help='A differentially expressed gene file in tab delimiate format, which is an export file from bpb3 differential_expression ')
  required.add_argument('-out_folder','--out_file_folder',required=True, help='Path of file folder to export results ')
  
  optional = parser.add_argument_group('Optional, has default values')
  optional.add_argument('-in_FTlist','--in_feature_list', default='TSS,gene,TES,5dist,enhancers',metavar='',type=str,
                         help='A comma separated string represents features of genomic regions that will be considered at here, default = TSS,gene,TES,5dist,enhancers, where each feature is separted by commoa')
  optional.add_argument('-in_cutoff','--in_dmr_minimum_cutoff', default=None, type=float, help='Minimum cutoff values for select blocks or mrs to export, default =None that means no minimum cutoff for exporting'\
					' data or there is not 9th column (p-value or other values can be used to filter by mimimum cutoff) in dataframe')

  return parser 

def run(args):
 in_blocks_genome_folder=args.in_blocks_genome_folder
 feature_list=args.in_feature_list.split(',')
 in_sorted_block_file=args.in_sortedBlock_file
 in_deg_file=args.in_deg_file
 out_folder=args.out_file_folder
 check_folder(out_folder)
 in_cutoff=args.in_dmr_minimum_cutoff

 print("Start to find DEG for mapped genomic regions ")
 print("In block mapped genomic regions: ", in_blocks_genome_folder)
 print("In block position file: ", in_sorted_block_file)
 print("In DEG file: ", in_deg_file)
 print("Feature list: ", feature_list)
 print("Output folder: ", out_folder)

 #in_blocks_genome_folder='out_blocks/'
 #in_blocks_genome_folder='out_blocks_single/'
 #feature_list=['TSS','gene','TES','5dist','enhancers']

 in_blocks_genome_files=glob.glob(os.path.join(in_blocks_genome_folder,'*.bed'))
 true_blocks_genome_files=[]
 for x in in_blocks_genome_files:
    is_feature=0
    for fe in feature_list:
        if '_' + fe + '_' in x and 'intergenic' not in x:
            is_feature +=1
    if is_feature>0:
       true_blocks_genome_files.append(x)
 #print(true_blocks_genome_files)

 #read_block_positon
 #block for minimum 2 patients
 #in_sorted_block_file='bp2_mussd_blocks_alireza/blocks_summary_sorted.bed'

 #block for minimum 1 pateint
 #in_sorted_block_file='bp2_mussd_blocks_single/blocks_summary_block_position.bed'

 #in_blocks_df=pd.read_csv(in_sorted_block_file,sep='\t',header=None)
 #in_blocks_df.columns=['chrs','start_pos','end_pos','block_id']

 #read differential expression genes
 #in_deg_file='bp2_diffexp_result_alireza/differentially_expressed_genes.txt'
 in_deg_df=pd.read_csv(in_deg_file,sep='\t')
 in_columns=list(in_deg_df.columns)
 new_columns=[]
 for ic in in_columns:
  if '#' in ic:
      new_columns.append(ic.replace('#',''))
  else:
      new_columns.append(os.path.basename(ic))
 in_deg_df.columns=new_columns

 #read all annoated blocks
 record_features={}
 for fi in true_blocks_genome_files:
   if  is_non_zero_file(fi) and os.stat(fi).st_size>0:  
     tmp_df=pd.read_csv(fi,sep='\t',header=None)
     if in_cutoff==None:
        tmp_df.columns=['chrs','start_pos','end_pos','feature','bk_chrs','bk_start_pos','bk_end_pos','block_id']
     else:
        tmp_df.columns=['chrs','start_pos','end_pos','feature','bk_chrs','bk_start_pos','bk_end_pos','block_id','block_value']
     #record df
     tmp_fi=set(fi.split('_'))
     feature=list(tmp_fi.intersection(set(feature_list)))
   
     if 'enhancer' in feature[0]:
       tmp_df['new_feature']=tmp_df[['chrs','start_pos','end_pos']].apply(lambda x: ':'.join(x.astype(str)),axis=1)+':enhancer'
       tmp_df['gene']=['enhancer']*tmp_df.shape[0]
     else:
       tmp_df['new_feature']=tmp_df.feature.str.split('\|\|')
       record_nf=[]
       record_gene=[]
       for nf in tmp_df.new_feature:
          record_nf.append('||'.join(nf[0:2]))
          record_gene.append(nf[2].split(':')[0])
       tmp_df['new_feature']=record_nf
       tmp_df['gene']=record_gene 
       #input('click')
     if len(feature)>0:
       print(feature)
       #find differential gene expression
       merged_tmp_df=pd.merge(tmp_df,in_deg_df,on='gene',how='left')
       record_features[feature[0]]=merged_tmp_df[['bk_chrs','bk_start_pos','bk_end_pos','block_id','chrs','start_pos','end_pos','new_feature','gene','differential_expression_T_test_p_value']].copy()
     else:
       print(fi, ' Genome Feature not find!')
   else:
       print(fi, ' not exist or empty, skip reat it !')

 #export data
 #out_folder=in_blocks_genome_folder+'/out_expression/'
 out_file_prefix=os.path.basename(in_sorted_block_file).replace('.bed','')
 print("\n")
 print("Export files at: ")
 for kk in record_features.keys():
     out_file=os.path.join(out_folder, out_file_prefix + '_'+kk + '.tsv')
     tmp_df=record_features[kk].copy()
     tmp_df.to_csv(out_file,sep='\t',index=False)
     print(out_file)


if __name__=='__main__':
 args= my_parser(argparse.ArgumentParser('python find_geneExp4block.py ')).parse_args()
 run(args)



