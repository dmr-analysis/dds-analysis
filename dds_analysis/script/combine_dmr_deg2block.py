#this script is used to combine DMR DEG results to blocks
import pandas as pd
import glob
import os
import argparse
#from bpb3.script.script_high.other.common import check_folder
from .script_high.functions_from_bpb3 import check_folder

#exec(open('combine_dmr_deg2block.py').read())
def my_parser(parser):
  required = parser.add_argument_group('Required')
  required.add_argument('-in_SBfile','--in_sortedBlock_patient_file', required=True,
                  help='Bpb3 exported mutation block summary position file in BED format and sorted by chromosome and position, where patient IDs are added by dds_analysis find_block_patientID ')
  required.add_argument('-in_DMRfile','--in_dmr_file', required=True, help='bpb3 block summary file with DMR information that exported by dds_analysis map_block2dmr')
  required.add_argument('-in_DEGfiles','--in_deg_folder_and_file_suffix',required=True, 
            help='Path of a file folder that contains results exported by dds_analysis find_geneExp4block, where DEG information is added to blocks (e.g., file_path/*.tsv) ')

  optional= parser.add_argument_group('Optional, has default values')
  optional.add_argument('-outFold','--out_file_folder', default='out_DmrDeg2block',metavar='', type=str, help='Path of output file folder ')

  return parser

def run(args):
 in_block_file=args.in_sortedBlock_patient_file
 in_dmr_file=args.in_dmr_file
 in_deg_folder_files=args.in_deg_folder_and_file_suffix
 out_folder=args.out_file_folder
 check_folder(out_folder)

 print("Start to combine DMR and DEG information in blocks ")
 print("In block position file with patieint ID: ", in_block_file)
 print("In DMR file with block information: ", in_dmr_file)
 print("In DEG file with blocks assigned to predefined genomic regions: ", in_deg_folder_files)
 print("Output folder: ", out_folder)

 #input data
 #SNP
 #minimum 2 patients
 #in_block_file='bp2_mussd_blocks_alireza/blocks_summary.tsv'
 #minimum 1 patient
 #in_block_file='bp2_mussd_blocks_single/blocks_summary.tsv'
 #in_block_file='bp2_mussd_blocks_single/blocks_summary_and_patientID.tsv'
 in_blocks_df=pd.read_csv(in_block_file,sep='\t',low_memory=False)

 #DMR
 #in_dmr_file='out_blocks_dmr/blocks_summary_sorted_500flank_0.7Proba_176blocks_73blocks2mr_30blocks2dmr.tsv'
 #in_dmr_file='out_blocks_dmr_single/blocks_summary_block_position_500flank_0.7Proba_66868blocks_28049blocks2mr_9478blocks2dmr.tsv'
 #in_dmr_file='out_blocks_dmr_single/blocks_summary_block_position_0flank_0.7Proba_66868blocks_13143blocks2mr_4604blocks2dmr.tsv'

 in_dmr_df= pd.read_csv(in_dmr_file,sep='\t')
 in_dmr_df['new_mr_sites']=in_dmr_df.mr_sites.apply(lambda x: ':'.join(x.split(':')[0:2]))

 #DEG
 #in_deg_files=glob.glob('out_blocks/out_expression/*.tsv')
 #in_deg_files=glob.glob('out_blocks_single/out_expression/*.tsv')
 in_deg_files=glob.glob(in_deg_folder_files)
 record_deg={}
 for fi in in_deg_files:
    tmp_df=pd.read_csv(fi,sep='\t')
    tmp_name=os.path.basename(fi).replace('.tsv','')
    #print(tmp_name)
    record_deg[tmp_name]=tmp_df.copy()

 #add DMR to block
 dmr2block_df=pd.merge(in_blocks_df,in_dmr_df,on='block_id',how='left')

 grouped = dmr2block_df.groupby(['block_id'])
 result = grouped.agg({'new_mr_sites':lambda x: ', '.join(tuple(x.astype(str).tolist())),
                      'mr_logReg_proba' : lambda x: ', '.join(tuple(x.astype(str).tolist())) }
                     )
 result=result.reset_index()
 tmp_result=result.copy()
 #add DEG to block
 for kk in record_deg.keys():
   tmp_df=record_deg[kk].copy()
   tmp_df=tmp_df.rename(columns={'differential_expression_T_test_p_value': 'deg_p_value'} )
   tmp_deg2block_df=pd.merge(tmp_result,tmp_df,on='block_id',how='left')
   tmp_grouped=tmp_deg2block_df.groupby(['block_id'])
   tmp_result2= tmp_grouped.agg({ 'gene': lambda x: ' ,'.join(tuple(list((x.astype(str).tolist())))),
                            'deg_p_value': lambda x: ','.join(tuple(list((x.astype(str).tolist())))) }
                             )
   tmp_result2=tmp_result2.reset_index()
   tmp_result2=tmp_result2.rename(columns={'gene': kk.split('_')[-1],'deg_p_value': 'deg_p_value2' + kk.split('_')[-1]})
   #merge deg with dmr
   #here need to remove gene with p-value == nan??
   result=pd.merge(result,tmp_result2,on='block_id',how='left')

 #export results
 out_result=pd.merge(in_blocks_df,result,on='block_id',how='left')
 out_file=os.path.basename(in_dmr_file).replace('.tsv','')+'_deg_info.tsv'
 out_file=os.path.join(out_folder, out_file)
 print("Export results at: ")
 print(out_file)
 out_result.to_csv(out_file,sep='\t',index=False)


if __name__== '__main__':
 args=my_parser(argparse.ArgumentParser('python combine_dmr_deg2block.py')).parse_args()
 run(args)

