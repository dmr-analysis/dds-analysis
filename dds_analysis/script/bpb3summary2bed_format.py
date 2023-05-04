#this script is used to convert mutation blocks from bpb3 to a bed format position file where block ID in the forth column
import pandas as pd
import os
import argparse

def my_parser(parser):
  required= parser.add_argument_group('Required')
  required.add_argument('-in_block','--in_block_summary_file',required=True, help='Input block summary file that is exported by bpb3 mussd analysis')
  optional= parser.add_argument_group('Optoinal, with default values')
  optional.add_argument('-out_path','--output_file_path', default='./', type=str, help='Output file path for exported new bed files, default= ./ ')
  return parser

def blocks_summary2bed_format(gene_block_df):
 '''convert bpb3 block summary file to a bed format postion file
    where the block ID in the forth column . Please not this demo works only
    for human genome, for other genome please modify the number chromosomes accordingly.
 '''
 #make a temp bed for for all block mutations
 all_blocks=[]
 for idx,row in gene_block_df.iterrows():
   tmp=[]
   tmp=row.block_id.strip().split(':')
   all_blocks += tmp
 uq_blocks=list(set(all_blocks))

 #find chr position for each block
 record_chr=[]
 record_st=[]
 record_ed=[]
 record_id=[]
 record_idx2chr=[]
 for bi in uq_blocks:
   tmp=bi.split('_')
   record_chr += ['chr'+tmp[2]]
   record_st += [tmp[3]]
   record_ed += [tmp[4]]
   record_id += [bi]
   if str.lower(tmp[2])=='x':
      tmp_chr=23
   elif str.lower(tmp[2])=='y':
      tmp_chr=24
   elif str.lower(tmp[2])=='m':
      tmp_chr=25
   else:
      tmp_chr=int(tmp[2])
   record_idx2chr += [tmp_chr]
 
 #make a dataframe for output data
 out_df=pd.DataFrame(columns=['chrs','start_pos','end_pos','block_id','idx2chr'])
 out_df.chrs=record_chr
 out_df.start_pos=record_st
 out_df.start_pos=out_df.start_pos.astype(int)
 out_df.end_pos=record_ed
 out_df.end_pos=out_df.end_pos.astype(int)
 out_df.block_id=record_id
 out_df.idx2chr=record_idx2chr

 #sort data by chromosome position
 out_df=out_df.sort_values(by=['idx2chr','start_pos'])
 out_df2=out_df[['chrs','start_pos','end_pos','block_id']]
 return out_df2.copy()

def run(args):
 #bpb3 exported mutation block summary file
 #in_gene_block_file='../../data/fl_12samples/out_data/demo3_fl_cohort_small/out/mussd_blocks/blocks_summary.tsv'
 in_gene_block_file=args.in_block_summary_file
 gene_block_df=pd.read_csv(in_gene_block_file,sep='\t')
 print('Read bpb3 block summary file: ', in_gene_block_file)

 out_df=blocks_summary2bed_format(gene_block_df)
 out_file=os.path.basename(in_gene_block_file.replace('.csv','_block_position.csv'))
 out_file=out_file.replace('.tsv','_block_position.bed')
 out_file2=os.path.join(args.output_file_path, out_file)
 print('Export bed format block summary file at: ', out_file2)
 out_df.to_csv(out_file2,sep='\t',header=None, index=False)
 return out_file2, out_df

if __name__== '__main__':
  args=my_parser(argparse.ArgumentParser('python bpb3summary2bed_format.py')).parse_args()
  run(args)

 #bpb3 exported mutation block summary file
 #in_gene_block_file='../../data/fl_12samples/out_data/demo3_fl_cohort_small/out/mussd_blocks/blocks_summary.tsv'
 #gene_block_df=pd.read_csv(in_gene_block_file,sep='\t')
 #print('Read bpb3 block summary file: ', in_gene_block_file)

 #out_df=blocks_summary2bed_format(gene_block_df)

# out_file=os.path.basename(in_gene_block_file.replace('.csv','_block_position.csv'))
# out_file=out_file.replace('.tsv','_block_position.bed')
# print('Export bed format block summary file at: ', out_file)
# out_df.to_csv(out_file,sep='\t',header=None, index=False)


