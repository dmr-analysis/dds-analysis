#this script is used to find enhancer target genes based on dds output gene list from final GO pathway analysis
import pandas as pd
import glob
import os
import numpy as np
import argparse
from .sort_filter_bed_file import sort_and_filter
#exec(open("find_enhancer_target_genes.py").read())

def my_parser(parser):
  required= parser.add_argument_group('Required')
  required.add_argument('-inEnhancerF','--in_enhancer_file_folder', required=True, help='file folder name for enhancer files where assuemes all files end with *.bed  ')
  required.add_argument('-inDDSFile','--in_DDS_file',required=True, help='file name of DDS final exported ranked gene list')
  required.add_argument('-inGene','--in_selected_gene',required=True, help='gene name in DDS export file that will need be searched for overlapping with enhancers')
  required.add_argument('-outFolder','--out_folder', required=True, help='Output file folder name')
  option = parser.add_argument_group('Optional, with default values')
  option.add_argument('-isAll','--is_all_blocks',help='Whether use all blocks or enhance only blocks for the enhancer-gene search, default is enhancer only blocks',
                      action="store_true", default=False)
  return parser

def find_target_genes(in_enhancer_file_folder, in_dds_file, select_gene, out_folder, is_all_blocks):
 '''find overlapping between mutation blocks assocated to selected genes and all enhancers
    . And also find the corresponding target genes for the enhancers '''
 #in_dds_file=args.in_DDS_file
 #select_gene=args.in_selected_gene
 #out_folder=args.out_folder
 in_enhancer_files=glob.glob(os.path.join(in_enhancer_file_folder,'*.bed'))

 #find compressed bed files
 if len(in_enhancer_files)==0:
     in_enhancer_files0=glob.glob(os.path.join(in_enhancer_file_folder,'*.bed.gz'))
     for efi in in_enhancer_files0:
         in_enhancer_files.append(efi.replace('.bed.gz','.bed'))
         os.system('gunzip ' + efi)

 if len(in_enhancer_files)==0:
     print('No file find in :', in_enhancer_file_folder)
     print('Please check enhancer file location, program stop!')
     exit()

 #select mutation blocks that are linked to a gene for searching their enhancers and target genes
 dds_df=pd.read_csv(in_dds_file,sep='\t')
 selected_df=dds_df[dds_df.gene_name==select_gene]

 enhancers=np.array(selected_df.enhancers.str.split('~').tolist()[0],dtype=str)
 is_enhancer=enhancers=='enhancer'
 block_ids=np.array(selected_df.block_id.str.split('~').tolist()[0],dtype=str)
 mr_sites=np.array(selected_df.new_mr_sites.str.split('~').tolist()[0],dtype=str)
 gene_type=np.array(selected_df.gene_type.str.split('~').tolist()[0],dtype=str)
 patients=np.array(selected_df.patients.str.split('~').tolist()[0],dtype=str)
 if not is_all_blocks:
   #find blocks overlap with enhancers
   print('Only consider blocks overlapping to 5dist enhancers ')
   enhancer_block_ids=block_ids[is_enhancer]
   enhancer_mr_sites=mr_sites[is_enhancer]
   enhancer_gene_type=gene_type[is_enhancer]
   enhancer_patients=patients[is_enhancer]
   enhancers2=enhancers[is_enhancer]
 else:
   print('Consider all available blocks')
   #use all available blocks
   enhancer_block_ids=block_ids
   enhancer_mr_sites=mr_sites
   enhancer_gene_type=gene_type
   enhancer_patients=patients
   enhancers2=enhancers

 #make output dataframe
 out_df=pd.DataFrame(columns=['block_ids','mr_sites','gene_type','patients','enhancers'])
 out_df['block_ids']=enhancer_block_ids
 out_df['mr_sites']=enhancer_mr_sites
 out_df['gene_type']=enhancer_gene_type
 out_df['patients']=enhancer_patients
 out_df['enhancers']=enhancers2
 #out_df.enhancers=out_df.enhancers.str.strip()

 out_df_chr_pos=pd.DataFrame(data=out_df.block_ids.str.split('_').to_list()).iloc[:,2:].copy()
 out_df_chr_pos.columns=['chrs','start_pos','end_pos']

 out_df2=pd.concat([out_df_chr_pos,out_df],axis=1)
 if out_df2.chrs[0]=='23':
  out_df2.chrs.replace('23','X')
 if out_df2.chrs[0]=='24':
  out_df2.chrs.replace('24','Y')
 if 'chr' not in out_df2.chrs:
  out_df2.chrs='chr'+out_df2.chrs

 if not os.path.exists(out_folder):
   os.makedirs(out_folder)
   print('Create folder: ', out_folder)
 out_folder2=os.path.join(out_folder,select_gene)
 if not os.path.exists(out_folder2):
   os.makedirs(out_folder2)
   print('Create folder: ', out_folder2)

 out_file=os.path.join(out_folder2,select_gene+'_blocks.bed')
 out_df2.to_csv(out_file,sep='\t',index=False, header=None)
 print('Export: ', out_file)

 #sort exported bed files
 #cmd='python sort_filter_bed_file.py -inFile '+ out_file + ' -noLength'
 #os.system(cmd)
 isLength=False
 isFilter=True
 ishuman=True
 sort_and_filter(out_file,isLength,isFilter, ishuman) 
 os.system('rm -f '+ out_file)

 #search for overlapping in enhancer target files
 min_overlap=1e-9
 in_folder=os.path.split(out_file)[0]
 in_gene_block_file=glob.glob(os.path.join(in_folder,select_gene+'*.bed'))
 record_enhancers=[]
 for gi in in_gene_block_file:
  for fi in in_enhancer_files:
   #fi = in_enhancer_files[0]
   out=os.path.join(in_folder,os.path.basename(fi).replace('.bed','')+'_'+os.path.basename(gi).replace('.bed','')+'.bed')
   out=out.replace('_sorted','')
   #print(out)
   cmd='bedtools intersect -a ' + gi + ' -b ' + fi + \
                  ' -wa -wb -f ' + str(min_overlap) + ' > ' + out
   os.system(cmd)

   outsize=os.path.getsize(out)
   if outsize==0:
     #print(out, ' is empty and remove !')
     cmd='rm -f ' + out
     os.system(cmd)
   else:
     print(out)
     record_enhancers.append(out)

 #compress enhancer files
 for ei in in_enhancer_files:
     os.system('gzip '+ ei)

 return record_enhancers

def run(args):
 in_dds_file=args.in_DDS_file
 select_gene=args.in_selected_gene
 out_folder=args.out_folder
 in_enhancer_folder=args.in_enhancer_file_folder
 is_all_blocks=args.is_all_blocks
 find_target_genes(in_enhancer_folder,in_dds_file, select_gene, out_folder,is_all_blocks)

if __name__=='__main__':
 #for enhancer target genes
 #in_dds_file='out_blocks_dmr_single3/out_dds/out_gene_DMR_or_DEG_in_GO_KEGG_BIOC_pathways_p_ls0.05.csv'
 #in_enhancer_files=glob.glob('in_data/human/in_enhancer/hg19_enhancer2gene_bed/*.bed')
 #select_gene='CTSO'
 #out_folder='out_enhancers'

 args=my_parser(argparse.ArgumentParser('python find_enhancer_target_genes.py ')).parse_args()
 run(args)


