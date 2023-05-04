#this script is used to map mutation blocks to genome regions
#exec(open("map_block2genome.py").read())
import os
import pandas as pd
import argparse
#from bpb3.script.script_high.other.common import check_folder
from .script_high.functions_from_bpb3 import check_folder, is_non_zero_file

def my_parser(parser):
  required= parser.add_argument_group('Required')
  required.add_argument('-inSBfile','--in_sortedBlock_file', required=True, help='Bpb3 exported mutation block summary position file in BED format and sorted by chromosome and position')
  required.add_argument('-inGRfile','--in_genomeRegion_file', required=True, help='a list of genome region files which cotains TSS,TES, gene, position et al ')
  required.add_argument('-inRGfile','--in_referenceGenome_file',required=True, help='sorted reference genome refFlat file in BED format')

  optional= parser.add_argument_group('Optional, has default values')
  optional.add_argument('-outFold','--out_file_folder', default='out_block_single',metavar='', type=str, help='output file folder name')
  optional.add_argument('-is_change','--change_name', action="store_true", help='whether to change the 4th column of name by reducing its length, columns 1,2,3 are chr, start_pos, and end_pos,\
                      but column 4 is the id information, default =False, does not change id information . If use this option, then program will split the name by : and only take the first two elements for exporting a new id name')
  optional.add_argument('-in_cutoff','--in_dmr_minimum_cutoff', default=None, type=float, help='Minimum cutoff values for select blocks or mrs to export, default =None that means no minimum cutoff for exporting data or there is not 9th column\
                                              (p-value or other values can be used to filter by mimimum cutoff) in dataframe')
  optional.add_argument('-isMorB', '--is_MR_or_Blocks', default=0, type=int, help="Is input bed position files are MR (methylation regions chr#:mr#) or mutation blocks block_#_chr_start_end, 0 for MR, 1 for Blocks, default=0 is MR")
  return parser

def is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def count_dmrs_not_mapped2genome(in_sorted_block_file,record_out_files,dmr_min_cutoff,is_MR_or_Blocks):
  #count how many DMRs/blocks are not mapped to annotated geneomic regions
  #coumt MR or DMR/blocks in genomic files
  if dmr_min_cutoff is not None:
     all_indata_df_columns=['chrs','start_pos','end_pos','genome_info','mr_chrs','mr_start_pos','mr_end_pos','mr_sites','mr_logReg_proba']
  else:
     all_indata_df_columns=['chrs','start_pos','end_pos','genome_info','mr_chrs','mr_start_pos','mr_end_pos','mr_sites']
  in_dmrs_df_columns=all_indata_df_columns[4:]
  #print(in_dmrs_df_columns)
  all_data_df=[]

  record_matched_block_in_file={}
  for fil in record_out_files:
    if is_non_zero_file(fil) and os.stat(fil).st_size>0:
      tmp_data_df=pd.read_csv(fil,header=None, sep='\t')
      all_data_df.append(tmp_data_df.copy())

      #for each record_out_file count how many dmr or blocks are matched
      tmp_data_df.columns=all_indata_df_columns
      uq_tmp_data_df=tmp_data_df.drop_duplicates(subset=['mr_sites'])
      record_matched_block_in_file[os.path.basename(fil)]=uq_tmp_data_df.copy()

  #combine all dataFrame
  all_indata_df=pd.concat(all_data_df)
  all_indata_df.columns=all_indata_df_columns
  uq_indata_df=all_indata_df.drop_duplicates(subset=['mr_sites'])
  total_uq_mrs=uq_indata_df
  if dmr_min_cutoff is not None:
    min_cutoff=dmr_min_cutoff
    total_uq_dmrs=uq_indata_df[uq_indata_df['mr_logReg_proba']>=min_cutoff]
  else:
    total_uq_dmrs=uq_indata_df

  in_dmrs_df=pd.read_csv(in_sorted_block_file,header=None,sep='\t')
  in_dmrs_df.columns=in_dmrs_df_columns
  total_in_mrs=in_dmrs_df
  if dmr_min_cutoff is not None:
     total_in_dmrs=in_dmrs_df[in_dmrs_df['mr_logReg_proba']>=min_cutoff]
  else:
     total_in_dmrs=total_in_mrs

  #print('Number of input Blocks do not find mapped genome information')
  #print(total_uq_mrs.shape[0]-total_in_mrs.shape[0])
  #print(total_uq_dmrs.shape[0]-total_in_dmrs.shape[0])

  #print('Perentage of input Blocks mapped to genome info')
  #print(total_uq_mrs.shape[0]/total_in_mrs.shape[0])
  #print(total_uq_dmrs.shape[0]/total_in_dmrs.shape[0])

  diff_in_mr=set(total_in_mrs.mr_sites.to_list())- set(total_uq_mrs.mr_sites.to_list()) 
  diff_in_dmr=set(total_in_dmrs.mr_sites.to_list())- set(total_uq_dmrs.mr_sites.to_list())
  percentage_mapped_block=total_uq_mrs.shape[0]/total_in_mrs.shape[0]

  #check the distribution of unmapped MRs in chromes
  chrs=[]
  for i in range(1,25):
    if i<23:
        chrs.append('chr'+str(i))
    elif i==23:
        chrs.append('chrX')
    elif i==24:
        chrs.append('chrY')

  diff_in_dmr_df=pd.DataFrame(data=list(diff_in_dmr),columns=['dmr_sites'])
  dict_chr={}
  total_not_mapped_block=0
  for ii in chrs:
    if is_MR_or_Blocks ==0 :
        #this is for DMR block
        #print('Mehytlation Regions ')
        dict_chr[ii]=diff_in_dmr_df[diff_in_dmr_df.dmr_sites.str.contains(ii+':')]
    else:
        #this is for mutation block
        #print('Mutation Bloks')
        #print(diff_in_dmr_df.dmr_sites)
        dict_chr[ii]=diff_in_dmr_df[diff_in_dmr_df.dmr_sites.apply(lambda x: x.split('_')[2]==ii.replace('chr',''))]
    total_not_mapped_block += dict_chr[ii].shape[0]

  #remove dictionary key with empty value/dataframe
  not_mapped_block_dict_chr={k: v for k,v in dict_chr.items() if not v.empty}
  return total_not_mapped_block, percentage_mapped_block, not_mapped_block_dict_chr, record_matched_block_in_file, diff_in_dmr_df

def bed_region_match4_list_genomic_regions(region_files, methylation_file, reference_file, out_folder, min_overlap):
  #this function is used to match a bed format region to a list of predifined genomic regions
  #for each region to find its Blocks
  record_out_files=[]
  for fil in region_files:
   if is_non_zero_file(fil) and os.stat(fil).st_size>0:
     region_file=fil
     region_name=os.path.basename(fil).split('_')[0].lower()
     print(region_name)
     out_methylation_name = out_folder + '/' + os.path.basename(methylation_file)[:-4]
     out_region_name = os.path.basename(region_file)[:-4]

     #print(methylation_file)
     #print(region_file)
     dist5_methylation_file = out_methylation_name + '_' + 'noGenes.bed'
     print(dist5_methylation_file)

     #print(out_methylation_name,out_region_name)

     if region_name == '5dist':
        # For 5distance we first remove genes(TSS, geneBody and TES) from the two methylation files
        os.system('bedtools intersect -a ' + methylation_file + ' -b ' + reference_file + \
                  ' -v > ' + dist5_methylation_file)

        out = dist5_methylation_file[:-4] + '_' + out_region_name + '_overlap' + str(min_overlap) + '.bed'

        os.system('bedtools intersect -a ' + region_file + ' -b ' + \
                  dist5_methylation_file + ' -wa -wb -f ' + str(min_overlap) + ' > ' + out)

        os.system('rm ' + dist5_methylation_file)  # removes temporary file
     else:
        out = out_methylation_name + '_' + out_region_name + '_overlap' + str(min_overlap) + '.bed'
        os.system('bedtools intersect -a ' + region_file + ' -b ' + methylation_file + \
                  ' -wa -wb -f ' + str(min_overlap) + ' > ' + out)
     print(out)
     record_out_files.append(out)
   else:
     print('Not find or empty file of: ', fil)
     print('Error in Program map_block2genome stop! ')
     exit()
  return record_out_files

def run(args):
  in_sorted_block_file=args.in_sortedBlock_file
  in_region_files=args.in_genomeRegion_file
  reference_file=args.in_referenceGenome_file
  out_folder=args.out_file_folder
  check_folder(out_folder)
  
  print("Start map blocks to predefined genomic regions ")
  print("In block postion file: ", in_sorted_block_file)
  print("In predefined regions file: ", in_region_files)
  print("In reference genoome file: ", reference_file)
  print("Output folder: ", out_folder)

  # input data
  #block of minimum 2 patients
  #in_sorted_block_file='bp2_mussd_blocks_alireza/blocks_summary_sorted.bed'

  #block of minimum 1 patients
  #in_sorted_block_file='bp2_mussd_blocks_single/blocks_summary_block_position.bed'
  #in_region_files='genome_data/list_region_files.txt'

  region_files=pd.read_csv(in_region_files,header=None)
  region_files=region_files.loc[:,0].to_list()
  #print(region_files)
  methylation_file=in_sorted_block_file
  dmr_min_cutoff=args.in_dmr_minimum_cutoff
  #change the fourth column name of input bed file 
  if args.change_name:
     print('Change the name of fourth column in bed file ', methylation_file)
     bed_df=pd.read_csv(methylation_file,sep='\t',header=None)
     bed_df[3]=bed_df[3].apply(lambda x: ':'.join(x.split(':')[0:2]) )
     out_file='_'.join([methylation_file,'shortName'])
     print("\n")
     print('Export name changed bed file in ', out_file)
     if dmr_min_cutoff == None:
       print('Only read the first four columns ! ')
       out_bed_df=bed_df.loc[:,0:3].copy()
     else:
       out_bed_df=bed_df.copy()

     out_bed_df.to_csv(out_file,sep='\t',index=False, header=None)
     methylation_file=out_file

  #reference_file='genome_data/data/hg19.refFlat_clean_sorted.bed'
  min_overlap=1E-9

  #output data
  #out_folder='out_blocks_single/'
  
  record_out_files=bed_region_match4_list_genomic_regions(region_files, methylation_file,reference_file, out_folder, min_overlap)

  #count how many DMRs are not mapped to annotated geneomic regions
  #coumt MR or DMR in genomic files
  #dmr_min_cutoff=0.8
  #dmr_min_cutoff=None
  #dmr_min_cutoff=args.in_dmr_minimum_cutoff 

  #remove an elemenet before checking
  #record_out_files.pop()
  record_out_files2=[x for x in record_out_files if 'intergenic' not in x  ]
  
  total_not_mapped_block,percentage_mapped_block, not_mapped_block_dict_chr, record_matched_block_in_file , diff_in_dmr_df =count_dmrs_not_mapped2genome(methylation_file,record_out_files2,dmr_min_cutoff,args.is_MR_or_Blocks)

  #print recorded matached genomic regions
  #for kk in record_matched_block_in_file.keys():
      #print(kk,record_matched_block_in_file[kk].shape)
  #print(not_mapped_block_dict_chr)
  
  if True:
    #export mapped statistic information 
    out_df=pd.DataFrame(columns=['ID','information'])
    list2regions=list(record_matched_block_in_file.keys())
    out_df_ID=['total_not_mapped_blocks','percentage_of_mapped_blocks']+list2regions+['not_mapped_blocks']
    out_df_info=[total_not_mapped_block, percentage_mapped_block]+ [record_matched_block_in_file[x].shape[0] for x in list2regions ]+[','.join(diff_in_dmr_df.dmr_sites.to_list())]
    out_df.ID=out_df_ID
    out_df.information=out_df_info
    out_file=os.path.join(out_folder,os.path.basename(methylation_file).replace('.bed','_mapped_genome_information.tsv'))
    print("\n")
    print("Export mutation blocks mapped to predefined genomic regions and its summary information at: ")
    print(out_file)
    out_df.to_csv(out_file,sep='\t',index=False)


if __name__=='__main__':
  args=my_parser(argparse.ArgumentParser('python map_block2genome.py ')).parse_args()
  run(args)


