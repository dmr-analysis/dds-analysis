#this script file is used to map block to predicted chromSegment data
#import map_block2genome
from .map_block2genome import count_dmrs_not_mapped2genome
import glob
import os
import pandas as pd
import argparse
import subprocess
#from bpb3.script.script_high.other.common import check_folder
from .script_high.functions_from_bpb3 import check_folder


#exec(open('map_block2chromSegment.py').read())
#map DMR to chromatin segmant files

def my_parser(parser):
    required= parser.add_argument_group('Required')
    required.add_argument('-in_BSpos','--in_sortedBlock_file', required=True, help='Bpb3 exported mutation block summary position file in BED format and sorted by chromosome and position')
    required.add_argument('-out_folder','--out_file_folder', required=True, help='Export file folder')
    required.add_argument('-in_chrSeg','--in_chromSegment_folder',required=True, help='File folder path of chromatin segmement files')

    optional = parser.add_argument_group('Optional , has default values')
    optional.add_argument('-in_prefStr','--in_prefix_string4chromSegment_file', default='combined*merged.bed.gz',metavar='',  type=str, help='Prefix string for chromatin segment file names , default= combined*merged.bed.gz')
    optional.add_argument('-is_change','--change_name', action="store_true", help='whether to change the 4th column of name by reducing its length, columns 1,2,3 are chr, start_pos, and end_pos,\
                      but column 4 is the id information, default =False, does not change id information . If use this option, then program will split the name by : and only take the first two elements for exporting a new id name')
    optional.add_argument('-in_cutoff','--in_dmr_minimum_cutoff', default=None, type=float, help='Minimum cutoff values for select blocks or mrs to export, default =None that means no minimum cutoff for exporting data or there is not 9th column\
                                              (p-value or other values can be used to filter by mimimum cutoff) in dataframe')
    optional.add_argument('-isMorB', '--is_MR_or_Blocks', default=0, type=int, help="Is input bed position files are MR (methylation regions chr#:mr#) or mutation blocks block_#_chr_start_end, 0 for MR, 1 for Blocks, default=0 is MR")
    return parser

def submit_job_in_subprocess(record_cmds):
    processes=[]
    for cmd in record_cmds:
        processes.append(subprocess.Popen(cmd,shell=True))
    for p in processes:
        p.communicate()

def run(args):
 in_sorted_block_file=args.in_sortedBlock_file
 out_folder= args.out_file_folder
 chromSegment_file_string=os.path.join(args.in_chromSegment_folder, args.in_prefix_string4chromSegment_file)
 check_folder(out_folder)

 print("Start map blocks to chromatin segmentations ")
 print("In block position file: ", in_sorted_block_file)
 print("In chromatin segmentation file: ", chromSegment_file_string)
 print("Output folder: ", out_folder)
 
 #block for minimum 2 patients
 #in_sorted_block_file='bp2_mussd_blocks_alireza/blocks_summary_sorted.bed'

 #block for minimum 1 patients
 #in_sorted_block_file='bp2_mussd_blocks_single/blocks_summary_block_position.bed'
 methylation_file=in_sorted_block_file
 dmr_min_cutoff=args.in_dmr_minimum_cutoff
 #change the fourth column name of input bed file 
 if args.change_name:
    print('Change the name of fourth column in bed file ', methylation_file)
    bed_df=pd.read_csv(methylation_file,sep='\t',header=None)
    bed_df[3]=bed_df[3].apply(lambda x: ':'.join(x.split(':')[0:2]) )
    out_file='_'.join([methylation_file,'shortName'])
    print('Export name changed bed file in ', out_file)
    if dmr_min_cutoff == None:
       print('Only read the first four columns ! ')
       out_bed_df=bed_df.loc[:,0:3].copy()
    else:
       out_bed_df=bed_df.copy()
    
    out_bed_df.to_csv(out_file,sep='\t',index=False, header=None)
    methylation_file=out_file
 
 #out_folder='out_blocks/chromatin_segment/'
 #out_folder='out_blocks_single/chromatin_segment/'
 min_overlap=1E-9
 #merged_out_files=glob.glob('genome_data/chromatin_segment/combined*merged.bed')
 merged_out_files=glob.glob(chromSegment_file_string)

 if len(merged_out_files)<1:
    print('Error: chromatin segmentation files are not found, please check file name: ', chromSegment_file_string)
    print('Program stop !\n\n') 
    exit()

 record_bed_files=[]
 processes=[]
 for fi in merged_out_files:
     if '.gz' in fi:
        cmd='gunzip -f ' + fi
        processes.append(subprocess.Popen(cmd, shell=True))
        record_bed_files.append(fi.replace('.gz',''))
     else:
        record_bed_files.append(fi)

 for P in processes:
     P.communicate()

 if len(record_bed_files)>0:
     merged_out_files=record_bed_files

 record_out_files=[]
 record_cmds=[]
 #print(merged_out_files)
 for fil in merged_out_files:
    region_file=fil
    out_methylation_name = out_folder + '/' + os.path.basename(methylation_file)[:-4]
    out_region_name = os.path.basename(fil)[:-4]
    out = out_methylation_name + '_' + out_region_name + '_overlap' + str(min_overlap) + '.bed'
    #os.system('bedtools intersect -a ' + region_file + ' -b ' + methylation_file + \
    #              ' -wa -wb -f ' + str(min_overlap) + ' > ' + out)
    print(out)
    cmd='bedtools intersect -a ' + region_file + ' -b ' + methylation_file + \
                 ' -wa -wb -f ' + str(min_overlap) + ' > ' + out
    record_cmds.append(cmd)
    record_out_files.append(out)
    #print(region_file, methylation_file, out, cmd)

 submit_job_in_subprocess(record_cmds)

 #count how many DMRs are not mapped to annotated geneomic regions
 #coumt MR or DMR in genomic files
 #dmr_min_cutoff=0.8
 #total, dict_chr=map_dmr2genome.count_dmrs_not_mapped2genome(in_sorted_dmr_file,record_out_files,dmr_min_cutoff)
 #dmr_min_cutoff=args.in_dmr_minimum_cutoff
 print('Minimum cutoff for 9th column or it is : ', str(dmr_min_cutoff) )
 #print(record_out_files)
 total_not_mapped_block,percentage_mapped_block, not_mapped_block_dict_chr, record_matched_block_in_file , diff_in_dmr_df =count_dmrs_not_mapped2genome(methylation_file,record_out_files,dmr_min_cutoff,args.is_MR_or_Blocks)

 #print recorded files
 #for kk in record_matched_block_in_file.keys():
 #     print(kk,record_matched_block_in_file[kk].shape)
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
    print("Export mutation blocks mapped chromatin segment files and its summary information at: ")
    print(out_file)
    out_df.to_csv(out_file,sep='\t',index=False)

 if len(record_bed_files)>0:
   for fi in record_bed_files:
        os.system('gzip -f ' + fi)

if __name__=='__main__':
 args= my_parser(argparse.ArgumentParser('python map_block2chromSegment.py ')).parse_args()
 run(args)


