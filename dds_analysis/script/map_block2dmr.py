#this script is used to map block mutation to all available DMR/MR , where we can put flank regions on both sides of blocks before the mapping!
import pandas as pd
import os
import subprocess
import argparse 
#exec(open("map_block2dmr.py").read())

def my_parser(parser):
  required = parser.add_argument_group('Required')
  required.add_argument('-inBKfile','--in_sortedBlock_file', required=True, help='Bpb3 exported mutation block summary position file in BED format and sorted by chromosome and position')
  required.add_argument('-inDMRfile','--in_dmr_file', required=True, help='DMRs of all chromosomes in BED format that ranked exported by dmr_analysis dmr_combine_multChrs4rank ')

  optional = parser.add_argument_group('Optional , has default values')
  optional.add_argument('-inFlank','--in_flank_region2block', default=0, metavar='',type=int, help='Add N bp flank regions on both sides of the blocks, default =0 donot add flank regions')
  optional.add_argument('-outFolder','--out_file_folder',default='out_blocks_dmr_single/', metavar='',type=str, help='Path of output file folder, default= out_blocks_dmr_single/')
  optional.add_argument('-min_cutoff','--dmr_min_cutoff', default=0.6, metavar='', type=float, help='minimum cutoff value for selecting DMR regions from dmr_analysis exported file, default=0.6')
  return parser  

def count_blocks_not_mapped2dmr(in_sorted_block_file, record_out_file, dmr_min_cutoff):
  #find how many blocks are mapped to DMR and MR
  in_blocks_df=pd.read_csv(in_sorted_block_file,sep='\t',header=None)
  #print(in_blocks_df)
  in_blocks_df.columns=['chrs','start_pos','end_pos','block_id']  
  total_blocks=in_blocks_df.shape[0]

  in_blocks2mr_df=pd.read_csv(record_out_file,sep='\t',header=None)
  #print(in_blocks2mr_df)
  in_blocks2mr_df.columns=['chrs','start_pos','end_pos','block_id','mr_chrs','mr_start_pos','mr_end_pos','mr_sites','mr_logReg_proba']
  total_blocks2mr=in_blocks2mr_df.drop_duplicates(subset=['block_id']).shape[0]

  in_blocks2DMR_df=in_blocks2mr_df[in_blocks2mr_df.mr_logReg_proba>=dmr_min_cutoff].copy()
  total_blocks2DMR=in_blocks2DMR_df.drop_duplicates(subset=['block_id']).shape[0]

  print('minimum cutoff for DMR> %5.2f, total blocks %d , total blocks mapped to MR %d, total blocks mapped DMR %d' % (dmr_min_cutoff, total_blocks, total_blocks2mr,total_blocks2DMR))
  #list the blocks mapped to DMRs
  #print(in_blocks2DMR_df.block_id)

  return total_blocks, total_blocks2mr, total_blocks2DMR, in_blocks2DMR_df


def run(args):
  in_block_file=args.in_sortedBlock_file
  flank_region2block=args.in_flank_region2block
  methylation_file=args.in_dmr_file
  out_folder=args.out_file_folder
  dmr_min_cutoff=args.dmr_min_cutoff 

  print("Start map blocks to DMRs ")
  print('Input block file: ', in_block_file)
  print('Input DMR file: ', methylation_file)
  print('Add flank region on both sides of block: ', flank_region2block)
  print('Minimum cutoff for DMR: ', str(dmr_min_cutoff))

  in_block_df=pd.read_csv(in_block_file,sep='\t',header=None)
  #add flank regions in both side of block
  out_block_df=in_block_df.copy()
  out_block_df.columns=['chr','pos_start','pos_end','block_id']
  out_block_df.pos_start=out_block_df.pos_start-flank_region2block
  #out_block_df.pos_start[out_block_df.pos_start<0]=0
  out_block_df.loc[out_block_df.pos_start<0,'pos_start']=0
  out_block_df.pos_end=out_block_df.pos_end+flank_region2block
  out_block_file=in_block_file.replace('.bed','_flank'+str(flank_region2block)+'bp2block.bed')
  #print(out_block_file)
  out_block_df.to_csv(out_block_file,sep='\t',header=None,index=False)

  #input sorted MR/DMR bed file
  #methylation_file='out_smooth_dmr/24_chroms_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_2_isModTest_0__all_dmrRanking_top_0.95_minLogReg_proba_0.7.bed'

  #out parameter
  #out_folder='out_blocks_dmr/'
  #out_folder='out_blocks_dmr_single/'
  min_overlap=1E-9
  is_bed=1
  if not os.path.exists(out_folder):
      print('Create folder: ', out_folder)
      os.makedirs(out_folder)

  #map block to DMR/MR
  out =os.path.join( out_folder , os.path.basename(out_block_file).replace('.bed','') + '_vs_' + os.path.basename(methylation_file).replace('.bed','') + '.bed')
  if is_bed==1:
     #here os doesnot wait for completing
     #os.system('bedtools intersect -a ' + out_block_file + ' -b ' + methylation_file + \
     #             ' -wa -wb -f ' + str(min_overlap) + ' > ' + out)
     command= 'bedtools intersect -a ' + out_block_file + ' -b ' + methylation_file + \
                ' -wa -wb -f ' + str(min_overlap) + ' > ' + out
     #print(command)
     #use subprocess to wait for command complete
     cmd=subprocess.Popen(command,shell=True)
     cmd.communicate()

     #here for wait for multiple processes to complete
     #processes=[]
     #processes.append(subprocess.Popen(command,shell=True))
     #while len(processes)>0 :
     #    processes[0].wait()
     #    del processes[0]
  #print(out)

  #check how many blocks are mapped to DMR 
  #in_block_file contains original block position
  #out_block_file contains block positions with added flank regions e.g., 500bp
  in_sorted_block_file=out_block_file
  record_out_files=out
  #dmr_min_cutoff=0.7 
  results=count_blocks_not_mapped2dmr(in_sorted_block_file,record_out_files,dmr_min_cutoff)

  out_block2dmr_file = os.path.join(out_folder ,  \
           os.path.basename(in_block_file).replace('.bed','_' + \
           str(flank_region2block)+'flank_')  +str(dmr_min_cutoff)+'Proba_'+str(results[0]) + \
           'blocks_' + str(results[1]) + 'blocks2mr_' + str(results[2]) + 'blocks2dmr'+'.tsv' )

  out_df=results[3]
  out_df.to_csv(out_block2dmr_file,sep='\t',index=False)
  print("\n")
  print('Export mutation block mapped to DMR file at: ')
  print(out_block2dmr_file)


if __name__=='__main__':
  #input sorted block mutation bed file
  #block of minimum 2 patients
  #in_block_file='bp2_mussd_blocks_alireza/blocks_summary_sorted.bed'

  #block of minimum 1 patients
  #in_block_file='bp2_mussd_blocks_single/blocks_summary_block_position.bed'
  #flank_region2block=1000

  #input parameters
  args=my_parser(argparse.ArgumentParser('python map_block2dmr.py ')).parse_args()
  run(args)






