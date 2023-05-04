#this script is used to find block mutation patient ids
#based on bayespi-bar2 mussd output results for each mutation block 
#exec(open("find_block_patientID_april.py").read())
import pandas as pd
import os
import argparse
#this function is not tested after changing in May 7

def my_parser(parser):
  required=parser.add_argument_group('Required')
  required.add_argument('-in_BSfile','--in_block_summary_file',required=True,help='Path to block summary file exported by bpb3 mussd ')
  required.add_argument('-in_BSfolder','--in_block_folder', required=True, help='Path to blocks exported by bpb3 mussd that contains sequence, mutation and patient information ')
  return parser

def run(args):
 in_file=args.in_block_summary_file
 block_folder=args.in_block_folder

 print("Start to find block patient IDs ")
 print("In block summary file: ", in_file)
 print("In block folder: ", block_folder)

 #sumary of mutation block from BB2
 #in_file='bp2_mussd_blocks_single/blocks_summary.tsv'
 in_df=pd.read_csv(in_file,sep='\t')

 #insert a new column to dataframe
 in_df.insert(7,'patient_id',['']*in_df.shape[0])

 #there are a huge number of files at here if the number of mutation block is large 
 #block_folder='bp2_mussd_blocks_single/out_bp2_mussd/out_bp2/'
 record_patients=[]

 #there are two ways to do it , list append, or datafram insertion
 #it seems that list append more quick then the datafram insertion ??
 #datafram insertion used 30 minutes
 #list append used 
 for indx,row in in_df.iterrows():
    tmp_block=row.block_id
    tmp_file=os.path.join(block_folder,tmp_block+'.bed')
    tmp_df=pd.read_csv(tmp_file,sep='\t',header=None)
    #assume patient information in the last column of dataframe and extract it
    tmp_patients=tmp_df.iloc[:,-1].apply(lambda x: '_'.join(x.split('_')[-2:]))
    tmp_patientsID=','.join(tmp_patients.to_list())
    record_patients.append(tmp_patientsID)
    #add new element in the new column
    #in_df.loc[indx,'patient_id']=tmp_patientsID

 in_df['patient_id']=record_patients

 out_file=in_file.replace('.tsv','')+'_and_patientID2.tsv'
 print("\n")
 print("Export mutation block with patient ID at:")
 print(out_file)
 in_df.to_csv(out_file,sep='\t',index=False)


if __name__=='__main__':
 args=my_parser(argparse.ArgumentParser('python find_block_patientID.py ')).parse_args()
 run(args)


