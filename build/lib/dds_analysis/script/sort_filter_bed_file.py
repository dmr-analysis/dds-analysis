#this script is used to filter non standard chroms then sort the bed file
import pandas as pd
import argparse
from dmr_analysis.script.script_high.map_dmr2chromSegment import read_bed_file_and_export2type_bed_file, sort_bed_file_df

def my_parser(parser):
  required_name=parser.add_argument_group("required arguments")
  required_name.add_argument("-inFile","--in_file_name", help="Input file name", required=True, type=str)

  optional_name=parser.add_argument_group("Optional arguments, have default value")
  optional_name.add_argument("-noLength","--not_cacl_length", help="Is calculate length between start and end postion, default=True", 
                              action="store_false",default=True)
  optional_name.add_argument("-noFilter","--not_filter_chr",help="Is filter nonstandard chromosomes before sorting, default=True",
                               action="store_false", default=True)
  optional_name.add_argument("-noHuman", "--not_human_genome", help="Is human genome, default=True",
                               action="store_false", default=True)
  return parser

def sort_and_filter(in_file, isLength, isFilter_chr, ishuman):
  '''sort and filter BED file '''
  tmp_in_pd=pd.read_csv(in_file, sep='\t',header=None)
  columns_name=['chrs','start_pos','end_pos','type','len_of_region']
  if isLength :
    print('Calculate region legnth')
    if tmp_in_pd.shape[1]==4:
       tmp_in_pd.insert(4,'4',list(tmp_in_pd.iloc[:,2]-tmp_in_pd.iloc[:,1]),True)
    elif tmp_in_pd.shape[1]==3:
       tmp_in_pd.insert(3,'3',list(tmp_in_pd.iloc[:,2]-tmp_in_pd.iloc[:,1]),True)
       tmp_in_pd.insert(4,'4',list(tmp_in_pd.iloc[:,2]-tmp_in_pd.iloc[:,1]),True)
    else:
       print('Input bed file column size >4 ')
       pass
    tmp_in_pd.columns=columns_name
    in_pd=tmp_in_pd.copy()
  else:
    print('Join from fifth columns to the end of dataframe into one column')
    tmp_joined_columns=tmp_in_pd.iloc[:,4:].apply(lambda x: '~'.join(x.astype(str)),axis=1).to_list()
    in_pd=tmp_in_pd.iloc[:,0:4].copy()
    in_pd.insert(4,'4',tmp_joined_columns, True)
    in_pd.columns=columns_name

  #isFilter_chr=True
  if isFilter_chr:
    print('Filter nonstandard Chromosomes before sorting')
    in_pd2=in_pd[~in_pd.chrs.str.contains('_')].copy()
  else:
    in_pd2=in_pd.copy()

  #ishuman=True
  tmp_sorted_pd=sort_bed_file_df(in_pd2,columns_name,ishuman)
  if not isLength:
    print('Expand sorted dataframe to original form')
    tmp_df=pd.DataFrame(tmp_sorted_pd.iloc[:,4].str.split('~').to_list())
    sorted_pd=pd.concat([tmp_sorted_pd.iloc[:,0:4],tmp_df],axis=1)
  else:
    sorted_pd=tmp_sorted_pd.copy()

  print('Input file: ', in_file)
  if '.bed' in in_file:
     out_file=in_file.replace('.bed','_sorted.bed')
  elif '.bdg' in in_file:
     out_file=in_file.replace('.bdg','_sorted.bdg')
  else:
     print('Not find bed or bdg in file: ', in_file)

  print('Export sorted bed files at: ', out_file )
  sorted_pd.iloc[:,-1]=sorted_pd.iloc[:,-1].replace('\t','')
  sorted_pd.to_csv(out_file,sep='\t',index=False,header=None)
  return sorted_pd,out_file

def run(args):
  in_file=args.in_file_name
  #in_file='out/43-A7v/43-A7v_S43_spmr_treat_pileup.bdg'
  #in_file='out/43-A7v/43-A7v_S43_spmr_peaks.no_encode_anomalies.main_chromosomes_only.narrowPeak.bed'

  #True input file with four columns and need to calculate region length
  #otherwise, combine all 5 plus columns into one column 
  #isLength=True
  isLength=args.not_cacl_length
  isFilter_chr=args.not_filter_chr
  ishuman=args.not_human_genome
  sorted_df,out_file_name=sort_and_filter(in_file,isLength, isFilter_chr, ishuman)

if __name__=='__main__':
 #in_file='hg38_genome_amp/hg38-blacklist.v2.bed'
 args=my_parser(argparse.ArgumentParser('python sort_filter_bed_file.py ')).parse_args()
 run(args)


