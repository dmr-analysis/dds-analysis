#this script is used to filter blocks with neither DMR and DEG
import pandas as pd
import numpy as np
import argparse

#exec(open('filter_blocks.py').read())

def my_parser(parser):
  required = parser.add_argument_group('Required')
  required.add_argument('-in_BKfile','--in_combined_DmrDegBlock_file',required=True, 
                 help='A bpb3 mutation block file combined with DMR DEG and genomic reigon information that is exported by dds_analysis combine_dmr_deg2block.' \
	' Here, we selecte blocks that are satisifed by condition 1: blocks associated with both DMR (not null) and DEG/enhancer (not null), and condition 2: blocks associated with either' \
	' DMR (not null) or DEG/enhancer (not null) ; where DMR and DEG/enhancer are already defined in previous analysis, if not then they are marked by nan in input file --in_combined_DmrDegBlock_file')

  optional= parser.add_argument_group('Optional arguments with default value')
  optional.add_argument("--not_use_enhancer", help='Whether use/consider enhancer information during the filtering , default is True, if use this option then will not use enhancer information during the filtering.',
                        action="store_true")
  return parser

def replace_nan(x):
   ''' replace element with full (nan,nan,...) in dataframe by a signle NaN
   '''
   tmp_x=x.split(',')
   out=[]
   out2=[]
   for ii in tmp_x:
      if ii.lower() == 'nan':
         out.append(1)
      else:
         out2.append(ii)

   if len(out)==len(tmp_x):
      out3=np.nan
   else:
      out3=x
   return out3

def sort_Xlist_by_Ylist(X,Y):
   #sort list X by list Y
   sorted_X_by_Y=sorted(zip(Y,X))
   sorted_x=[x for _,x in sorted_X_by_Y ]
   sorted_y=[y for y,_ in sorted_X_by_Y ]
   return sorted_x,sorted_y

def find_uq_based_on_list1(list1,list2):
  #find unique value in np_array1 and extract the corresponding firt value in np_array2 too
  uq_list1=list(set(list1))
  uq_list2=[]
  for i in uq_list1:
      uq_list2.append(list2[np.where(list1==i)][0] )
  sorted_uq_list1,sorted_uq_list2=sort_Xlist_by_Ylist(uq_list1,uq_list2)
  return sorted_uq_list1,sorted_uq_list2

def remove_gene_pval_eq_nan(features, in_block_df0):
  #function to remove gene with p-values is nan 
  #features=['TSS','TES','gene','5dist']
  in_block_df=in_block_df0.copy()
  for fe in features:
    tmp_col1=fe
    tmp_col2='deg_p_value2'+tmp_col1
    record_name=[]
    record_pval=[]
    for idx,rows in in_block_df[[tmp_col1,tmp_col2]].iterrows():
      if rows.isnull().sum()<2:
         tmp_name=np.array(str(rows[tmp_col1]).replace(' ','').split(','))
         tmp_pval=np.array(str(rows[tmp_col2]).split(','),float)

         list1=np.array(list(tmp_name[~np.isnan(tmp_pval)]))
         list2=np.array(list(np.array(tmp_pval[~np.isnan(tmp_pval)],dtype=str)))
         uq_list1,uq_list2=find_uq_based_on_list1(list1,list2)
         #print(uq_list1,uq_list2) 
         #new_name=', '.join(list(tmp_name[~np.isnan(tmp_pval)]))
         #new_pval=', '.join(list(np.array(tmp_pval[~np.isnan(tmp_pval)],dtype=str)))
         #print(new_name,new_pval)
         new_name=','.join(uq_list1)
         new_pval=','.join(uq_list2)

         record_name.append(new_name)
         record_pval.append(new_pval)
      else:
         record_name.append(rows[tmp_col1])
         record_pval.append(rows[tmp_col2])
    new_col1='filtered_'+tmp_col1
    new_col2='filtered_'+tmp_col2
    in_block_df[new_col1]=record_name
    in_block_df[new_col2]=record_pval
  return in_block_df

def find_feature_index(feature_array,feature_str):
  '''find index of a feature string that contained in an numpy array'''
  idx2_feature=np.where(np.char.find(feature_array,feature_str)>=0)[0]
  return idx2_feature

def generate_filter_conditions(columns_df,not_use_enhancer):
  '''check the existence of features (logReg, TSS, TES, gene, 5dist, and enhancer) in a numpy array , then 
     then generate the filtering conditions for DEG_or_DMR, DEG_and_DMR, DEG_and_DMR_and_block_in_enhancer
     based on existed features.
  '''
  idx2_logReg=find_feature_index(columns_df,'_logReg_proba')
  idx2_TSS=find_feature_index(columns_df,'_p_value2TSS')
  idx2_TES=find_feature_index(columns_df,'_p_value2TES')
  idx2_gene=find_feature_index(columns_df,'_p_value2gene')
  idx2_5dist=find_feature_index(columns_df,'_p_value25dist')
  idx2_enhancer=find_feature_index(columns_df,'_p_value2enhancers')

  if len(idx2_logReg)>0:
    filter1_logReg = '(~in_block_df.mr_logReg_proba.isnull()) '
  else:
    filter1_logReg=''

  if len(idx2_TSS)>0:
    filter1_TSS =  ' (~in_block_df.deg_p_value2TSS.isnull()) '
  else:
    filter1_TSS=''

  if len(idx2_TES)>0:
    filter1_TES = ' (~in_block_df.deg_p_value2TES.isnull()) '
  else:
    filter1_TES=''

  if len(idx2_gene)>0:
    filter1_gene = '(~in_block_df.deg_p_value2gene.isnull()) '
  else:
    filter1_gene=''

  if len(idx2_5dist) >0:
    filter1_5dist = ' (~in_block_df.deg_p_value25dist.isnull()) '
  else:
    filter1_5dist=''

  if len(idx2_enhancer)>0:
    filter1_enhancer=' (in_block_df.enhancers==\'enhancer\') '
  else:
    filter1_enhancer=''

  #DMR or DEG
  filter_str=''
  filter2_str=''
  filter3_str=''
  if len(filter1_logReg)>0:
    filter_str = filter1_logReg

  if len(filter1_TSS)>0:
    filter_str += ' | '+ filter1_TSS
    filter2_str = filter1_TSS
    filter3_str = filter1_TSS

  if len(filter1_TES)>0:
    filter_str += ' | ' + filter1_TES
    filter2_str += ' | ' +filter1_TES
    filter3_str += ' | ' + filter1_TES

  if len(filter1_gene)>0:
    filter_str += ' | ' + filter1_gene
    filter2_str += ' | ' + filter1_gene
    filter3_str += ' | ' + filter1_gene

  if len(filter1_5dist)>0:
    filter_str += ' | ' + filter1_5dist
    filter2_str += ' | ' + filter1_5dist
    if len(filter1_enhancer)>0 : #and not not_use_enhancer :
       #print('Consider Enhancer information during filtering')
       filter3_str += '| ( '+  filter1_5dist + ' & ' + filter1_enhancer + ' )'
    else:
       #print('Enhancer is not considered during filtering')
       filter3_str +=  ' | ' + filter1_5dist

  #DMR or DEG
  condition2filter= '( ' + filter_str + ' )'
  #DMR and DEG
  condition2filter2= '('+filter1_logReg + ' & (' + filter2_str + '))'
  #DMR and DEG and block in enhancer
  condition2filter3= '('+filter1_logReg + ' & (' + filter3_str + '))'
  final_features=[len(filter1_TSS) , len(filter1_TES), len(filter1_gene), len(filter1_5dist)]
  return condition2filter, condition2filter2, condition2filter3, final_features


def run(args):
  in_block_file=args.in_combined_DmrDegBlock_file

  print("Start to filter blocks based on DMR , DEG information")
  print("For example, blocks with the two sides flank regions do not overlap to DMR and no DEG in assigned genes will be removed.")
  print("Only blocks with two-sides of flank regions overlap to DMR or has a DEG in assigend gene will be exported.")
  print("In combined DMR DEG , and Block information file: ", in_block_file)
  #input data
  #minimum 2 patients
  #in_block_file='blocks_summary_sorted_500flank_0.7Proba_176blocks_73blocks2mr_30blocks2dmr_deg_info.tsv'
  #minimum 1 pateint
  #in_block_file='out_blocks_gene_single/blocks_summary_block_position_500flank_0.7Proba_66868blocks_28049blocks2mr_9478blocks2dmr_deg_info.tsv'
  #in_block_file='out_blocks_gene_single/blocks_summary_block_position_0flank_0.7Proba_66868blocks_13143blocks2mr_4604blocks2dmr_deg_info.tsv'

  in_block_df=pd.read_csv(in_block_file,sep='\t',dtype={'mr_logReg_proba': str, 'deg_p_value2enhancers' :str, 'mutation_distribution': str})
  in_block_df0=in_block_df.copy()
    
  #Filtering conditions
  #check available features before generating featue filtering conditins 
  columns_df=np.array(in_block_df.columns.to_list())

  #first to replace an element of p.values with all nan,nan by a single nan in dataframe
  idx2_logReg=find_feature_index(columns_df,'_logReg_proba')
  idx2_TSS=find_feature_index(columns_df,'_p_value2TSS')
  idx2_TES=find_feature_index(columns_df,'_p_value2TES')
  idx2_gene=find_feature_index(columns_df,'_p_value2gene')
  idx2_5dist=find_feature_index(columns_df,'_p_value25dist')
  idx2_enhancer=find_feature_index(columns_df,'_p_value2enhancers')

  if len(idx2_logReg)>0:
    in_block_df[columns_df[idx2_logReg][0] ]=eval('in_block_df.'+ columns_df[idx2_logReg][0] +'.astype(str).apply(lambda x: replace_nan(x))')
  if len(idx2_TSS)>0:
    in_block_df[columns_df[idx2_TSS][0] ]=eval('in_block_df.'+ columns_df[idx2_TSS][0] +'.astype(str).apply(lambda x: replace_nan(x))')
  if len(idx2_TES)>0:
    in_block_df[columns_df[idx2_TES][0] ]=eval('in_block_df.'+ columns_df[idx2_TES][0] +'.astype(str).apply(lambda x: replace_nan(x))')
  if len(idx2_gene)>0:
    in_block_df[columns_df[idx2_gene][0] ]=eval('in_block_df.'+ columns_df[idx2_gene][0] +'.astype(str).apply(lambda x: replace_nan(x))')
  if len(idx2_5dist)>0:
    in_block_df[columns_df[idx2_5dist][0] ]=eval('in_block_df.'+ columns_df[idx2_5dist][0] +'.astype(str).apply(lambda x: replace_nan(x))')
  if len(idx2_enhancer)>0:
    in_block_df[columns_df[idx2_enhancer][0] ]=eval('in_block_df.'+ columns_df[idx2_enhancer][0] +'.astype(str).apply(lambda x: replace_nan(x))')

  #filter condtions for DMR_or_DEG, DMR_and_DEG, DMR_and_DEG_and_block in enhancer, respectively
  condition2filter, condition2filter2, condition2filter3, final_features=generate_filter_conditions(columns_df,args.not_use_enhancer)

  #DMR or DEG
  #condition2filter=( (~in_block_df.mr_logReg_proba.isnull()) | (~in_block_df.deg_p_value2TSS.isnull()) |
  #                (~in_block_df.deg_p_value2TES.isnull()) | (~in_block_df.deg_p_value2gene.isnull())  |
  #                (~in_block_df.deg_p_value25dist.isnull()) )

  #DMR and DEG
  #condition2filter2=( (~in_block_df.mr_logReg_proba.isnull()) & ( (~in_block_df.deg_p_value2TSS.isnull()) |
  #                (~in_block_df.deg_p_value2TES.isnull()) | (~in_block_df.deg_p_value2gene.isnull())  |
  #                (~in_block_df.deg_p_value25dist.isnull()) ))

  #DMR and DEG and block in enhancer
  #condition2filter3=( (~in_block_df.mr_logReg_proba.isnull()) & 
  #                  (( ~in_block_df.deg_p_value2TSS.isnull()) | (~in_block_df.deg_p_value2TES.isnull()) | (~in_block_df.deg_p_value2gene.isnull())  |
  #                    ((~in_block_df.deg_p_value25dist.isnull()) & (in_block_df.enhancers=='enhancer')) )) 

  #Filtering blocks by conditions
  #print(condition2filter)
  #print(condition2filter2)
  dd=in_block_df0[eval(condition2filter)].copy()
  dd2=in_block_df0[eval(condition2filter2)].copy()
  dd3= in_block_df0[eval(condition2filter3)].copy()
  #dd[dd['5dist'].apply(lambda x : 'BCL6' in str(x))]['new_mr_sites']

  #generate exist features
  tmp_features=['TSS','TES','gene','5dist']
  tmp_col2=[['filtered_TSS','filtered_deg_p_value2TSS'],
            ['filtered_TES','filtered_deg_p_value2TES'],
            ['filtered_gene','filtered_deg_p_value2gene'],
            ['filtered_5dist', 'filtered_deg_p_value25dist']]
  features=[]
  col2=[]
  for i in range(0,len(tmp_features)):
      if final_features[i]>0:
         features.append(tmp_features[i])
         col2 += tmp_col2[i]

  #remove genes with p-value equals nan, and get unique gene name
  dd_filtered= remove_gene_pval_eq_nan(features, dd)
  dd2_filtered= remove_gene_pval_eq_nan(features, dd2)        
  dd3_filtered= remove_gene_pval_eq_nan(features, dd3)        
  
  #export dd
  col1=['block_id' ,'chrom'   ,'start_pos' ,'end_pos' ,'number_of_mutations','number_of_patients','mutation_distribution','patient_id','new_mr_sites','mr_logReg_proba']
  col3=['enhancers','deg_p_value2enhancers']
  #col2=['filtered_TSS','filtered_deg_p_value2TSS','filtered_TES','filtered_deg_p_value2TES','filtered_gene','filtered_deg_p_value2gene','filtered_5dist', 'filtered_deg_p_value25dist' ]
  out_file1=in_block_file.replace('.tsv','_filtered_DMR_or_DEG.tsv')
  new_col_names=col1+col2+col3
  out_dd_df=eval('dd_filtered[[\'' +'\',\''.join(new_col_names)+ '\']].copy()')
  out_dd_df.to_csv(out_file1,sep='\t',index=False)
  print("Export mutation blocks are associated with either DMR or DEG ")
  print(out_file1)

  #export dd2 
  out_file=in_block_file.replace('.tsv','_filtered_DMR_and_DEG.tsv')
  out_dd2_df=eval('dd2_filtered[[\'' +'\',\''.join(new_col_names)+ '\']].copy()')
  out_dd2_df.to_csv(out_file,sep='\t',index=False)
  print("Export mutation blocks are associated with both DMR and DEG")
  print(out_file)
 

if __name__=='__main__':
  args= my_parser(argparse.ArgumentParser('python filter_blocks.py ')).parse_args()
  run(args)



