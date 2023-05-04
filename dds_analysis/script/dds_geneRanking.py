#this script is used to rank final uniq gene list from dds_analysis in out_dds
import pandas as pd
import numpy as np
import argparse
#exec(open('dds_geneRanking.py').read())

def my_parser(parser):
  required = parser.add_argument_group('Required')
  required.add_argument('-inUqGene','--in_unique_gene_file', required=True, help='input a unique gene list file that exported from dds_analysis')
  required.add_argument('-inDEGFile','--in_DEG_file', required=True, help='input differential expressed gene file list that exported from BayesPI-BAR2')
  required.add_argument('-inDMRFile','--in_DMR_file',required=True, help='input differential methylation region file that exported drom dmr_analysis')
  optional = parser.add_argument_group('Optional , with default values')
  optional.add_argument('-inCutoff','--in_cutoff_pval4score',default= 0.5, metavar='', type=float, help='cutoff probability value for selecting top ranked genes, default = 0.5')
  optional.add_argument('-TSS', '--TSS_score', default=4.0, metavar='',type=float, help='Weight score for TSS, default= 4.0')
  optional.add_argument('-gene','--gene_score', default=2.0, metavar='',type=float, help='Weight score for gene, default= 2.0')
  optional.add_argument('-TES', '--TES_score', default=2.0, metavar='',type=float, help='Weight score for TES, default= 2.0')
  optional.add_argument('-dist', '--dist_score', default=1.0, metavar='',type=float, help='Weight score for 5 distance regions, default= 1.0')
  optional.add_argument('-enhancer', '--enhancer_score', default=3.0, metavar='',type=float, help='Weight score for enhancer, default= 3.0')
  return parser

def separate_patient_id(strs):
  tmp_strs=strs.split('~')
  out_strs=[]
  for ti in tmp_strs:
      out_strs.extend(ti.split(','))
  out_set=set(out_strs)
  return len(out_set)

def find_gene_pval(strs,in_diffExp_df):
    #return signifncat P-value of DEG
    tmp_pval=in_diffExp_df.loc[(in_diffExp_df.gene.str.strip().str.lower()==strs.strip().lower()),in_diffExp_df.columns[1]]
    return tmp_pval.to_numpy()[0]

def find_dmr_pval(strs,selected_dmr_df2):
  #return maximum probabilty of DMRs
  if str(strs) != 'nan':
    tmp_strs=strs.split('~')
    tmp_strs2=[]
    for ts in tmp_strs:
        tmp_strs2.extend(ts.split(','))
    tmp_strs=list(set(tmp_strs2))
    tmp_pval=[]
    for ts in tmp_strs:
         if ts.lower() != 'nan':
           #print(selected_dmr_df.loc[selected_dmr_df['dmr_id'] ==ts.strip() ,'logReg_pval'  ].to_numpy())
           tmp_pval.append(selected_dmr_df2.loc[selected_dmr_df2['dmr_id'] ==ts.strip() ,'logReg_pval'  ].to_numpy()[0] )
    #print(tmp_pval)
    if tmp_pval==[]:
       max_pval=0
    else:
       max_pval=np.max(tmp_pval)
  else:
    max_pval=0
  #
  return max_pval

def find_geneType_score(strs,weight_scores):
  #here return maximum geneType score
  tmp_strs=strs.split('~')
  tmp_strs=list(set(tmp_strs))
  if 'TSS' in tmp_strs:
      out_score= weight_scores[0]
  elif 'gene' in tmp_strs:
      out_score=weight_scores[1]
  elif 'enhancer' in tmp_strs:
      out_score=weight_scores[2]
  elif 'TES' in tmp_strs:
      out_score=weight_scores[3]
  elif '5dist' in tmp_strs:
      out_score=weight_scores[4]
  else:
      out_score=0
  return out_score

def find_geneType_meanScore(strs, weight_scores):
  #caclulate average geneType score for each gene based
  #on its annotated muation blocks in various genomic regions
  gi=strs
  [gt,eh]=gi.split(':')
  gt_list=gt.split('~')
  eh_list=eh.split('~')
  joined_list=list(map(':'.join,zip(gt_list,eh_list)))
  total_len=len(joined_list)
  uq_list=list(set(joined_list))
  uq_score=[]
  uq_count=[]
  for ul in uq_list:
      if 'TSS' in ul:
          out_score=weight_scores[0]
      elif 'gene' in ul:
          out_score=weight_scores[1]
      elif 'enhancer' in ul:
          out_score=weight_scores[2]
      elif 'TES' in ul:
          out_score=weight_scores[3]
      elif '5dist' in ul:
          out_score=weight_scores[4]
      else:
          out_score=0
      uq_score.append(out_score)
      uq_count.append(np.count_nonzero(np.array(joined_list)==ul))
  out_mscore=np.sum(np.array(uq_score)*np.array(uq_count))/total_len
  return out_mscore

def max_min_normalization(in_val):
   #maximum and minimum normalization of data
   max_val=np.max(in_val)
   min_val=np.min(in_val)
   if max_val==min_val and len(in_val)==1:
     out_val=[1.0]
   else:
     out_val=(in_val-min_val)/(max_val-min_val)
   return out_val


def check_enhancers(sorted_in_df, cutoff_pval):
  #check enhancer within block
  s2=sorted_in_df[sorted_in_df.mean_scores>=cutoff_pval].copy()
  #s3=s2[s2.geneType_enhancers2.apply(lambda x: 'enhancer' in x)]
  #
  sid=s2.block_id.apply(lambda x: str(x).split('~'))
  sen=s2.enhancers.apply(lambda x: str(x).split('~'))
  sty=s2.gene_type.apply(lambda x: str(x).split('~'))
  #
  all_pairs=[]
  for i in range(0,sid.shape[0]):
    all_pairs.extend(list(zip(sid[i], sen[i],sty[i])))
  #
  all_blocks=set(all_pairs)
  num_of_blocks=len(all_blocks)
  #
  is_enhancer=0
  is_tss_tes_or_gene=0
  for ii in list(all_blocks):
     if 'enhancer' in ii:
        is_enhancer +=1
     elif ('TSS' in ii) or ('TES' in ii) or ('gene' in ii) : 
        is_tss_tes_or_gene +=1
  print('Cutoff value is: '+ str(cutoff_pval))
  if num_of_blocks>0:
    print('Percentage of blocks in enhancers: '+ str(is_enhancer/num_of_blocks) +'; Percentage of blocks in TSS/TES/Gene: ' + str(is_tss_tes_or_gene/num_of_blocks))
  else:
    print('Percentage of blocks in enhancers/number of blocks: '+ str(num_of_blocks) +'; Percentage of blocks in TSS/TES/Gene, number of blocks: ' + str(num_of_blocks))


def run(args):
 in_file=args.in_unique_gene_file
 in_diffExp_file=args.in_DEG_file
 in_dmr_file=args.in_DMR_file
 score_cutoff=args.in_cutoff_pval4score
 #weight scores to regions
 weight_scores=[args.TSS_score, args.gene_score, args.enhancer_score, args.TES_score, args.dist_score]
 print('Weight scores for TSS, gene, enhancer, TES, 5Dist are: ', [str(i) for i in weight_scores])

 in_df=pd.read_csv(in_file,sep='\t')
 in_diffExp_df=pd.read_csv(in_diffExp_file,sep='\t')
 in_columns=in_diffExp_df.columns.to_list()
 in_columns[0]='gene'
 in_diffExp_df.columns=in_columns
 in_dmr_df=pd.read_csv(in_dmr_file,sep='\t',header=None)
 in_dmr_df.columns=['chrom','start_pos','end_pos','dmr_info','logReg_pval']
 #logReg_cutoff=0.7
 #selected_dmr_df=in_dmr_df[in_dmr_df.logReg_pval>=logReg_cutoff].copy()
 #selected_dmr_df=selected_dmr_df.reset_index(drop=True)
 in_dmr_df['dmr_id']=in_dmr_df.dmr_info.apply(lambda x: ':'.join(x.split(':')[0:2]))
 #
 in_df['geneType_enhancers2']=in_df[['gene_type','enhancers']].apply(lambda x: ':'.join(x.astype(str)),axis=1)
 in_df['geneType_meanScore']=in_df.geneType_enhancers2.apply(lambda x: find_geneType_meanScore(x,weight_scores))
 #
 #in_df['geneType_enhancers']=in_df[['gene_type','enhancers']].apply(lambda x: '~'.join(x.astype(str)),axis=1)
 in_df['num_of_patients']=in_df.patient_id.apply(separate_patient_id)
 in_df['diffExp_pval']=in_df.gene_name.apply(lambda x: find_gene_pval(x,in_diffExp_df))
 in_df['dmr_pval']=in_df.new_mr_sites.apply(lambda x: find_dmr_pval(x, in_dmr_df))
 #
 #tss=4 gene=4, enhancer=3, TES=2 , 5dist=1
 #in_df['geneType_score']=in_df.geneType_enhancers.apply(lambda x: find_geneType_score(x)) 
 #in_df['normalized_geneType_score']=list(max_min_normalization(in_df['geneType_score'].to_numpy()))
 #
 in_df['normalized_geneType_meanScore']=list(max_min_normalization(in_df['geneType_meanScore'].to_numpy()))
 in_df['normalized_num_of_patients']=list(max_min_normalization(in_df['num_of_patients'].to_numpy()))
 in_df['normalized_dmr_pval']=list(max_min_normalization(in_df['dmr_pval'].to_numpy()))
 in_df['normalized_diffExp_pval']=list(max_min_normalization(in_df['diffExp_pval'].apply(lambda x: np.abs(np.log10(x)))))
 #
 in_df['median_scores']=in_df[['normalized_geneType_meanScore','normalized_num_of_patients','normalized_dmr_pval','normalized_diffExp_pval']].median(axis=1).to_list()
 in_df['mean_scores']=in_df[['normalized_geneType_meanScore','normalized_num_of_patients','normalized_dmr_pval','normalized_diffExp_pval']].mean(axis=1).to_list()
 in_df.sort_values(by='mean_scores')
 sorted_in_df=in_df.sort_values(by='mean_scores',ascending=False).copy()
 sorted_in_df=sorted_in_df.reset_index(drop=True)
 selected_in_df=sorted_in_df[sorted_in_df.mean_scores>=score_cutoff].copy()

 #remove some of columns from export file 
 selected_in_df.drop(['geneType_enhancers2','normalized_geneType_meanScore','normalized_num_of_patients','normalized_dmr_pval','normalized_diffExp_pval'],axis='columns',inplace=True)
 out_file=in_file+'_selectedGenes_gt_'+str(score_cutoff)+'.txt'
 out_file=out_file.replace('.tsv','')
 print('Export : ', out_file)
 selected_in_df.to_csv(out_file,sep='\t',index=False)
 check_enhancers(sorted_in_df, score_cutoff)
 return sorted_in_df 

if __name__=='__main__':
  args= my_parser(argparse.ArgumentParser('python dds_geneRanking.py ')).parse_args()
  run(args)





 
