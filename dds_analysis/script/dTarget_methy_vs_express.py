#this script is used to do methylation vs. expression analysis based on generated data , the method is adopted from paper
#https://pubmed.ncbi.nlm.nih.gov/25994056/
#for example, to extract methylation levels for MRs
import os
import re
import pandas as pd
import glob
import multiprocessing as mp
from dmr_analysis.script.script_high.dmr_utility import divide_data2chunk
import numpy as np
from scipy import stats
import json
#from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
import matplotlib.pyplot as plt
import matplotlib as mlp
from dmr_analysis.script.dmr_exportData import get_mr_in_blocks
import warnings
import argparse
warnings.filterwarnings("ignore")
np.seterr(all='ignore')
mlp.use('Agg')

def my_parser(parser):
  required= parser.add_argument_group('Required')
  required.add_argument('-inGeneMRfile','--in_gene_mr_file_from_dds', required=True, type=str,help='Input file name for gene-MR relationships that exported by dds_analysis check_block_gene_inTAD, \
                                             file format: gene_name; gene_type; block_id; new_mr_sites; patients; isTAD; enhancers; patient_id')
  required.add_argument('-inGeneEXPfile','--in_gene_expression_file_from_bpb3', required=True, type=str, help='Input file name for gene expression profiles that exported by bpb3 differential_expression,\
                                             file format: gene_name; P-values, sample1, sample2 , ....  ')
  required.add_argument('-inMRfolder','--in_methylation_folder_from_dmr', required=True, type=str, help='Input file folder path for methylation data, the data are extracted by \
                                              dmr_analysis dmr_exportData based on predifined blocks or regions')
 
  required.add_argument('-inBackgroundList','--in_mr_background_file_list', required=True, type=str, help='A list of files will be used in background sampling, one file one line')

  optional= parser.add_argument_group('Optional, with default values ')
  optional.add_argument('-cutoff','--expected_pval_cutoff', default=0.05 , metavar='', type=float, help='P-value cutoff for expected P-value after sampling, default= 0.05')
  optional.add_argument('-reg_cutoff','--regression_pval_cutoff', default=0.05, metavar='', type=float, help='P-value cutoff for regression model when fitting expression data against methylation data, default=0.05')
  optional.add_argument('-totalSamples','--total_samples_in_random_selection', default=1000, metavar='', type=int, help='Total number of samples will draw from random backgroun samples, default=1000')
  optional.add_argument('-numOfprocesses', '--number_of_parallel_processes', default=10, metavar='', type=int, help='Number of parallel processes will be used in calculation, default= 10 ')
  optional.add_argument('-expStartCol','--gene_expression_start_column', default=2, metavar='',type=int, help='Start column number (e.g., 0,1,2,...) for samples in gene expression profile file, \
                  default =2 for bpb3 exported file')
  optional.add_argument('-sampleName','--in_sampleName_for_replace', default='sample_name4replace.tsv', type=str, help='Input file name for rename the old sample names in gene expression data, first column old_name is old sample name, second \
                                                column new_name is the new name that used for replacement, default file name is sample_name4replace.tsv which can be an empty file if no replacement is needed for sample names.')
  optional.add_argument('-pathDepth','--in_diffExp_pathDepth4sampleName',default=12, type=int, help='file path depth for column sample name, for example, if columna name is /path/folder/sample_name/.., then the \
                                                                    depth of sample name is 3, in other words, it is the number of "/" before the sample name, default=12')
  optional.add_argument('-expTAB','--is_tab_delimated_file4diffExp', help="Whether input differential expression gene file is a bpb3 differential_expression exported file format or a common tab delimated file format \
                                                           (e.g., the first column is the gene name and the rest of columns is gene expression data in samples) , default=False, if use this optional then input gene expression file \
                                                            will be treated as common tab delimated file.", action="store_true" )
  optional.add_argument('-mrTAB','--is_tab_delimated_file4MR', help='Whether input MR regions files is in dds_analysis check_block_gene_inTAD exported format or a common tab delimated file format \
                                                           (e.g., the first column is the gene_name , then is the gene_type, new_mr_sites), decault= False, if use this optional then input MR region file will \
                                                            be treated as common tab delimated file.', action="store_true")
  optional.add_argument('-no_dmr_format','--not_dmr_analysis_format', help='Input MR data folder is in dmr_analysis exprt data format, default= False, if use this optional then the data folder format is not in dmr_analysis export format', action="store_true")
  optional.add_argument('-outName','--output_file_name', help='Prefix file name for out data exported by this program, default=output_result_of_*', default='output_result_of_',type=str)
  optional.add_argument('-mrColName','--mr_id_column_name',help='Column name of MR ids where contains string ~ joined mr IDs such as chr1:mr1~chr1:mr11, default=new_mr_sites, if mr IDs in other column then pleaes input its new column name',
                                    default='new_mr_sites',type=str )
  optional.add_argument('-output_path','--output_files_path', help='Output files path for data exported by program, default=./ that export files in local directory', default='./', type=str)
  return parser

  
def read_geneexp_data(in_gene_expression, gene_exp_start_col,pathDepth4sampleName):
  '''input gene expression data from bpb3 exported differentially expressed genes'''
  #in_gene_expression='new_diffExp_genes/differentially_expressed_genes.txt'
  in_geneexp_df=pd.read_csv(in_gene_expression,sep='\t')
  tmp_exp_cols=in_geneexp_df.columns[gene_exp_start_col:].to_list()
  tmp_exp_sample_cols=[ i.split('/')[pathDepth4sampleName].replace('.just_counts.tsv_log_of_quantile_normalized_RPKM','') for i in tmp_exp_cols]
  return in_geneexp_df, tmp_exp_cols, tmp_exp_sample_cols

def combine_mr_and_expression_data(tmp_mr_file, tmp_out_folder , in_geneexp_data, exp2mr_normal_id,tmp_exp_sample_cols,
                                  mr_start_col=1  ):
   '''link gene to its methylation and expression data'''
   #here the first three input are for mr data,
   #read MR data
   #tmp_chr=tmps[0]
   #tmp_mr_id=tmps[1]
   #tmp_out_folder='out_mr_in_block/mr_in_block_out'
   #fine gene for its expression and the methylation levels

   #read MR methylation data
   tmp_mr_data_folder=tmp_out_folder
   #tmp_mr_file='_'.join(tmps)+'_raw.dat.gz'
   not_find_files=[]
   mean_methy_data=[]
   exp_data=[]
   tmp_mr_idx=[]
   tmp_exp_idx=[]
   if os.path.exists(os.path.join(tmp_mr_data_folder,tmp_mr_file)):
         #print(tmp_mr_file)
         tmp_methylation_data=pd.read_csv(os.path.join(tmp_mr_data_folder,tmp_mr_file),sep='\t',compression='gzip')
         tmp_mean_methylation=tmp_methylation_data.iloc[:,mr_start_col:].mean().to_numpy()/100
         #tmp_expression_data=in_geneexp_df[in_geneexp_df['#gene']==ki].iloc[:,geneexp_start_col:].to_numpy()
         #if tmp_expression_data.shape[0]==0:
         #   print('No expression data find, use zero instead')
         #   tmp_expression_data=np.zeros((1,tmp_expression_data.shape[1]))
         tmp_expression_data=in_geneexp_data

         #change columns for both methylation and expression data
         #jbw may
         #simple replace
         #tmp_mr_cols=[ i.replace('meth_','') for i in tmp_methylation_data.columns[1:].to_list()]
         #tmp_exp_cols=[i.replace('meth_','') for i in tmp_exp_sample_cols]
 
         #replace with wildcard
         #tmp_mr_cols=[ re.sub(r'meth(.?)_','',i) for i in tmp_methylation_data.columns[1:].to_list()]
         #tmp_exp_cols=[re.sub(r'meth(.?)_','',i) for i in tmp_exp_sample_cols]
       
         #jbw may
         #no replace and assumen Gene Expression label will use the same label as the methylation samples!
         tmp_mr_cols=[ i for i in tmp_methylation_data.columns[1:].to_list()]
         tmp_exp_cols=[ i for i in tmp_exp_sample_cols]

         #print(tmp_exp_cols)
         #print(tmp_mr_cols)
         #bug here if find a gene with two record in expression profiles??
         #print(tmp_expression_data.shape)
         tmp_expression_data=tmp_expression_data.reshape(len(tmp_exp_cols),)

         #extract samples appear in both exp and mr
         new_exp_cols=[]
         if len(exp2mr_normal_id)>0:
           for tmi in tmp_exp_cols:
             if tmi in exp2mr_normal_id.keys() :
                new_exp_cols.append(exp2mr_normal_id[tmi])
             else:
                new_exp_cols.append(tmi)
         else:
           new_exp_cols=tmp_exp_cols

         #find correct column index for samples in mr and geneexp
         loop=0
         tmp_exp_idx=[]
         tmp_mr_idx=[]
         for ni in new_exp_cols:
           mr_idx=np.where(np.char.find(tmp_mr_cols,ni)>=0)
           if len(mr_idx[0])>0:
             tmp_mr_idx.append(mr_idx[0][0])
             tmp_exp_idx.append(loop)
           loop +=1
         tmp_mr_cols=np.array(tmp_mr_cols)
         new_exp_cols=np.array(new_exp_cols)

         #sort samples in both methy and expresion with the same order before calculating mean methylation in the same MR
         mean_methy_data=tmp_mean_methylation[tmp_mr_idx]
         exp_data=tmp_expression_data[tmp_exp_idx]
   else:
         #print('Not find: ', os.path.join(tmp_mr_data_folder, tmp_mr_file))
         not_find_files.append(tmp_mr_file)
   return mean_methy_data, exp_data,  not_find_files,tmp_mr_idx, tmp_exp_idx


def do_utest_and_regression4dmr(mean_methy_data, len_of_mr_idx, exp_data,tmp_exp_idx, 
       isPlot=False, out_fig='fig_mr_vs_exp.jpg',xlabels='DNA Methylation', ylabels='Gene Expression', 
       gene_exp_cols=[],fig_dpi=100 ):
   '''do test or regression by using methylation level against the expression profiles across all samples'''
   #sort expression data based on mean of methyation in MR, then divide them into half 
   #represents regulation in unmethylated and methylated groups, respectively,
   #then do one-side utest for gene expresssions between u and m groups
   sorted_exp_data=[x for _,x in sorted(zip(mean_methy_data ,exp_data )) ]
   u_group=sorted_exp_data[0:len(sorted_exp_data)//2]
   m_group=sorted_exp_data[len(sorted_exp_data)//2:]
   tmp_utest=stats.mannwhitneyu(m_group,u_group, alternative='less')

   #regression analysis
   #X=mean_methy_data
   X=np.concatenate(( np.ones([len_of_mr_idx,1]), mean_methy_data.reshape([len_of_mr_idx,1]) ),axis=1)
   Y=exp_data.reshape([len(tmp_exp_idx),1])
   #reg = LinearRegression().fit(X, y)
   mod = sm.OLS(Y,X)
   fii = mod.fit()
   #pvalue to regression analysis
   p_values = fii.summary2().tables[1][['t','P>|t|']]
   reg_const_x1=fii.summary2().tables[1][['Coef.']].to_numpy()

   #plot fit
   if isPlot: 
      plt.clf()
      # p_values.loc['x1','P>|t|']<0.05 :
      plt.plot(X[:,1],Y,'o',color='red')
      if len(gene_exp_cols)>0:
          tmp_x=X[:,1]
          #add color for gcb or normal samples
          plt.plot(tmp_x[gene_exp_cols>=0],Y[gene_exp_cols>=0], 'o',color='green')
          plt.legend(['Tumor','Normal'])
 
      plt.plot(X[:,1],reg_const_x1[0]+ reg_const_x1[1]*X[:,1],'b-')
      plt.xlabel(xlabels, fontweight='bold', fontsize=10)
      plt.ylabel(ylabels, fontweight='bold',fontsize=10)
      plt.title( 'T-value =' +  '{:{width}.{prec}f}'.format(p_values['t'].x1,width=5,prec=1) + \
                 ', P-value <' + '{:{width}.{prec}e}'.format(p_values['P>|t|'].x1,width=5,prec=1), fontweight='bold',fontsize=12 )
      plt.savefig(out_fig,dpi=fig_dpi)
      print(out_fig)
   #  input('Click ')
   return tmp_utest, fii, p_values


def do_parallel4random_samples(args):
  '''do parallel computation for regression analysis in all genes '''
  (randomly_selected_mrs, p_values, in_geneexp_df, exp2mr_normal_id,
         mr_start_col, tmp_exp_idx,tmp_exp_sample_cols) =args
  rd_ls_p=0
  for ri in randomly_selected_mrs:
      rd_mr_path=ri
      rd_out_folder, rd_mr_file =os.path.split(rd_mr_path)
      rd_mean_methy_data, rd_exp_data, rd_not_find_files, rd_mr_idx,rd_exp_idx=combine_mr_and_expression_data(rd_mr_file,rd_out_folder ,
                                  in_geneexp_df, exp2mr_normal_id, tmp_exp_sample_cols, 
                                  mr_start_col)
      rd_len_of_mr_idx=len(rd_mr_idx)
      rd_utest, rd_fii, rd_p_values= do_utest_and_regression4dmr(rd_mean_methy_data, rd_len_of_mr_idx, rd_exp_data,tmp_exp_idx)
      if  rd_p_values.loc['x1','P>|t|']< p_values.loc['x1','P>|t|'] and p_values.loc['x1','t']*rd_p_values.loc['x1','t'] >0:
               rd_ls_p +=1
  return rd_ls_p


def run_parallel(args):
  ''' split data into blocks then do parallel computation for all of them'''
  (record_gene2mr, in_geneexp_df,in_mrs_not_in_enhancer_tss_df, 
      mr_in_block_folder,total_samples,num_of_processes,  
      exp2mr_normal_id, tmp_exp_sample_cols,  cutoff,reg_cutoff, 
      gene_exp_start_col,tmp_gene_col_name, not_dmr_folder_format,output_file_path)= args
  all_ttest={}
  not_find_files=[]
  gene_find_data=0
  not_find_geneexp=[]
  for ki in record_gene2mr.keys():
    #record_gene2mr is a dictionary , each key contains : mr_id, enhancer, gene_type, mutation block
    #where only mr_id, enhancer and gene_type will be used 
    (tmp_mr,tmp_enhancers,tmp_genetype,tmp_block_id)=record_gene2mr[ki]
    #find expression data
    print(ki)
    tmp_exp=in_geneexp_df[in_geneexp_df[tmp_gene_col_name]==ki]
    tmp_geneexp0=tmp_exp.iloc[:,gene_exp_start_col:].to_numpy()
    if tmp_exp.shape[0]==0:
       print('Gene expression not find for: ', ki, ' use zero instead!' )
       tmp_geneexp0=np.zeros((1,tmp_exp.shape[1]-gene_exp_start_col))
       not_find_geneexp.append(ki)

    tmp_gene_find_data=0
    #find methylation data
    tmp_list=[]
    len_of_tmp=len(tmp_mr)
    #loop in gene exp
    for ti in range(0,tmp_geneexp0.shape[0]):
     tmp_geneexp=tmp_geneexp0[ti]
     #loop in mr
     for i in range(0,len_of_tmp):
       ti=tmp_mr[i]
       if len(tmp_enhancers)==len_of_tmp:
         ei=tmp_enhancers[i]
       else:
         ei=tmp_enhancers[0]
       if len(tmp_genetype)==len_of_tmp:
         gi=tmp_genetype[i]
       else:
         gi=tmp_genetype[0]

       tmps0=ti.split(',')
       #subset of mrs
       for tii in tmps0:
         tmps=tii.strip().split(':')
         tmp_chr=tmps[0]
         tmp_mr_id=tmps[1]
         tmp_out_folder=mr_in_block_folder
         #fine gene for its expression and the methylation levels
         mr_start_col=1
         tmp_not_find_files=[]
   
         #added jbw 2023
         #here is the new export data formt from dmr_analysis
         if not not_dmr_folder_format:
           tmp_out_folder2= os.path.join(tmp_out_folder, tmp_chr,'data')
           #below is the old export data format by dmr_analysis
         else: 
           tmp_out_folder2=tmp_out_folder+tmp_chr
           print('MR data folder is not in dmr_analysis export format: ', tmp_out_folder2 )
         mean_methy_data, exp_data, tmp_not_find_files, tmp_mr_idx,tmp_exp_idx=combine_mr_and_expression_data('_'.join(tmps)+'_raw.dat.gz',
                                  tmp_out_folder2 , tmp_geneexp, exp2mr_normal_id, tmp_exp_sample_cols,
                                  mr_start_col   )
         #this is for debug purpose
         #print('Find methy data:', len(tmp_mr_idx))
         #np_exp_cols=np.array(tmp_exp_sample_cols)
         #print('Find expression data:', len(tmp_exp_idx), np_exp_cols[tmp_exp_idx])
         if len(tmp_not_find_files)==0:
           len_of_mr_idx=len(tmp_mr_idx)
           #do regression and uttest, here p_values are regressio analysis pval
           tmp_utest, fii, p_values= do_utest_and_regression4dmr(mean_methy_data, len_of_mr_idx, exp_data,tmp_exp_idx)

           #do random sampling for 1000 mrs genome widely
           #total_samples=1000
           randomly_selected_mrs=in_mrs_not_in_enhancer_tss_df.sample(total_samples,random_state=1)
           #num_of_processes=25
           block_chunks, num_of_processes=divide_data2chunk(num_of_processes, randomly_selected_mrs[0].tolist())

           pool = mp.Pool(processes=num_of_processes)
           rd_ls_ps= pool.map(do_parallel4random_samples,
             [ (block_chunks[loop], p_values,tmp_geneexp, exp2mr_normal_id, mr_start_col,tmp_exp_idx, tmp_exp_sample_cols) for loop in range(0, num_of_processes) ], chunksize=1)
           rd_ls_p=sum(rd_ls_ps)
           pool.close()

           ex_p=rd_ls_p/total_samples
           tmp_gene_find_data +=1

           if p_values.loc['x1','P>|t|']<reg_cutoff and  ex_p< cutoff : #and gi in ['TSS','5dist']:
             tmp_list.append((tii,p_values.loc['x1','t'],p_values.loc['x1','P>|t|'],ei,gi,tmp_utest.pvalue, ex_p))

           #record not find files
         else:
            #print(tmp_not_find_files)
            not_find_files += tmp_not_find_files

    #record gene with mrs passed filtering condition
    if len(tmp_list)>0:
      all_ttest[ki]=tmp_list
      #export temparty file
      json.dump( all_ttest, open( os.path.join(output_file_path, "_all_ttest_from_methy_vs_express.json"), 'w' ) )
 
    if tmp_gene_find_data>0:
      gene_find_data +=1

  print('Genes find data: ', gene_find_data )
  print('Genes find data and passed filtering: ', len(all_ttest))

  return all_ttest, not_find_files

def export_results(all_ttest, not_find_files,output_file,output_file_path):
  '''Export results from regresssion analysis '''
  import json
  out_file1=os.path.join(output_file_path, output_file+'_not_find_mr.json')
  print('Export at: ', out_file1 )
  #remove duplicated results
  uq_not_find_files=list(set(not_find_files))
  with open(out_file1,'w') as fp:
     json.dump(uq_not_find_files,fp)

  #fp= open('test_data.json') 
  #all_ttest= json.load(fp)
  all_ttest_df=pd.DataFrame(data=all_ttest.items())
  new_gene=[]
  new_row=[]
  for idx, row in all_ttest_df.iterrows():
       tmp_gene=row[0]
       tmp_row=row[1]
       for tr in tmp_row:
          new_gene.append(tmp_gene)
          new_row.append(tr)

  out_file2=os.path.join(output_file_path, output_file+'.csv')
  print('Export: ', out_file2)
  new_df1=pd.DataFrame(data=new_row,columns=['mr_id','rg_Tvalue','rg_Pvalue','enhancer','gene_type','utest_Pvalue','expectPvalue'])
  new_df1.insert(0, 'gene',new_gene)
  #out_file='test_result_of_1000sampling.csv'
  new_df1.to_csv(out_file2, sep='\t',index=False)
  return out_file1, out_file2

def run(args):
  '''Input parameters then run the calculations'''
  #cutoff=0.05
  cutoff=args.expected_pval_cutoff
  reg_cutoff=args.regression_pval_cutoff
  print('P-value cutoff for expected P: ', str(cutoff))
  print('P-value cutoff for regression model when fitting expression data to methylation data:', str(reg_cutoff))
  output_file_name=args.output_file_name
  output_file_path=args.output_files_path
  not_dmr_folder_format=args.not_dmr_analysis_format 
  
  #new results in novembe 2021 from dds_analysiss check_block_gene_inTAD: gene_name; gene_type; block_id; new_mr_sites; patients; isTAD; enhancers; patient_id 
  #in_gene_mr_file='out_blocks_dmr_single3_nov2021/out_dds/blocks_summary_block_position_0flank_0.7Proba_66868blocks_13143blocks2mr_4603blocks2dmr_deg_info_filtered_DMR_or_DEG_uniqGene_commonTAD_Boundary_list2UqGene.tsv'
  #in_gene_mr_file='dds_out.tsv'
  in_gene_mr_file=args.in_gene_mr_file_from_dds
  mr_id_column_name=args.mr_id_column_name
  print('Column name for MR IDs: ', mr_id_column_name)
  in_df=pd.read_csv(in_gene_mr_file, sep='\t')
  #keep only unique mr id in each row
  in_df[mr_id_column_name]=in_df[mr_id_column_name].apply(lambda x: '~'.join(list(set(x.split('~')))))

  if not args.is_tab_delimated_file4MR :
    #in_df=pd.read_csv(in_gene_mr_file,sep='\t')
    print('Input file of MR/block regions from dds_analysis exported file: ',  in_gene_mr_file)

    #count how many genes have DMRs
    mr_df=in_df[in_df[mr_id_column_name].apply(lambda x:   str(x).find(':mr') >=0)].copy()
    mr_df.reset_index(inplace=True)
    #mr_df[mr_df.gene_name.apply(lambda x: x.find('EZH2')>=0)]

    #check how many mr located in TSS or 5distance regions
    #out is a dictionary with gene as key, and 4 element hash : new_mr_sites, enhancers, gene_type, block_id
    record_gene2mr,out_mr_in_blocks=get_mr_in_blocks(mr_df)
  else:
    print('Input file of MR/block regions from a common tab delimated file: ', in_gene_mr_file) 
    #for a common tab file , it must have gene_name, gene_type, new_mr_sites, and enhancers columns
    record_gene2mr={}
    in_df.fillna('nan',inplace=True)
    if 'enhancers' in in_df.columns:
      for index,row in in_df.iterrows():
        if row.gene_name in record_gene2mr.keys():
           tmp1,tmp2,tmp3,tmp4=record_gene2mr[row.gene_name]
           record_gene2mr[row.gene_name]=(tmp1+row[mr_id_column_name].split('~') , tmp2+row.enhancers.split('~'), tmp3+row.gene_type.split('~'), tmp4+[''])
        else:
           record_gene2mr[row.gene_name]=(row[mr_id_column_name].split('~'),row.enhancers.split('~'),row.gene_type.split('~'),[''])
    else:
      for index,row in in_df.iterrows():
        if row.gene_name in record_gene2mr.keys():
           tmp1,tmp2,tmp3,tmp4=record_gene2mr[row.gene_name]
           record_gene2mr[row.gene_name]=(tmp1+row[mr_id_column_name].split('~') , tmp2+['nan'], tmp3+row.gene_type.split('~'), tmp4+[''])
        else:
           record_gene2mr[row.gene_name]=(row[mr_id_column_name].split('~'),['nan'],row.gene_type.split('~'),[''])
 

  #gene expression data
  #new results in november 2021 from bpb3 
  #in_gene_expression='new_diffExp_genes/differentially_expressed_genes.txt_12Kdiffgenes'
  #in_gene_expression='bpb3_diffExp.tsv'
  in_gene_expression=args.in_gene_expression_file_from_bpb3
  #gene_exp_start_col=2
  gene_exp_start_col=args.gene_expression_start_column
  pathDepth4sampleName=args.in_diffExp_pathDepth4sampleName

  if args.is_tab_delimated_file4diffExp:
    print('Input of gene expression file is a common tab delimiated file:', in_gene_expression)
    in_geneexp_df=pd.read_csv(in_gene_expression,sep='\t')
    tmp_exp_cols=in_geneexp_df.columns[1:].to_list()
    tmp_exp_sample_cols=tmp_exp_cols
    tmp_gene_col_name=in_geneexp_df.columns[0]
    gene_exp_start_col=1
  else:
    in_geneexp_df, tmp_exp_cols, tmp_exp_sample_cols=read_geneexp_data(in_gene_expression, gene_exp_start_col,pathDepth4sampleName)
    tmp_gene_col_name=in_geneexp_df.columns[0]
    print('Input file of bpb3 exported differentially expressed genes: ', in_gene_expression)

  #mr list for sampling in background model
  #mrs_not_in_enhancer_tss_file='mr4sampling_list.csv'
  mrs_not_in_enhancer_tss_file=args.in_mr_background_file_list
  print('Input file of background regions: ', mrs_not_in_enhancer_tss_file)

  in_mrs_not_in_enhancer_tss_df=pd.read_csv(mrs_not_in_enhancer_tss_file, sep='\t',header=None)

  #perform gene expression versus methyltion test for blocked mr and 10000 randomly selected mr from genome
  #here for 4 normal samples just simply assign normal mr sample id to the normal exp
  #exp2mr_normal_id={'SRR834983': 'gcb_4118819', 'SRR834984':  'gcb_4122131' , 'SRR834985': 'gcb_4160735', 'SRR834986':'gcb_4174884' }
 
  in_sampleName=args.in_sampleName_for_replace
  if os.path.exists(in_sampleName)  and os.stat(in_sampleName).st_size>0 :
     in_sample_df=pd.read_csv(in_sampleName, sep='\t',index_col=0, header=None)
     exp2mr_normal_id=in_sample_df.to_dict()[1]
     print('Input file of sample name replacement: ', in_sampleName)
  else:
     print('Input file of sample name replacement is not available or empty, skip it!', in_sampleName)
     exp2mr_normal_id=[]


  #methylation data in each region: extracted mr data for testing in the regions/block
  #mr_in_block_folder= 'out_mr_in_block_november/mr_in_block_out'
  mr_in_block_folder=args.in_methylation_folder_from_dmr
  print('Input folder of methylation data: ', mr_in_block_folder)

  #total_samples=1000
  #num_of_processes=15
  
  total_samples=args.total_samples_in_random_selection
  num_of_processes=args.number_of_parallel_processes
  print('Total number of random samples draw from the background: ', str(total_samples))
  print('Number of parallel processes: ', str(num_of_processes))

  all_ttest, not_find_files=run_parallel((record_gene2mr, in_geneexp_df,in_mrs_not_in_enhancer_tss_df,
                                mr_in_block_folder,total_samples,num_of_processes,exp2mr_normal_id, tmp_exp_sample_cols, 
                                cutoff,reg_cutoff, gene_exp_start_col, tmp_gene_col_name, not_dmr_folder_format,output_file_path))
  export_results(all_ttest, not_find_files, output_file_name+str(total_samples)+'sampling',output_file_path)


if __name__=='__main__':
  args=my_parser(argparse.ArgumentParser('python dTarget_methy_vs_express.py')).parse_args()
  run(args)

