#this script is used to plot and do linear regression between mr and expression
import pandas as pd
import os
import glob
import numpy as np
from dds_analysis.script.dTarget_methy_vs_express import read_geneexp_data, combine_mr_and_expression_data, do_utest_and_regression4dmr
from dds_analysis.script.plot_tss_enhancer_mrs import check_folder
import argparse
# exec(open('plot_mr_vs_exp.py').read())

def my_parser(parser):
  required= parser.add_argument_group('Required')
  required.add_argument('-inGeneEXPfile','--in_gene_expression_file_from_bpb3', required=True, type=str, \
                         help='Input file name for gene expression profiles that exported by bpb3 differential_expression,\
                                             file format: gene_name; P-values, sample1, sample2 , ....  ')
  required.add_argument('-inMRfolder','--in_methylation_folder_from_dmr', required=True, type=str, \
                         help='Input file folder path for methylation data, the data are extracted by \
                                              dmr_analysis dmr_exportData based on predifined blocks or regions')
  optional= parser.add_argument_group('Optional, with default values ')
  optional.add_argument('-inGene','--in_gene_name', type=str, help='target gene name for a MR/DMR region, default=None',default=None, metavar='')
  optional.add_argument('-inMR','--in_mr_id', type=str, help='MR/DMR ID of a target gene such as chr1:mr1234, default=None', default=None,metavar='')
  optional.add_argument('-inFile', '--in_gene_mr_file', type=str, help="Input file contains both DMR/MR ID and the target gene name, which is a tab \
                         delimated texf file, and the column lablels are mr_id and gene, respectively. The file format is the same as export file from dTarget_methy_vs_express.py default=None", default=None, metavar='' )

  optional.add_argument('-expTAB','--is_tab_delimated_file4diffExp', \
                         help="Whether input differential expression gene file is a bpb3 differential_expression exported file format or a common tab delimated file format \
  (e.g., the first column is the gene name and the rest of columns is gene expression data in samples) , default=False, if use this optional then input gene expression file \
                                                            will be treated as common tab delimated file.", action="store_true" )
  optional.add_argument('-expStartCol','--gene_expression_start_column', default=2, metavar='',type=int, \
                 help='Start column number (e.g., 0,1,2,...) for samples in gene expression profile file, \
                  default =2 for bpb3 exported file')
  optional.add_argument('-pathDepth','--in_diffExp_pathDepth4sampleName',default=12, type=int, \
                 help='file path depth for column sample name, for example, if columna name is /path/folder/sample_name/.., then the \
                       depth of sample name is 3, in other words, it is the number of "/" before the sample name, default=12')
  optional.add_argument('-sampleName','--in_sampleName_for_replace', default='sample_name4replace.tsv', type=str, \
                 help='Input file name for rename the old sample names in gene expression data, first column old_name is old sample name, second \
   column new_name is the new name that used for replacement, default file name is sample_name4replace.tsv which can be an empty file if no replacement is needed for sample names.')
  optional.add_argument('-dpi', '--figure_resolution_dpi', type=int, default=60, metavar='', help='Exported figure resolution in dpi, default dpi=60')
  optional.add_argument('-wtStr','--wildType_fileString', default='gcb_', type=str, metavar='',
                 help='First few character string that labeled in file name as Wide Type condition, for example, if file name start with gcb_meht1_* as wild type control sample'
                              ', then --wildType_fileString is gcb, which is the default setting in the program ') 
  optional.add_argument('-output_path','--output_file_path', default='./', metavar='',
                 help='Output path for exported figures, default = ./ at local directory.', type=str)
  return parser

def replace_sample_name(exp2mr_normal_id, tmp_geneexp_cols2):
   '''replace sample name based a dictioanry from sample_name_replace {old_name: new_name}
   '''
   new_exp_cols=[]
   if len(exp2mr_normal_id)>0:
      for tmi in tmp_geneexp_cols2:
           if tmi in exp2mr_normal_id.keys() :
             new_exp_cols.append(exp2mr_normal_id[tmi])
           else:
             new_exp_cols.append(tmi)
   else:
      new_exp_cols=tmp_exp_cols
   return new_exp_cols

def plot_mr_vs_expression_figure(in_gene_name, in_mr_name, in_geneexp_df,
     gene_exp_start_col,exp2mr_normal_id,tmp_exp_sample_cols, args):
  '''plot mr vs. expression figure,
  '''
  #find the gene needs to be plot
  tmp_gene_col_name=in_geneexp_df.columns[0]
  tmp_exp=in_geneexp_df[in_geneexp_df[tmp_gene_col_name]==in_gene_name]
  tmp_geneexp=tmp_exp.iloc[:,gene_exp_start_col:].to_numpy()
  tmp_geneexp_cols=tmp_exp.columns[gene_exp_start_col:]
  if tmp_exp.shape[0]==0:
     print('Gene expression not find for: ', in_gene_name, ' use zero instead!' )
     tmp_geneexp=np.zeros((1,tmp_exp.shape[1]-gene_exp_start_col))

  #split mr name to chr and mr id
  tmp_mr=in_mr_name.split(':')
  tmp_chr=tmp_mr[0]
  tmp_mr_id=tmp_mr[1]
  mr_start_col=1
  mr_in_block_folder=args.in_methylation_folder_from_dmr
  print('Input methylation data: ', mr_in_block_folder)

  #combine mr data and gene expression data
  #print('_'.join([tmp_chr,tmp_mr_id,'raw.dat.gz']))

  mean_methy_data, exp_data, tmp_not_find_files, tmp_mr_idx,tmp_exp_idx=combine_mr_and_expression_data('_'.join([tmp_chr,tmp_mr_id,'raw.dat.gz']),
                                  os.path.join(mr_in_block_folder,tmp_chr,'data') , tmp_geneexp, exp2mr_normal_id, tmp_exp_sample_cols,
                                  mr_start_col )

  #do regression and plot figure
  len_of_mr_idx=len(tmp_mr_idx)
  tmp_geneexp_cols1=tmp_geneexp_cols[tmp_exp_idx].to_list()
  #print(tmp_exp_idx,tmp_geneexp_cols1, args.wildType_fileString)
  #print(exp2mr_normal_id) 

  #replace gene express sample name by using inputted sampeName_replacement
  tmp_geneexp_cols2=replace_sample_name(exp2mr_normal_id, tmp_geneexp_cols1)
  #print(tmp_geneexp_cols2)
  #jbw may 
  if len(tmp_geneexp_cols2)>0:
    tmp_geneexp_cols_label=np.char.find(tmp_geneexp_cols2,args.wildType_fileString)
    tmp_utest, fii, p_value= do_utest_and_regression4dmr(mean_methy_data, len_of_mr_idx, exp_data,tmp_exp_idx
        ,True, os.path.join(args.output_file_path, in_gene_name+'_'+in_mr_name.replace(':','_') + '.jpg'), 
         in_mr_name+' (Methylation)',in_gene_name+' (Expression)',tmp_geneexp_cols_label,args.figure_resolution_dpi)
  else:
    tmp_utest=[]
    fii=[]
    p_value=[]
    print('Not find ',in_gene_name)
    print('No mr vs exp plot!')
  return tmp_utest, fii, p_value

def run(args):
  #input gene files
  #in_gene_exp_file='data/fl_14samples_dmr_april22/new_data/differentially_expressed_genes_tab.tsv'
  in_gene_exp_file=args.in_gene_expression_file_from_bpb3

  #input gene name and mr id
  #in_gene_name='IGHJ4'
  #in_mr_name='chr14:mr34378:hyper:D' #'chr18:mr22253:hypo:D' #'chr14:mr34378:hyper:D' #'chr14:mr34361:hyper:D'
  in_gene_name=args.in_gene_name
  in_mr_name=args.in_mr_id
  in_gene_mr_file=args.in_gene_mr_file

  #jbw may 
  check_folder(args.output_file_path)

  #read gene expression data
  in_geneexp_df=pd.read_csv(in_gene_exp_file,sep='\t')
  tmp_gene_col_name=in_geneexp_df.columns[0]
  gene_exp_start_col=args.gene_expression_start_column
  pathDepth4sampleName=args.in_diffExp_pathDepth4sampleName
  in_sampleName=args.in_sampleName_for_replace
  if args.is_tab_delimated_file4diffExp:
    print('Input of gene expression file is a common tab delimiated file:', in_gene_exp_file)
    in_geneexp_df=pd.read_csv(in_gene_exp_file,sep='\t')
    tmp_exp_cols=in_geneexp_df.columns[1:].to_list()
    tmp_exp_sample_cols=tmp_exp_cols
    tmp_gene_col_name=in_geneexp_df.columns[0]
    gene_exp_start_col=1
  else:
    in_geneexp_df, tmp_exp_cols, tmp_exp_sample_cols=read_geneexp_data(in_gene_exp_file, gene_exp_start_col,pathDepth4sampleName)
    tmp_gene_col_name=in_geneexp_df.columns[0]
    print('Input file of bpb3 exported differentially expressed genes: ', in_gene_expression)

  ##find sample name replacement file
  in_sampleName=args.in_sampleName_for_replace
  if os.path.exists(in_sampleName) and os.stat(in_sampleName).st_size>0 :
     in_sample_df=pd.read_csv(in_sampleName, sep='\t',index_col=0, header=None)
     exp2mr_normal_id=in_sample_df.to_dict()[1]
     print('Input file of sample name replacement: ', in_sampleName)
  else:
     print('Input file of sample name replacement is not available or empty, skip it !', in_sampleName)
     exp2mr_normal_id=[]

  if in_gene_name==None or in_mr_name==None:
     #try read input gene and mr id from a tab delimiated file
      if in_gene_mr_file != None:
         print('Read DMR/MR ID and gene name from a file:', in_gene_mr_file)
         gene_mr_df=pd.read_csv(in_gene_mr_file,sep='\t')
         for index, rows in gene_mr_df.iterrows():
             in_gene_name=rows.gene
             if 'dmr_id' in gene_mr_df.columns.to_list():
               in_mr_name=rows.dmr_id 
             else:
               in_mr_name=rows.mr_id
             print(in_gene_name , ' -> ', in_mr_name)
             plot_mr_vs_expression_figure(in_gene_name, in_mr_name, in_geneexp_df,gene_exp_start_col, exp2mr_normal_id, tmp_exp_sample_cols, args)
      else:
         print('No input information for Gene name and DMR/MR ID , Please try again! ')
         exit() 
  else:
     #gene and mr id are available from input option
     plot_mr_vs_expression_figure(in_gene_name, in_mr_name, in_geneexp_df,gene_exp_start_col, exp2mr_normal_id, tmp_exp_sample_cols, args)

if __name__=='__main__':
  args=my_parser(argparse.ArgumentParser('python plot_mr_vs_exp.py')).parse_args()
  run(args)



