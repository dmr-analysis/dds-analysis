#this script is used to plot the average methylatoin levels in promoter and enhancer for a selceted enhancer target gene!
#which is predicted by dTarget_methy_vs_express
import pandas as pd
import os
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import scipy.interpolate as sc_i
import numpy as np
from .script_high.plot_TSSdistanceRegion import main as plt_tgt_main
import argparse
#from hmst_seq_analyzer.scripts_high.scripts.plot_TSSgeneTES import main as plt_tgt_main
#exec(open('plot_tss_enhancer_mrs.py').read())

def my_parser(parser):
  required= parser.add_argument_group('Required')
  #required.add_argument('-folder','--in_data_folder', required=True, type=str , \
  #                      help='Input data folder where conains files for DEG, DMR, predicted tss_file ,enc_file ')
  required.add_argument('-exp_file','--differential_gene_expression_file', required=True, type= str, \
                        help='A tab delimited file with at least four columns (e.g., gene, group1_mean, group2_mean, rratio/log2folder change) '\
                        'contains differential gene expressions information, which is exported by (bpb3 filterDEG4bpb3). If more than four columns are ' \
                        'included in this file , then the program will assume the first column is gene, and the last 3 columns are group1_mean, '\
                        'group2_mean, and rratio/log2folder change, respectively. ')
  required.add_argument('-dmr_file','--differential_methylation_file',type=str, \
                        help='A tab delimited file with only five columns (e.g., chrorm, start_pos, end_pos, infos, pval) contains differential methylation '\
                        'regions which is exported by (dmr_analysis  dmr_combine_multChrs4rank). For the column infos, the data structure is chrom:mr_id:anything:hyper/hypo:anything, '\
                        'where the information in the first, second,and forth strings separted by : will be used by the program for finding DMRs . ')
  required.add_argument('-tss_file','--tss_target_file',type=str, \
                        help='A tab delimited file with nine columns (index, gene, mr_id, rg_Tvalue, rg_Pvalue, enhancer, gene_type, utest_Pvalue, expectPvalue ) conains ' \
                        ' predicted information of TSS target genes by using ( dds_analyssís dTarget_methy_vs_express)')
  required.add_argument('-enc_file','--enhancer_target_file',type=str, \
                         help='A tab delimited file with nine columns (index, gene, mr_id, rg_Tvalue, rg_Pvalue, enhancer, gene_type, utest_Pvalue, expectPvalue ) contains ' \
                        'predicted information of enhancer target genes by using ( dds_analyssís dTarget_methy_vs_express)' ) 
  required.add_argument('-genes','--selected_genes4plot', required=True, type=str, \
                        help='a string of gene symbols that will be selected for ploting the average methylation levles in their predicted regulatory ' \
                        ' regions such as TSS and Enhancer. Gene symbols are separeted by comma(,) such as gene1,gene2,gene3 .' )
  required.add_argument('-mr_folder','--mr_data_folder',type=str, \
                        help='A file folder contains all DNA methylation data of methylated regions that exported by dmr_analysis dmr_exportData')
  required.add_argument('-folder_name','--prefix_foldername4chrom',type=str, \
                        help='A prefix folder name for exported methylation data of each chromosome that under -mr_folder ' \
                        '(e.g., prefix_foldername4chrom + chr1 means all methylation data for chr1)')

  optional= parser.add_argument_group('Optional, with default values')
  optional.add_argument('-out_folder','--out_data_folder',type=str, default='data', metavar='', \
                        help='Out folder name for exported data from the current analysis, if this folder does ' \
                              'not exist in in_data_folder then it will be created. Default=data')
  optional.add_argument('-out_plots','--out_plot_folder', type=str, \
                        help="Out folder name for all plots, default=plots ", default= 'plots',metavar='')
  optional.add_argument('-is_negative','--is_negative_target', type=int, \
                        help='Select plot for negative, positive or both responses of predicted DMR targets, 0,1,2, for positive, negative, and all responses , respectively. Default=1 for negative targets',
                          default=1, metavar='')
  optional.add_argument('-gX','--gene_upstream_X', type=int , \
                        help='set lenght of upstream to transcription start side - TSS in plot, default=500', default=500 , metavar='') 
  optional.add_argument('-gY','--gene_downstream_Y', type=int, \
                        help='set length of downstream to transcripton start side -TSS in plot, default=500', default=500, metavar='')
  optional.add_argument('-G','--gene_or_enhancer_length',type=int, \
                        help='set gene or enhancer length in plot, default= 2000 ', default=2000, metavar='' )
  optional.add_argument('-w','--step_size_in_region', type=int, \
                        help='Step size in a region for data smoothing, default=100', default=100, metavar='' )
  optional.add_argument('-window','--window_size_in_combined_regions', type=int, \
                        help='Window size of combined regions such as TSS and Enhancer, default= 600. ' \
                              'To increase this parameter >500 if there is a gap between two combined regions', default=600, metavar='')
  optional.add_argument('-sigma','--data_smooth_sigma', type=int, \
                        help='data smooth paramter sigma, default=50', default=50, metavar='')
  optional.add_argument('-flank','--flank_region', type=int, \
                        help='Flank region to add on the two-sides of a promoter or enhancer region , default=100', default=100, metavar='')
  optional.add_argument('-MRs_txt','--MRs_txt_description', type=str, \
                        help='Text description of MRs. Default=DMRs',default='DMRs', metavar='') 
  optional.add_argument('-dmr_compressed','--dmr_file_not_compressed', action="store_false",
                        help='If dmr file is not compressed file then use this option, default=False, dmr_file is compressed !')
  optional.add_argument('-flip_rr','--not_flip_sign_of_rratio', action='store_false',
                        help='If not flip sign of rratio then use this option, default=False, sign of rratio is flipped!')
  optional.add_argument('-wtStr','--wildType_fileString', default='gcb_',type=str, metavar='',
                        help='First few character string or strings that labeled in file name as Wide Type condition, for example, if file name start with gcb_meht1_* as wild type control sample'
                              ', then --wildType_fileString is gcb, which is the default setting in the program ')

  return parser


def unique_list(in_list):
  ''' remove duplicates in a list and keep 
      all elements in the same order as the original list
  '''
  new_list=[]
  for i in in_list:
      if i not in new_list:
          new_list.append(i)

  return new_list


def check_folder(out_folder):
 if not os.path.exists(out_folder):
    print("Create , ", out_folder)
    os.makedirs(out_folder)
 else:
    print("Exists , ", out_folder)

def find_deg_dmr4target_genes(out_folder, tss_file, enc_file,exp_file,dmr_file, is_dmr_compressed =True, is_flip_sign_of_rr_ratio=True):
  ''' find DEG and DMR information for predicted target genes,
     folder: file folder contains input data of tss_file, enc_file and exp_file
     exp_file: DEG file from bpb3 by rr-ratio filtered genes:  bpb3 filterDEG4bpb3 --help
     dmr_file: DMR file from dmr_analysis: dmr_analysis  dmr_combine_multChrs4rank 
     tss_file and enc_file: predicted tss, enhancer target genes from DDS: dds_analyssís dTarget_methy_vs_express

  '''
  #read data
  #jbw read data new, no index in input file
  tss_df=pd.read_csv( tss_file,sep='\t') #,index_col=0)
  enhancer_df=pd.read_csv( enc_file ,sep='\t') #,index_col=0)
  exp_df=pd.read_csv( exp_file ,sep='\t')
  if is_dmr_compressed:
    dmr_df=pd.read_csv(dmr_file ,sep='\t',header=None, compression='gzip')
  else:
    dmr_df=pd.read_csv( dmr_file ,sep='\t',header=None)
  dmr_df.columns=['chr','start_pos','end_pos','infos','pval']


  #remove duplicated rows from predicted tss/enhancer target genes !!
  tss_df.drop_duplicates(['gene','mr_id'],inplace=True)
  enhancer_df.drop_duplicates(['gene','mr_id'],inplace=True)

  #export DEG DMR information for predicted target genes
  #exp_data_df=exp_df[['gene','HAP1_KO1','HAP1_P1','rratio']].copy()
  #assume the dataframe the first column is gene id and the last three coumns are group means and rr-ratios that compluted by bpb3 filterDEG4bpb3
  exp_data_df=exp_df.iloc[:,[0,-2,-3,-1]].copy()
  exp_data_df.columns=['gene','group1_exp_mean','group2_exp_mean','rratio']
  #change sign of DEG comparison to KO vs. Wild
  if is_flip_sign_of_rr_ratio:
    exp_data_df.rratio=exp_data_df.rratio*(-1)
  #add DEG info to tss and enhancer 
  tss_df2= pd.merge(tss_df,exp_data_df,on='gene',how='left').copy()
  enhancer_df2= pd.merge(enhancer_df,exp_data_df,on='gene',how='left').copy()
  #add DMR info to tss and enhancer
  info_df=dmr_df.infos.str.split(':',expand=True).copy()
  dmr_df['mr_id']=info_df[0]+':'+info_df[1]
  dmr_df['dmr_KO_vs_P']=info_df[3]
  dmr_data_df=dmr_df[['mr_id','dmr_KO_vs_P']]
  tss_df3=pd.merge(tss_df2, dmr_data_df,on='mr_id',how='left').copy()
  enhancer_df3=pd.merge(enhancer_df2, dmr_data_df,on='mr_id',how='left').copy()

  #jbw export new file to out_folder
  out_tss_file=os.path.join(out_folder,os.path.basename( tss_file.replace('.csv','_deg_dmr_info.tsv')) )
  out_enc_file=os.path.join(out_folder, os.path.basename( enc_file.replace('.csv','_deg_dmr_info.tsv')))
  tss_df3.to_csv(out_tss_file, sep='\t',index=False)
  enhancer_df3.to_csv(out_enc_file,sep='\t',index=False)
  return out_tss_file, out_enc_file, tss_df3.copy(), enhancer_df3.copy()
 

def extract_data4selected_genes(out_folder, tss_file, enhancer_file,selected_genes):
  '''
    Extract data for selected genes for ploting their mean methylation levels in tss and enhancers.
    tss_file and enhancer_file are predicted target gene files from dds_analyssís dTarget_methy_vs_express  and added
    deg and dmr information by using find_deg_dmr4target_genes
    folder: file folder contains tss_file and enhancer_file
  '''
  #jbw 
  tss_df=pd.read_csv( tss_file ,sep='\t')
  enhancer_df=pd.read_csv( enhancer_file ,sep='\t')


  #remove duplicates in rows!!
  tss_df.drop_duplicates(['gene','mr_id'],inplace=True)
  enhancer_df.drop_duplicates(['gene','mr_id'],inplace=True)

  record_tss_df=pd.DataFrame()
  record_enhancer_df=pd.DataFrame()
  for gi in selected_genes:
    gi_tss_df=tss_df[tss_df.gene==gi].copy()
    gi_enc_df=enhancer_df[enhancer_df.gene==gi].copy()
    record_tss_df=pd.concat([record_tss_df,gi_tss_df],axis=0)
    record_enhancer_df=pd.concat([record_enhancer_df,gi_enc_df],axis=0)
  #export file
  out_tss_file=os.path.join(out_folder,'selected_genes_tss.tsv')
  out_enc_file=os.path.join(out_folder,'selected_genes_enhancer.tsv')
  #export
  record_tss_df.to_csv(out_tss_file,sep='\t',index=False)
  record_enhancer_df.to_csv(out_enc_file,sep='\t',index=False)
  return record_tss_df.copy(), record_enhancer_df.copy()

def generate_hmtseq_data(wtype_str, tmp_gene, flank, feature, enc_file,isNegative, mr_folder, folder_name,out_folder ):
  ''' Genreate hmst-seq-analyzer format DMRs files for making plots in tss and enhancer
      tmp_gene: gene_name/symbol to extract methylation data 
      flank: flank region to add on the two-side of promoter or enhancer regions
      isNegative: negative 1 , positive 0 , or all 2 response of predicted target gene with regulatory regions.
      feature: tss or enhancer 
      enc_file: enhancer file
      mr_folder: file folder contains all methylation data exported by dmr_analysis
      folder_name: prefix folder name in mr_folder that represent each chromosome
      out_folder: export file folder
  '''
  #feature='tss'
  #enc_file=os.path.join(folder, 'selected_genes_'+feature+'.tsv')
  #target gene with negative or positive correlation between expression and methylation levels
  #isNegative=False
  check_folder(out_folder)
  enhancer_df=pd.read_csv(enc_file,sep='\t')
  uniq_genes=enhancer_df.gene.unique()
  record_gene={}
  for gi in uniq_genes:
    if isNegative==1:
      record_gene[gi]=enhancer_df.mr_id[(enhancer_df.gene==gi) & (enhancer_df.rg_Tvalue<0) ].to_list()
      #out_string='neg'
    elif isNegative==0:
      record_gene[gi]=enhancer_df.mr_id[(enhancer_df.gene==gi) & (enhancer_df.rg_Tvalue>0) ].to_list()
      #out_string='pos'
    else:
      record_gene[gi]=enhancer_df.mr_id[(enhancer_df.gene==gi)  ].to_list()
      #out_string='all'

  if isNegative==1: 
      out_string='neg'
  elif isNegative==0:
      out_string='pos'
  else:
      out_string='all'

  #get selected genes' MR data and convert it to hmst-seq-analyzer foramt
  #exported mr data from dmr_analysis dmr_exportData  
  #mr_folder='mr_data/dmr_in_deg_tss_5dist/'
  #folder_name='out4dmr_in_deg_tss_5dist'
  ##export selected genes MR data
  #out_folder='demo/data'

  #tmp_gene='ULK1'
  tmp_meth_points_all=[]
  tmp_power=[]
  tmp_tissue=[]
  tmp_df=pd.DataFrame()
  #jbw
  if tmp_gene in record_gene.keys():
    tmp_gi=record_gene[tmp_gene]
    for mi in tmp_gi:
      chrom, mr_id=mi.split(':')
      tmp_data_file=os.path.join(mr_folder, folder_name+chrom,'data',chrom+'_'+mr_id+'_raw.dat.gz')
      tmp_df=pd.read_csv(tmp_data_file,sep='\t')
      tmp_tissue.append(mi)
      tmp_power.append(tmp_df.iloc[:,1:].to_numpy())
      tmp_meth_points_all.append(tmp_df.position.to_list())

  enhancers_mc=[]
  tmp_columns=tmp_df.columns[1:].to_list()
  #add 100bp flank regons at two-sides of MR
  #flank=100
  if tmp_df.shape[0] !=0:
    for i in range(0,len(tmp_columns)):
      #jbw may, find wtype str index
      #if wtype_str in tmp_columns[i]:
      #   print(tmp_columns[i].split(wtype_str))
      #   tmp_file_str=[wtype_str + tmp_columns[i].split(wtype_str)[1]]
      #else:
      wtype_str_index=np.where(np.array(tmp_columns[i].split('_'))==wtype_str)
      if len(wtype_str_index[0])>0:
           start_idx=wtype_str_index[0][0]
      else:
           start_idx=1
      end_idx=start_idx+3
      tmp_file_str=unique_list(tmp_columns[i].split('_')[start_idx:end_idx])
      #jbw
      #print(tmp_file_str)
      out_file='_'.join(tmp_file_str) +'.'+tmp_gene+ '_'+out_string+'_5mC_'+feature+'_allDMRs.csv'
      out_file=os.path.join(out_folder,out_file)
      enhancers_mc.append(out_file)
      out_dmr_df=pd.DataFrame(columns=['region_start','region_end','meth_points_all','powers','tissue'])
      convert_dict = {'meth_points_all': object,
                'powers': object
                }
      out_dmr_df=out_dmr_df.astype(convert_dict)
      #out_dmr_df['power']=out_dmr_df['power'].astype(object)
      #out_dmr_df['meth_points_all']=out_dmr_df['meth_points_all'].astype(object)
      #tmp_regions= tmp_region.start_pos.to_list()+ tmp_region.end_pos.to_list()
      for j in range(0, len(tmp_tissue)):
        tmp_line=[]
        #tmp_regions= tmp_region.start_pos.to_list()+ tmp_region.end_pos.to_list()
        tmp_meth_pa=tmp_meth_points_all[j]
        tmp_regions=[min(tmp_meth_pa)-flank,max(tmp_meth_pa)+flank]
        tmp_pw=list(np.array(tmp_power[j][:,i]/100,dtype=object))
        tmp_line=tmp_regions
        tmp_line.append(tmp_meth_pa)
        tmp_line.append(tmp_pw)
        tmp_line.append(tmp_tissue[j])
        out_dmr_df.loc[len(out_dmr_df)]= tmp_line
      out_dmr_df.to_csv(out_file,index=False)
  else:
    print('Not find', tmp_gene)
  
  out_file_df=pd.DataFrame(data=enhancers_mc)
  out_file=os.path.join(out_folder, 'list_of_'+tmp_gene+'_'+out_string+ '_'+feature+'_files.txt')
  out_file_df.to_csv(out_file,sep='\t',header=None, index=False)
  print(out_file)
  return out_file_df.copy()

def plot_tss_enhancer_mr(wtype_str, tmp_gene, response, tss_mc, gbody_mc,tes_mc, gX,gY,G, w, window, sigma, MRs_txt,folder_out, is_plot2regions =True):
  '''
    tmp_gene: selected gene for ploting
    response: positive or negative respnose between gene and predicted regulatory regions śuch tss and enhancers
    tss_mc: tss methylation data
    gbody_mc: gene body methylatoin data or ennhancer methylaton data if only two regions are ploted
    tes_mc: tes methylation data
    gX , gY: gene upstream and downstram, or tss range, default=500
    G: gene length or enhancer range if only two regions are ploted , default= 2000 .
    window: window for combined regions, increase this parameter >500 if there is a gap between two combined regions, default=600
    sigma: data smoothing step, default=50
    MRs_txt: output file postprefix
    folder_out: output folder for plots
    is_plot2regions: whether to plot 2 or more regions
   
  '''
  feature=response
  check_folder(folder_out)
  folder_out=os.path.join(folder_out,tmp_gene+'_'+ feature)
  check_folder(folder_out)

  folder_out2=os.path.join(folder_out,'plotData')
  check_folder(folder_out2)
  #plot2regions=True
  plt_tgt_main(wtype_str, tss_mc, tes_mc, gbody_mc, gX, gY, G, w, window, sigma, MRs_txt, folder_out,is_plot2regions)


def run(args):
  #########################
  #1. gather dmr, deg information for predicted dmr vs target gene from dds package
  #jbw
  out_folder=args.out_data_folder
  check_folder(out_folder)

  #DEG file from bpb3 by rr-ratio filtered genes:  bpb3 filterDEG4bpb3 --help
  #exp_file='HAP1_P1_vs_HAP1_KO1_differentially_expressed_genes_min1.1Fd_min1RPKM_rratio_filtered.csv'
  #DMR file from dmr_analysis: dmr_analysis  dmr_combine_multChrs4rank 
  #dmr_file='24_chroms_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_2_isModTest_0__all_dmrRanking_top_0.92_minLogReg_proba_0.7_DMRs.bed.gz'
  #predicted tss, enhancer target genes from DDS: dds_analyssís dTarget_methy_vs_express
  #tss_file='dmr_in_deg_tss_10000sampling.csv'
  #enc_file='dmr_in_deg_enhancer_10000sampling.csv'
  
  exp_file=args.differential_gene_expression_file
  dmr_file=args.differential_methylation_file
  tss_file=args.tss_target_file
  enc_file=args.enhancer_target_file
  wtype_str=args.wildType_fileString  

  #1. find DEG and DMR information for predicted target genes
  #read data
  print('\nRead DEG, DMR inoformation for predicted target genes')
  is_dmr_compressed=args.dmr_file_not_compressed
  is_flip_sign_of_rr_ratio=args.not_flip_sign_of_rratio
  #print(is_dmr_compressed, is_flip_sign_of_rr_ratio)
  #jbw
  out_tss_file, out_enc_file, out_tss_df,out_enc_df=find_deg_dmr4target_genes(out_folder, tss_file, enc_file,exp_file,dmr_file, is_dmr_compressed , is_flip_sign_of_rr_ratio)

  #2. Extract data for selected genes for ploting their mean methylation levels in tss and enhancers.
  print('\nPrepare data for selected target for plotting their mean methylation levels in TSS/Enhancers')
  selected_genes=args.selected_genes4plot.split(',')

  #jbw
  print(out_folder,out_tss_file, out_enc_file)
  extract_data4selected_genes(out_folder, out_tss_file, out_enc_file,selected_genes)

  #3. genreate hmst-seq-analyzer format DMRs files for making plots in tss and enhancer
  print('\nGenerate hmst-seq-analyzer format of data for ploting .... ')
  tss_file3=os.path.join(out_folder, 'selected_genes_'+'tss' +'.tsv')
  enc_file3=os.path.join(out_folder, 'selected_genes_'+'enhancer' +'.tsv')
  #target gene with negative or positive correlation between expression and methylation levels
  isNegative=args.is_negative_target

  if isNegative==1:
    response='neg'
  elif isNegative==0:
    response='pos'
  else :
    response='all'
  print('Plot ' + response + ' response targets')

  mr_folder=args.mr_data_folder
  folder_name=args.prefix_foldername4chrom
  flank=args.flank_region

  gX=args.gene_upstream_X
  gY=args.gene_downstream_Y
  G=args.gene_or_enhancer_length
  w=args.step_size_in_region
  window=args.window_size_in_combined_regions
  sigma=args.data_smooth_sigma
  MRs_txt=args.MRs_txt_description
  MRs_txt= response + '_' + MRs_txt 
  plot_folder_out=os.path.join(out_folder, 'plots')
  data_folder_out=os.path.join(out_folder,'data')

  print('\nStart to plot figures ... ')
  for tmp_gene in selected_genes:
    #jbw
    generate_hmtseq_data(wtype_str,tmp_gene, flank, 'tss', tss_file3, isNegative, mr_folder, folder_name,data_folder_out )
    generate_hmtseq_data(wtype_str,tmp_gene, flank, 'enhancer', enc_file3, isNegative, mr_folder, folder_name,data_folder_out )

    print('\n Plot -> ',tmp_gene)
    #4. plot methylation level at TSS and Enhancer at selected gene 
    list_file_tss=os.path.join(data_folder_out, 'list_of_'+tmp_gene +'_'+response +'_tss_files.txt')
    if os.stat(list_file_tss).st_size==0:
      tss_mc=[]
    else:
      file_df=pd.read_csv(list_file_tss,header=None)
      tss_mc=file_df[0].to_list()

    list_file_enhancer=os.path.join(data_folder_out, 'list_of_'+tmp_gene +'_'+response+'_enhancer_files.txt')
    if os.stat(list_file_enhancer).st_size==0:
      gbody_mc=[]
    else:
      file_df=pd.read_csv(list_file_enhancer,header=None)
      gbody_mc=file_df[0].to_list()

    tes_mc=[]
    plot_tss_enhancer_mr(wtype_str, tmp_gene, response, tss_mc, gbody_mc,tes_mc, gX,gY,G, w, window, sigma, MRs_txt,plot_folder_out, is_plot2regions =True)
    print('\n')
  return

if __name__ == '__main__' :
  args=my_parser(argparse.ArgumentParser('python plot_tss_enhancer_mrs.py')).parse_args()
  run(args)

