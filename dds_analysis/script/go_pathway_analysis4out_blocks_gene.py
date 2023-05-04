#this script is useed to count how many patients are affecting a GO category or signalling pathway
import pandas as pd
import numpy as np
import argparse
#in_gene_file='p_ls_05/gene_count_in_patients_geneType_blockID.tsv'
#in_GO_file='p_ls_05/out_DAVID/DAVID_GO_3555_genes_p_ls05.csv'
#exec(open("go_pathway_analysis4out_blocks_geneApril.py").read())
#this function is not tested after changing in may 7!

def my_parser(parser):
  required=parser.add_argument_group('Required')
  required.add_argument('-in_GEfile','--in_gene_file',required=True, help='unique gene file after combing DMR DEG and TAD and mutaiton block information ')
  required.add_argument('-in_GOfile','--in_GO_file', required=True, help='DAVID export file')
  optional = parser.add_argument_group('Optional, with default values')
  optional.add_argument('-cutoff','--cutoff_Pvalue4GO',default=0.05, metavar='', type=float, help='cutoff p-value for GO enrichment terms, default= 0.05')
  return parser

def run(args):
 typeStr='DMR_or_DEG'
 #in_gene_file='out_blocks_gene_single/blocks_summary_block_position_0flank_0.7Proba_66868blocks_13143blocks2mr_4604blocks2dmr_deg_info_filtered_'+ typeStr +'_uniqGene_commonTAD_Boundary_list2UqGene.tsv'
 #in_GO_file='out_blocks_DAVID/block_single_0flank_newDEG/tad_'+typeStr.lower() +'/DAVID_out_TAD_'+ typeStr+ '_mussd_single.csv'
 in_gene_file=args.in_gene_file
 in_GO_file=args.in_GO_file
 cut_off=args.cutoff_Pvalue4GO

 gene_list_df=pd.read_csv(in_gene_file,sep='\t')
 go_list_df=pd.read_csv(in_GO_file, sep='\t')
 gene_list_df.gene_name=gene_list_df.gene_name.str.upper()

 #for each row in GO list to count number of patieints are involved in it.
 loop=0
 record_all_genes=[]
 record_pathways=[]
 record_patients=[]
 #cut_off=0.05

 print("Enrichment in GOTERM_BP, KEGG, and BIOCARTA")
 for idx,row in go_list_df.iterrows():
  #if ('GOTERM_BP' in row.Category or 'KEGG_' in row.Category or 'BIOCARTA' in row.Category) and ( 'signaling' in str.lower(row.Term) or 'pathway' in str.lower(row.Term) ) and (row.PValue<cut_off):
  if ( 'GOTERM_BP' in row.Category or 'KEGG_' in row.Category or 'BIOCARTA' in row.Category) and (row.PValue<cut_off):
    tmp_gene=str.upper(row.Genes)
    tmp_df=pd.DataFrame(data=tmp_gene.split(','),columns=['gene_name'])
    tmp_df.gene_name=tmp_df.gene_name.str.replace(' ','')
 #may need change of counting patieints??
    tmp_df2=pd.merge(gene_list_df,tmp_df,on='gene_name',how='inner')
    s1=set(tmp_df.gene_name.to_list())
    s2=set(tmp_df2.gene_name.to_list())
 #    sum_of_gene_in_patient=np.array(tmp_df2.iloc[:,1:15].sum().to_list())
 #    total_patients=sum_of_gene_in_patient[sum_of_gene_in_patient>0].shape[0]
    #count number of patients are affected by these genes
    tmp_counts=tmp_df2.patient_id.apply(lambda x : set(x.split('~'))).to_list()
    y=set()
    [y.update(x) for x in tmp_counts]
    y2=set()
    [y2.update(set(x.split(','))) for x in y ]
    total_patients=len(y2)

    #if 'GO:0008630~' in row.Term:
    #    print(row.Term, total_patients) 
    #    input(' ')
 #may end change
    #print(s1-s2)
    #print(row.Term)
    #print(total_patients)
    #input('Click')
    loop +=1
    record_all_genes += list(s2)
    record_pathways.append(row.Term)
    record_patients.append(total_patients)

 record_all_genes=np.array(record_all_genes)
 record_pathways=np.array(record_pathways)
 record_patients=np.array(record_patients)

 tmp_gene_df=pd.DataFrame(data=list(set(record_all_genes)),columns=['gene_name'])
 tmp_GO_df=pd.DataFrame(data=list(record_pathways),columns=['Term'])
 #bug here after set and list the np.array again, elements unordered ?
 tmp_GO_df['total_patients']=record_patients

 out_gene_df=pd.merge(gene_list_df,tmp_gene_df,on='gene_name',how='inner')
 out_GO_df=pd.merge(go_list_df,tmp_GO_df,on='Term',how='inner')

 out_gene_file='out_gene_'+ typeStr + '_in_GO_KEGG_BIOC_pathways_p_ls'+str(cut_off)+'.csv'
 out_GO_file='out_GO_'+ typeStr +'_in_GO_KEGG_BIOC_pathways_p_ls'+ str(cut_off)+'.csv'

 out_gene_df.to_csv(out_gene_file,sep='\t',index=False)
 out_GO_df=out_GO_df.sort_values(by='total_patients', ascending=False)
 out_GO_df.to_csv(out_GO_file,sep='\t',index=False)

 out_GO_df[(out_GO_df.total_patients>=10) & (out_GO_df.PValue<0.05)]


if __name__=='__main__':
 args=my_parser(argparse.ArgumentParser('python go_pathway_analysis4out_blocks_gene.py ')).parse_args()
 run(args)


