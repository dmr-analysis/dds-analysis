#this script is used to check whether a 5dist gene and a mutation block is located in the same TAD or not
#if not then check whether the mutation block is located in the boundaries regions or not
import pandas as pd
import numpy as np
import os
import subprocess
#exec(open("check_block_dist_gene_in_TAD2.py").read())

def bedtools_intersect2files(file1,file2,min_overlap,out_file,isMultiple=False):
    command1='bedtools intersect -a ' + file1 + ' -b ' + file2 + ' -wa -wb -f ' + str(min_overlap) + '> '+ out_file

    #os may not wait for the process to finish
    #os.system(command1)

    if not isMultiple :
      #use process to wait for command to finish before run the next one
      cmd=subprocess.Popen(command1, shell=True)
      cmd.communicate()
    else:
      #use process to wait for multipe ones
      cmd=subprocess.Popen(command1,shell=True)

    return cmd

def read_file2df(in_file): 
       #read tmp_tad_out2gene, block, tmp_tod_out2gene,block
       if os.path.getsize(in_file) >0:
          tmp_out_tad_df2gene=pd.read_csv(in_file,sep='\t',header=None)
          tmp_out_tad_df2gene.columns=['chrs','start_pos','end_pos','type','length','bl_chrs','bl_start_pos','bl_end_pos','bl_type']
       else:
          tmp_out_tad_df2gene=pd.DataFrame()
       
       return tmp_out_tad_df2gene

def add_element2new_df_column(index,tmp_out_tad_df2gene,in_block_gene_df,new_column_label):
       if not tmp_out_tad_df2gene.empty:
          tmp_tad=tmp_out_tad_df2gene.loc[:,['chrs','start_pos','end_pos','type']].copy()
          record_tad=[]
          for index0, row0 in tmp_tad.iterrows():
             record_tad.append(':'.join(map(str,row0.to_list())))
          in_block_gene_df.loc[index,[new_column_label]] ='|'.join(record_tad)
       else:
          in_block_gene_df.loc[index,[new_column_label]]='na'
       return in_block_gene_df

def intersect2columns_from_df(index,in_block_gene_df,column1_label,column2_label):
      tmp_tad2gene=in_block_gene_df.loc[index,[column1_label]].copy()
      tmp_tad2block=in_block_gene_df.loc[index,[column2_label]].copy()

      tmp_tad2gene=tmp_tad2gene[column1_label].split('|')
      tmp_tad2block=tmp_tad2block[column2_label].split('|')
      if 'na' in tmp_tad2gene:
              tmp_tad2gene.remove('na')
      if 'na' in tmp_tad2block:
              tmp_tad2block.remove('na')
      intersect_gene2block=set(tmp_tad2block).intersection(set(tmp_tad2gene))
      return intersect_gene2block

def find_TADorBoundary_for_genes(in_block_gene_file,in_block_gene_df,in_gene_pos_df,in_tad_pos_df,in_boundary_df,min_overlap):
    for index,row in in_block_gene_df.iterrows():
  #    print(row)
       #check 5dist with TAD
       tmp_gene=row.gene_name
       tmp_gene_pos=in_gene_pos_df[in_gene_pos_df.name.str.contains('\|'+tmp_gene+':')].loc[:,['chrs','start_pos','end_pos']]
       tmp_gene_pos=tmp_gene_pos.to_numpy().tolist()[0]
       
       tmp_block=row.block_id.split('_')
       tmp_block_pos=list(map(str,tmp_block[-3:]))
      
       #TAD bed file 
       tmp_chr=tmp_gene_pos[0]
       tmp_tad_pos=in_tad_pos_df[in_tad_pos_df.chrs==tmp_chr]
       tmp_out_1='tmp_tad_pos.bed'
       tmp_tad_pos.to_csv(tmp_out_1,sep='\t',header=None,index=False)   

       #TAD boundary bed file
       tmp_boundary_pos=in_boundary_df[in_boundary_df.chrs==tmp_chr]
       tmp_out_3='tmp_boundary_pos.bed'
       tmp_boundary_pos.to_csv(tmp_out_3,sep='\t', header=None, index=False)

       #gene and block position bed file
       tmp_out_df_gene=pd.DataFrame(columns=['chrs','start_pos','end_pos','name'])
       tmp_out_df_block=pd.DataFrame(columns=['chrs','start_pos','end_pos','name'])
       tmp_out_df_gene.loc[0]=tmp_gene_pos+['gene']
       tmp_out_df_block.loc[0]=tmp_block_pos+['block']
       tmp_out_2gene='tmp_gene.bed'
       tmp_out_2block='tmp_block.bed'
       tmp_out_df_gene.to_csv(tmp_out_2gene,sep='\t',header=None,index=False)
       tmp_out_df_block.to_csv(tmp_out_2block,sep='\t',header=None,index=False)
 
       #TAD vs gene
       isMultiple=True
       processes=[]
       processes.append(bedtools_intersect2files(tmp_out_1,tmp_out_2gene,min_overlap,'tmp_tad_out2gene',isMultiple))

       #TAD vs block
       processes.append(bedtools_intersect2files(tmp_out_1,tmp_out_2block,min_overlap,'tmp_tad_out2block',isMultiple))

       #boundary vs gene
       processes.append(bedtools_intersect2files(tmp_out_3,tmp_out_2gene,min_overlap,'tmp_bod_out2gene',isMultiple))

       #boundary vs block
       processes.append(bedtools_intersect2files(tmp_out_3,tmp_out_2block,min_overlap,'tmp_bod_out2block',isMultiple))

       if isMultiple:
         while len(processes)>0:
           for i in range(len(processes)):
              if processes[i].poll() is not None:
                 del processes[i]
                 break
         for p in processes:
           p.wait()

       #read tmp_tad_out2gene, block, tmp_tod_out2gene,block
       tmp_out_tad_df2gene=read_file2df('tmp_tad_out2gene')
       tmp_out_tad_df2block=read_file2df('tmp_tad_out2block')
       tmp_out_bod_df2gene=read_file2df('tmp_bod_out2gene')
       tmp_out_bod_df2block=read_file2df('tmp_bod_out2block')

       #remove all tempary files
       os.system('rm -f tmp_gene.bed tmp_block.bed tmp_tad_pos.bed tmp_boundary_pos.bed tmp_tad_out2gene tmp_tad_out2block tmp_bod_out2gene tmp_bod_out2block')
 
#       input('Click') 
       #assign values
       add_element2new_df_column(index,tmp_out_tad_df2gene,in_block_gene_df,'TAD2gene')
       add_element2new_df_column(index,tmp_out_bod_df2gene,in_block_gene_df,'Boundary2gene')
       add_element2new_df_column(index,tmp_out_tad_df2block,in_block_gene_df,'TAD2block')
       add_element2new_df_column(index,tmp_out_bod_df2block,in_block_gene_df,'Boundary2block')

       #check gene and Block with the same TAD/Boundary or not.
       if (in_block_gene_df.TAD2gene[index]=='na') & (in_block_gene_df.TAD2block[index]=='na') \
               & (in_block_gene_df.Boundary2gene[index]=='na') & (in_block_gene_df.Boundary2block[index]=='na') :
               in_block_gene_df.loc[index,['isTAD']]=-1
       else:

          intersect_gene2block=intersect2columns_from_df(index,in_block_gene_df,'TAD2gene','TAD2block')
          intersect_bod_gene2block=intersect2columns_from_df(index,in_block_gene_df,'Boundary2gene','Boundary2block')
 
          if len(intersect_gene2block) >0:
            in_block_gene_df.loc[index,['isTAD']]=1
          elif len(intersect_bod_gene2block)>0:
            in_block_gene_df.loc[index,['isTAD']]=2
          else:
            in_block_gene_df.loc[index,['isTAD']]=0

    #export
    if '.csv' in in_block_gene_file:
       out_file=in_block_gene_file.replace('.csv','_commonTAD_Boundary.tsv')
    elif '.tsv' in in_block_gene_file:
       out_file=in_block_gene_file.replace('.tsv','_commonTAD_Boundary.tsv')
    print(out_file)
    in_block_gene_df.to_csv(out_file,sep='\t',index=False)
    return in_block_gene_df, out_file

if __name__=='__main__':
  in_gene_pos_file='in_tad/gene_Up1000_Down1000removedShort_adjustedChrs.bed'
  in_tad_pos_file='in_tad/Table4_TAD_annotations_sorted.bed'
  in_boundary_pos_file='in_tad/Table1_common_boundaries_merged_sorted.bed'
  #in_block_gene_file='out_gene_info/patient_13_gene_info.tsv'
  #in_block_gene_file='out_gene_info_p05/patient_1_gene_info.csv'

  in_gene_pos_df=pd.read_csv(in_gene_pos_file, sep='\t',header=None)
  in_gene_pos_df.columns=['chrs','start_pos','end_pos','name']

  in_tad_pos_df=pd.read_csv(in_tad_pos_file,sep='\t',header=None)
  in_tad_pos_df.columns=['chrs','start_pos','end_pos','chrom_type','tad_length']

  in_boundary_df=pd.read_csv(in_boundary_pos_file,sep='\t',header=None)
  in_boundary_df.columns=['chrs','start_pos','end_pos','name','length']
  
  for i in range(0,14):
     in_block_gene_file='out_gene_info_p05/patient_' + str(i) + '_gene_info.csv'
     print(in_block_gene_file)
     in_block_gene_df=pd.read_csv(in_block_gene_file,sep='\t')


     min_overlap=1e-9
     #isTAD=-1,0,1,2,not find, not same TAD, same TAD, same TAD boundary
     in_block_gene_df, out_file_name=find_TADorBoundary_for_genes(in_block_gene_file,in_block_gene_df,in_gene_pos_df,in_tad_pos_df,in_boundary_df,min_overlap)

     #export gene list within the same TAD or Boundary
     #tmp_gene=sorted(set(in_block_gene_df[in_block_gene_df.isTAD>0].gene_name.to_list()))
     #gene_list_df=pd.DataFrame(data=tmp_gene)

     tmp_gene_df=in_block_gene_df.loc[in_block_gene_df.isTAD>0,['gene_name','gene_type','block_id','isTAD']].copy()
     tmp_gene_df=tmp_gene_df.sort_values(by='gene_name').drop_duplicates(['gene_name','gene_type'])
     out_file2=in_block_gene_file.replace('.csv','_commonTAD_Boundary_list.tsv')
     tmp_gene_df.to_csv(out_file2,sep='\t',index=False)
     print(out_file2)

#sorted(set(in_block_gene_df[in_block_gene_df.isTAD==0].gene_name.to_list()))
