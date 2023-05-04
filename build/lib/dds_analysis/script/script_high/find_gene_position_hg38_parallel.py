#this script is used to find gene postion based on gtf file
import pandas as pd
import csv
import numpy as np
import multiprocessing as mp
import math
#exec(open("find_gene_position_hg38_parallel.py").read())

def find_feature_column(x,feature):
    tmp=x.split(';')
    out_f=[i.replace('"','') for i in tmp if feature in i]
    if len(out_f)>0:
       out_f=out_f[0].strip().split(" ")[1]
    else:
       #print(out_f, tmp)
       out_f='na'
      
    return out_f

def make_bins(row):
  #may be need to add 1 at the end of position , the old version did not??
  bin_pos=[i for i in range(row['start'] , row['end'].astype(int)+1)]
  #bin_len=len(bin_pos)
  return bin_pos

def find_gene_bins_position(gene_df,gene_id,bin_size):
  selected_df=gene_df[gene_df.gene_id==gene_id].copy()
  selected_df1=selected_df.drop_duplicates(['seqname','start','end']).copy()
  #uq_transcript=selected_df1.transcript_id.unique()
  #record_uq=[]
  #for i in uq_transcript:
  #    record_uq.append(selected_df1[selected_df1.transcript_id==i].shape[0])
  #record_uq=np.array(record_uq)
  #max_transcript=uq_transcript[record_uq==record_uq.max()][0]
  #selected_df2=selected_df1[selected_df1.transcript_id==max_transcript].copy()
  selected_df2=selected_df1.copy()
  selected_df2['bin_pos'] = selected_df2.apply(make_bins, axis=1)
  selected_df2['bin_len'] = selected_df2.bin_pos.apply(lambda x: len(x))
  record_bin=[]
  for idx, row in selected_df2.iterrows():
    record_bin=record_bin+ row.bin_pos
  uq_record_bin=list(set(record_bin))
  uq_record_bin.sort()
  record_bins=np.array(uq_record_bin)
  bins=[i for i in range(1,len(record_bins),bin_size)]
  record_bins[bins[0:]]
  gene_bin_pos=record_bins[bins[0:]]
  gene_bins=np.array(bins[0:])
  return gene_bins, gene_bin_pos,selected_df2


def divide_data2chunk(num_of_processes, list_of_ids):
 #divide all id to blocks based on num of processes
 #num_of_processes=15
 #all_blocks= gene_df.gene_id.unique()
 all_blocks=list_of_ids
 num_in_chunks=int(math.ceil(len(all_blocks)/num_of_processes))
 block_chunks= [all_blocks[x:x+num_in_chunks] for x in range(0, len(all_blocks), num_in_chunks)]

 #do find_gene_bins_position in each chunk of blocks
 if len(block_chunks) > num_of_processes:
      num_of_processes=len(block_chunks)
 elif len(block_chunks) < num_of_processes:
      print('Number of blocks smaller than the number of available processes ')
      num_of_processes=len(block_chunks)
 elif len(block_chunks)==1:
      num_of_processes= 1
 return block_chunks, num_of_processes

def parallel_find_gene_bins_position(args):
  gene_df,gene_ids,bin_size=args
  record_bins=[]
  record_bins_pos=[]
  record_bins_gene=[]
  record_bins_chrom=[]
  for gene_id in gene_ids:
    tmp_bins, tmp_bins_pos,tmp_df2=find_gene_bins_position(gene_df,gene_id,bin_size)
    record_bins.append(','.join(tmp_bins.astype(str)))
    record_bins_pos.append(','.join(tmp_bins_pos.astype(str)))
    record_bins_gene.append(gene_id)
    record_bins_chrom.append(tmp_df2.seqname.iloc[0])
 
  #make a dataframe for gene bins
  gene_id2bins_position_df=pd.DataFrame()
  gene_id2bins_position_df['gene_id']=record_bins_gene
  gene_id2bins_position_df['bins']=record_bins
  gene_id2bins_position_df['bins_position']=record_bins_pos
  gene_id2bins_position_df['seqname']=record_bins_chrom
  return gene_id2bins_position_df

def find_gene_bins_chrom_position(x,gene_id2bins_position_df):
  #x=new_ip_count_df.ids.iloc[0].strip()
  x_id, x_pos=x.split(',')
  tmp_bin_pos= gene_id2bins_position_df[gene_id2bins_position_df.gene_id==x_id]
  tmp_bins_chrom=tmp_bin_pos.seqname.to_list()[0]
  tmp_bins=np.array( tmp_bin_pos.bins.str.split(',').to_list()[0],dtype=int)
  tmp_bins_pos=np.array(tmp_bin_pos.bins_position.str.split(',').to_list()[0],dtype=int)
  delt_bins= np.abs(tmp_bins-int(x_pos))
  delt_bins_pos= tmp_bins_pos[np.where(delt_bins==delt_bins.min())]
  return tmp_bins_chrom, delt_bins_pos[0],x

if __name__=='__main__':
 #read gtf file
 #f1='hg38.ensGene.gtf.gz'
 f1='hg38_v25_ek12.gtf.gz'
 #f1='test_hg38.gtf.gz'
 in_df=pd.read_csv(f1,sep='\t',skiprows=5,compression='gzip',  header=None,quoting=csv.QUOTE_NONE)
 in_df.columns=['seqname','source','feature','start','end','score','strand','frame','attribute']

 #extract feature from gtf file
 #type_of_features=in_df.feature.unique()
 gene_df=in_df[in_df.feature=='exon'].copy()

 bin_size=100

 #find gene id for transcript_df and gene_df
 feature_name='gene_id'
 feature_id2=gene_df.attribute.apply(lambda x: find_feature_column(x,feature_name)).to_list()
 gene_df.insert(0,feature_name, feature_id2)

 feature_name='transcript_id'
 feature_id3=gene_df.attribute.apply(lambda x: find_feature_column(x,feature_name)).to_list()
 gene_df.insert(0,feature_name, feature_id3)

 #for each exon divide to bin
 #for each gene-id or transcript id find its bins for all concated exons 
 #divide all id to blocks based on num of processes
 num_of_processes=15
 all_blocks= gene_df.gene_id.unique()
 block_chunks,num_of_processes= divide_data2chunk(num_of_processes, all_blocks)

 pool = mp.Pool(processes=num_of_processes)
 files= pool.map(parallel_find_gene_bins_position, [ (gene_df,block_chunks[loop], bin_size) for loop in range(0, num_of_processes) ],1)
 #start to merge results
 gene_id2bins_position_df=pd.concat(files)
 pool.close()

 #export gene bin position file
 out_file=f1.replace('.gtf','_bins_positions.tsv')
 gene_id2bins_position_df.to_csv(out_file,sep='\t',index=False)
 print('Export at: ', out_file)

 #out_file='hg38_gene_id_bins_positions_old.tsv'
 #gene_id2bins_position_df=pd.read_csv(out_file,sep='\t')

 #read count table file
 #find gene id for IP count_df
 f2='IP-counts_normalized_adjusted_filtered.tab'
 ip_count_df=pd.read_csv(f2,sep='\t')
 ip_count_df.insert(0,'ids',ip_count_df.index.to_list())
 new_info_df=ip_count_df.ids.str.split(',',expand=True).rename(columns={0:'gene_id', 1:'bin_start_pos'})
 new_ip_count_df=pd.merge(ip_count_df,new_info_df,left_index=True, right_index=True).copy()

 pos2ip_count_df=pd.DataFrame(new_ip_count_df.ids.apply(lambda x: find_gene_bins_chrom_position(x,gene_id2bins_position_df)).to_list(), columns=['chrom','start_pos','ids'])
 pos2ip_count_df.insert(2,'end_pos',list(pos2ip_count_df.start_pos+bin_size),True)

 merged_ip_count_df=pd.merge(pos2ip_count_df,new_ip_count_df,on='ids')

 out_df=merged_ip_count_df.copy()

 #export dataframe
 out_file=f2.replace('.tab','_with_position.bed')

 #sort dataframe
 #out_df['num4chr']= out_df.seqname.apply(lambda x: x.replace('chr','')).to_list()
 out_df.insert(0,"num4chr",  out_df.chrom.apply(lambda x: x.replace('chr','')).to_list(), True)
 sorted_out_df=out_df.sort_values(by=['num4chr','start_pos']).copy()
 sorted_out_df.pop('num4chr')
 sorted_out_df2= sorted_out_df.drop_duplicates(['chrom','start_pos','end_pos','ids'])
 sorted_out_df2.to_csv(out_file,sep='\t',index=False)





