import argparse
import sys
import os

class Main(object):
  def __init__(self):
     parser= argparse.ArgumentParser(
             description ='DDS-Analysis: An integrated data analysis pipeline by considering both Differential '
                          'Methylation Region (DMR) and Differential Expression Gene (DEG) in SNP analysis. '
                          'This package dependents on some of functinos from bpb3 package! ',
             usage= ''' dds_analysis <task> [<args>]

     Tasks available for using:
        bpb3summary2bed_format	Convert bpb3 block summary file to a bed format file
        map_block2genome 	Map mutation block to genomic regions
        map_block2chromSegment 	Map mutation block to chromatin Segmenent regions
        map_block2dmr 	Map mutation block to differential methylated regions
        find_geneExp4block 	Find differential expressed genes for mutation blocks
        find_block_patieintID 	Find patient ID for mutatin blocks
        combine_dmr_deg2block 	Combine DMR, DEG, and mutation block information together
        filter_blocks 	Filter mutation blocks by using DMR or/and DEG condition 
        collect_gene_names4blocks 	Collect unique gene names for mutation blocks with DMR and/or DEG
        check_block_gene_inTAD 	Check whether block and gene are in the same TAD or TAD boundary
        dds_geneRanking 	Select top ranked genes from final prediction
        go_pathway_analysis4out_blocks_gene 	GO pathway analysis of genes
        find_enhancer_target_genes 	find enhancer and its target genes overlapping with mutation bocks that associated with selected gene
        chromSegment_test4blocks	Enrichment test of mutation blocks or methylation regions that associated with genes in 7 chromatin segmentations of human genome  
        dTarget_methy_vs_express	predict long distance target gene for a specific region (e.g., mutation block) based on coupling of methylation and gene expression across samples 
        plot_mr_vs_exp			Plot DMR/MR methylation level and Gene expression for a pair of DMR and its target gene
        plot_tss_enhancer_mrs		Plot the average methylation level of predicted DMRs at TSS and enhancer regions by the target genes predicted from dTarget_methy_vs_express 
	filterDEG4bpb3		Filter Differential Expressed Genes (DEG) by rratio based on exported file from bpb3 differential_expression then export it with group mean and rratio
	preprocess          This module first find DEG in TSS, 5dist regions then preprocess data for dds_analysis.	
  ''')
     parser.add_argument('task', help= 'Pipeline task to run ')
     args= parser.parse_args(sys.argv[1:2])
     if not hasattr(self, args.task):
        print('****Error: Unrecognized task ******')
        parser.print_help()
        exit(1)
     getattr(self, args.task)()
  
  def bpb3summary2bed_format(self):
      from .script.bpb3summary2bed_format import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='dds_analysis bpb3summary2bed_format', 
                description='Convert bpb3 block summary file to a bed format file, where the forth column is block ID, only works for human genome.'))
      run(parser.parse_args(sys.argv[2:]))
    
  def map_block2genome(self):
      from .script.map_block2genome import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='dds_analysis map_block2genome',
            description='Map mutatin blocks to genomic regions '))
      run(parser.parse_args(sys.argv[2:]))

  def map_block2chromSegment(self):
      from .script.map_block2chromSegment import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='dds_analysis map_block2chromSegment',
            description='Map mutatin blocks to predicted chromation segment regions '))
      run(parser.parse_args(sys.argv[2:]))

  def map_block2dmr(self):
      from .script.map_block2dmr import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='dds_analysis map_block2dmr',
            description='Map mutatin blocks to differential methylated regions '))
      run(parser.parse_args(sys.argv[2:]))

  def find_geneExp4block(self):
      from .script.find_geneExp4block import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='dds_analysis find_geneExp4block',
            description='Find differential gene expression for mutation blocks '))
      run(parser.parse_args(sys.argv[2:]))

  def find_block_patientID(self):
      from .script.find_block_patientID import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='dds_analysis find_block_patientID',
            description='Find mutation blocks patieint IDs '))
      run(parser.parse_args(sys.argv[2:]))

  def combine_dmr_deg2block(self):
      from .script.combine_dmr_deg2block import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='dds_analysis combine_dmr_deg2block',
            description='Combine information from DMR, DEG, mutation blocks '))
      run(parser.parse_args(sys.argv[2:]))

  def filter_blocks(self):
      from .script.filter_blocks import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='dds_analysis filter_blocks',
            description='Filter mutation blocks by using DMR and/or DEG information '))
      run(parser.parse_args(sys.argv[2:]))

  def collect_gene_names4blocks(self):
      from .script.collect_gene_names4blocks import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='dds_analysis collect_gene_names4blocks',
            description='Collect unique gene names for mutatino blocks '))
      run(parser.parse_args(sys.argv[2:]))

  def check_block_gene_inTAD(self):
      from .script.check_block_gene_inTAD import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='dds_analysis check_block_gene_inTAD',
            description='Check whether mutatino blocks and genes are in the same TAD or TAD boundary '))
      run(parser.parse_args(sys.argv[2:]))
  
  def dds_geneRanking(self):
      from .script.dds_geneRanking import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='dds_analysis dds_geneRanking',
               description='Select top ranked genes and blocks from final prediction ')) 
      run(parser.parse_args(sys.argv[2:]))

  def go_pathway_analysis4out_blocks_gene(self):
      from .script.go_pathway_analysis4out_blocks_gene import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='dds_analysis go_pathway_analysis4out_blocks_gene',
            description='GO pathway analysis for genes associated to mutaiton blocks '))
      run(parser.parse_args(sys.argv[2:]))

  def find_enhancer_target_genes(self):
      from .script.find_enhancer_target_genes import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='dds_analysis find_enhancer_target_genes',
                description='Find enhancer and its target genes are overlapping with mutation blocks that assocaited with a selected gene'))
      run(parser.parse_args(sys.argv[2:]))
   
  def chromSegment_test4blocks(self):
      from .script.chromSegment_test4blocks import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='dds_analysis chromSegment_test4blocks',
                         description='Enrichment test of mutation blocks or methylation regions assocated with genes in 7 chromatin segmentations of human genome'))
      run(parser.parse_args(sys.argv[2:]))

  def dTarget_methy_vs_express(self):
      from .script.dTarget_methy_vs_express import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='dds_analysis dTarget_methy_vs_express',
                description='Predict long distance target gene for a sepcific region (e.g., mutation block) based on coupling of methylation and gene expression across samples'))
      run(parser.parse_args(sys.argv[2:]))
  
  def plot_mr_vs_exp(self):
      from .script.plot_mr_vs_exp import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='dds_analysis plot_mr_vs_exp',
               description='Plot DMR/MR methylation levels and gene expression for a pair of DMR and its target gene'))
      run(parser.parse_args(sys.argv[2:]))
  
  def plot_tss_enhancer_mrs(self):
     from .script.plot_tss_enhancer_mrs import my_parser, run
     parser = my_parser(argparse.ArgumentParser(prog='dds_analysis plot_tss_enhancer_mrs',
              description='Plot the average DMR methylation levels at TSS and enhancer regions based on their predicted target genes from dTarget_methy_vs_express, '\
                          'and export differential gene and differential methylation information for each pair of putative targets and TSS or enhancers.'))
     run(parser.parse_args(sys.argv[2:]))

  def filterDEG4bpb3(self):
     from .script.filterDEG4bpb3 import my_parser, run
     parser = my_parser(argparse.ArgumentParser(prog='dds_analysis filterDEG4bpb3',
		description= 'Filter Differential Expressed Genes (DEG) by rratio based on exported file from bpb3 differential_expression, then export it with group mean and rratio.' ))
     run(parser.parse_args(sys.argv[2:]))

  def preprocess(self):
     from .script.preprocess import my_parser, run
     parser = my_parser(argparse.ArgumentParser(prog='dds_analysis preprocess',
		description= 'This module first find DEG in TSS, 5dist regions then preprocess data for dds_analysis.' ))
     run(parser.parse_args(sys.argv[2:]))

 
def main():
   Main()
   #jbw 2025
   os._exit(os.EX_OK)

if __name__ == '__main__':
   Main()






