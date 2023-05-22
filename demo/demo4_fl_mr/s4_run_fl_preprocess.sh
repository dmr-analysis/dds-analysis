#!/bin/bash
#Before running following steps, it assumes that DMRs are already predicted by dmr_analysis
#
#the script is used to prepare files for dds_analysis dTarget_methy_vs_express

#A descripton of parameters
#in_folder:	is path of DMRs mapped to genome 
#in_string:	string in exported file
#in_tss_file_mr:	DMRs mapped to TSS regions
#in_dist_file:		DMRs mapped to 5distance regions
#in_deg_file:		differential expression file (bpb3 output format)
#out_folder:		output file path
#tss_file:		a predfined TSS regions for further analysis
#full_mr_file:		a list of all ranked dmr from dmr_analysis
#in_genom_file:		a bed formated enhancer postion file
#gene_col_name:		column name of gene ID/name


#Run start
dds_analysis preprocess \
      -in_folder  '../../data/fl_mr/out_data/out_map2genome/' \
      -in_string '_fl' \
      -in_tss_file_mr '../../data/fl_mr/out_data/out_map2genome/2_chroms_all_mr_data_range_dmrRanking_TSS_Up10000_Down10000_overlap1e-09.bed' \
      -in_dist_file '../../data/fl_mr/out_data/out_map2genome/2_chroms_all_mr_data_range_dmrRanking_noGenes_5dist_Up1000000_Up10000_overlap1e-09.bed' \
      -in_deg_file '../../data/fl_mr/in_data/in_DEG/differentially_expressed_genes_nonfolder.txt' \
      -out_folder '../../data/fl_mr/out_data/' \
      -tss_file '../../data/fl_mr/in_data/in_genome_regions/data/TSS_Up10000_Down10000.bed' \
      -full_mr_file '../../data/fl_mr/in_data/DMR_CpG_context/2_chroms_all_mr_data_range_dmrRanking.bed' \
      -in_genome_file '../../data/fl_mr/in_data/in_genome_regions/data/hg19_all_enhancers_merged_4dmr.bed' \
      -gene_col_name '#gene'


echo "To find DMR regions that are overlapping with TSS or 5distance regions of DEG - and preprocess Done"
 
