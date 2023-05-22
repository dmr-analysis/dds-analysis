#!/bin/bash
#Before running following steps, it assumes that DMRs are already predicted by dmr_analysis
#
#the script is used to prepare files for dds_analysis

dds_analysis preprocess \
      -in_folder ../../data/hap1_cell/in_data/final_demo_data/hap1_cell/out_data/DMR_CpG_context/out_map2genome/ \
      -in_string '_hap1' \
      -in_tss_file_mr ../../data/hap1_cell/in_data/final_demo_data/hap1_cell/out_data/DMR_CpG_context/out_map2genome/3_chroms_all_mr_data_range_dmrRanking_TSS_Up5000_Down1000_removedShort_overlap1e-09.bed \
      -in_dist_file ../../data/hap1_cell/in_data/final_demo_data/hap1_cell/out_data/DMR_CpG_context/out_map2genome/3_chroms_all_mr_data_range_dmrRanking_noGenes_5dist_Up1000000_Up5000removedShort_overlap1e-09.bed \
      -in_deg_file ../../data/hap1_cell/in_data/final_demo_data/hap1_cell/in_data/DEG/HAP1_P1_vs_HAP1_KO1_differentially_expressed_genes_min1.1Fd_min1RPKM.txt \
      -out_folder ../../data/hap1_cell/out_data/ \
      -tss_file  ../../data/hap1_cell/in_data/final_demo_data/hap1_cell/out_data/DMR_CpG_context/data/TSS_Up5000_Down1000_removedShort.bed \
      -full_mr_file ../../data/hap1_cell/in_data/final_demo_data/hap1_cell/out_data/DMR_CpG_context/3_chroms_all_mr_data_range_dmrRanking.bed \
      -in_genome_file ../../data/hap1_cell/in_data//final_demo_data/genome/hg38/hg38_all_enhancers_merged_hglft_genome_327b3_4dmr.bed \
      -gene_col_name '#gene'


echo "To find DMR regions that are overlapping with TSS or 5distance regions of DEG - and preprocess Done"
 
