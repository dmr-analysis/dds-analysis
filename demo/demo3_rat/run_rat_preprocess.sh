#!/bin/bash
#Before running following steps, it assumes that DMRs are already predicted by dmr_analysis
#
#the script is used to prepare files for dds_analysis

dds_analysis preprocess \
      -in_folder ../../data/rat_data/in_data/final_demo_data/rat_data/out_data/DMR_CpG_context/out_map2genome/ \
      -in_string '_rat' \
      -in_tss_file_mr ../../data/rat_data/in_data/final_demo_data/rat_data/out_data/DMR_CpG_context/out_map2genome/5_chroms_all_mr_data_range_dmrRanking_TSS_Up5000_Down1000_removedShort_overlap1e-09.bed \
      -in_dist_file ../../data/rat_data/in_data/final_demo_data/rat_data/out_data/DMR_CpG_context/out_map2genome/5_chroms_all_mr_data_range_dmrRanking_noGenes_5dist_Up1000000_Up5000removedShort_overlap1e-09.bed\
      -in_deg_file ../../data/rat_data/in_data/final_demo_data/rat_data/in_data/DEG/Adrenal1vsAdrenal2_DEG_genes_zscores.tsv\
      -out_folder ../../data/rat_data/out_data/ \
      -tss_file ../../data/rat_data/in_data/final_demo_data/rat_data/out_data/DMR_CpG_context/data/TSS_Up5000_Down1000_removedShort.bed \
      -full_mr_file ../../data/rat_data/in_data/final_demo_data/rat_data/out_data/DMR_CpG_context/5_chroms_all_mr_data_range_dmrRanking.bed \
      -in_genome_file ../../data/rat_data/in_data//final_demo_data/genome/rn6/rn6.enhancers_all_rn5_merged_rn6liftOvered_4dmr.bed \
      -gene_col_name 'gene_name'


echo "To find DMR regions that are overlapping with TSS or 5distance regions of DEG - and preprocess Done"
 
