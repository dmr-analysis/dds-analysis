demoHAP1
========

Step 1: dds_analysis preprocess
_______________________________
The `dds_analysis preprocess` step prepares the input files for further analysis using the `dds_analysis` module. It requires specifying various parameters and input file paths. Here is the code:

.. code-block:: bash

   # Before running following steps, it assumes that DMRs are already predicted by dmr_analysis

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

Set paths:
----------

 .. code-block:: bash
    #main path of input data
    IN_DATA_PATH='../../data/hap1_cell/in_data/final_demo_data/hap1_cell/'

    #path of DMR results from dmr_analysis
    IN_MR_PATH=${IN_DATA_PATH}'/out_data/DMR_CpG_context/'

    #path of DEG results from bpb3
    IN_DEG_PATH=${IN_DATA_PATH}'/in_data/DEG/'

    #DEG file name from bpb3 differential_analysis, the original DEF file from bpb3 that was used to convert Zscores in dds_analysis preprocess
    IN_DEG_FILE='HAP1_P1_vs_HAP1_KO1_differentially_expressed_genes_min0Fd_min0RPKM.txt '
    in_data_str='_hap1'

    #path to output data
    OUT_PATH='../../data/hap1_cell/out_data/'

    #path to exported MRs that are not located in TSS or enhancer regions
    FILE_FOLD=${OUT_PATH}/out4mr_not_in_tss_enhancer
    #file name for background sample list that contain all MRs not located in TSS or enhancers
    BACK_FILE=${OUT_PATH}/background_samples_list.tsv

    #whether to skip below two steps in the pipeline
    is_run_dmr_export=1 # 1 for exporting, other values for skipping this step
    is_run_dtarget=1    # 1 for run dTarget prediction , other values for skipping this step



Step 2: DMR data export
_______________________
The DMR data export step involves exporting DMR data that is either located in TSS or 5'distance regions. It uses the `dmr_exportData` command from the `dmr_analysis` module. Here is the code:

.. code-block:: bash

   # Before running this step, ensure DMRs are already predicted by dmr_analysis

   dmr_analysis dmr_exportData \
         --input_mr_data_folder ${IN_MR_PATH} \
         --output_file_folder ${OUT_PATH}/out4dmr_in_deg_tss_5dist \
         --input_file_format 0 \
         --number_of_processes 10 --input_file ${OUT_PATH}/uqdmr_regions_in_deg_tss_5dist${in_data_str}.bed -wtStr 'HAP1_P_'

   echo "Export data of DMRs overlapping to TSS or 5distance - Done"

   # Export data of MRs not in TSS or enhancers
   dmr_analysis dmr_exportData  \
         --input_mr_data_folder ${IN_MR_PATH} \
         --output_file_folder ${OUT_PATH}/out4background \
         --input_file_format 0 \
         --number_of_processes 10 --input_file ${OUT_PATH}/uqdmr_regions_not_in_deg_tss_5dist${in_data_str}.bed -wtStr 'HAP1_P_'

   echo "Export data of MRs not in TSS or enhancers - Done"

Step 3: Background file list creation
_____________________________________
The background file list creation step involves creating a background file list if it doesn't already exist. The background file list contains all MRs that are not located in TSS or enhancer regions. Here is the code:

.. code-block:: bash

   background_filelist="${OUT_PATH}/backgroundFileList.txt"

   if [[ ! -f "$background_filelist" ]]; then
       cd ${OUT_PATH}/out4background
       ls | grep ".bed" > $background_filelist
       cd -
       echo "Creating a list of background files - Done"
   else
       echo "Background file list already exists."
   fi

   echo "Background file list creation - Done"

Step 4: dTarget prediction
__________________________

The dTarget prediction step predicts putative target genes for DMRs based on gene expression profiles. It performs the prediction separately for DMRs associated with TSS regions and 5'distance regions. Here is the code:

.. code-block:: bash

   # Predict target genes for DMRs overlapping with TSS regions
   dds_analysis dTarget_methy_vs_express \
         -in_folder ${OUT_PATH}/out4dmr_in_deg_tss_5dist \
         -out_folder ${OUT_PATH}/out4dmr_in_deg_tss_5dist_out_targetGenes \
         -out_file_targetgenes ${OUT_PATH}/out4dmr_in_deg_tss_5dist_dTargetGenes.txt \
         -background_filelist ${OUT_PATH}/backgroundFileList.txt \
         -cpg_col_name '3' -fold_diff_threshold 1 -wt_str HAP1_P_

   echo "Target gene prediction for DMRs overlapping with TSS regions - Done"

   # Predict target genes for DMRs associated with 5'distance regions
   dds_analysis dTarget_methy_vs_express \
         -in_folder ${OUT_PATH}/out4dmr_in_deg_tss_5dist \
         -out_folder ${OUT_PATH}/out4dmr_in_deg_tss_5dist_out_targetGenes \
         -out_file_targetgenes ${OUT_PATH}/out4dmr_in_deg_tss_5dist_dTargetGenes.txt \
         -background_filelist ${OUT_PATH}/backgroundFileList.txt \
         -cpg_col_name '3' -fold_diff_threshold 1 -wt_str HAP1_P_

   echo "Target gene prediction for DMRs associated with 5'distance regions - Done"

Step 5: Plotting selected target gene and DMR associations
__________________________________________________________
This step involves plotting the associations between selected target genes and DMRs based on gene expression profiles. Here is the code:

.. code-block:: bash

   dds_analysis plot_mr_vs_exp \
         -in_file ${OUT_PATH}/out4dmr_in_deg_tss_5dist_dTargetGenes.txt \
         -cutoff_pval 1 -cutoff_FDR 1 -cutoff_abs_diff_methy 0 -output_file ${OUT_PATH}/plot_MR_vs_exp_dTargetGenes.pdf

   echo "Plotting selected target gene and DMR associations - Done"

Step 6: Plotting average methylation pattern
____________________________________________

The final step involves plotting the average methylation pattern for the selected target genes and DMRs. Here is the code:

.. code-block:: bash

   dds_analysis plot_avg_methylation_pattern \
         -in_file ${OUT_PATH}/out4dmr_in_deg_tss_5dist_dTargetGenes.txt \
         -output_file ${OUT_PATH}/plot_avg_methylation_pattern_dTargetGenes.pdf

   echo "Plotting average methylation pattern - Done"

.. image:: 5mC_Enhancer_X2000_Y1000_G2000_binsize100_sigma50_neg_DMRs_2023-05-24.jpg
   :alt: Enhancer vs Methylation level
