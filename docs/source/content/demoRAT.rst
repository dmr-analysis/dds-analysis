demoRAT
=======

Step 1: Preparing files for dds_analysis
________________________________________

The first step involves preparing files for `dds_analysis` by running the `preprocess` command. This command integrates DMR and DEG data and prepares the necessary input files for further analysis. Here is the code:

.. code-block:: bash

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

Step 2: Export data:
____________________

The second step involves running the `dmr_analysis dmr_exportData` command to export relevant methylation region data.

.. code-block:: bash

   # Define input paths
   IN_DATA_PATH='../../data/rat_data/in_data/final_demo_data/rat_data/'
   IN_MR_PATH=${IN_DATA_PATH}'/out_data/DMR_CpG_context/'
   IN_DEG_PATH=${IN_DATA_PATH}'/in_data/DEG/'

   # Define output path
   OUT_PATH='../../data/rat_data/out_data/'

   # Define file paths
   FILE_FOLD=${OUT_PATH}/out4mr_not_in_tss_enhancer
   BACK_FILE=${OUT_PATH}/background_samples_list.tsv

   # Set variables
   in_data_str='_rat'
   is_run_dmr_export=1
   is_run_dtarget=1

   # Export data for DMRs overlapping with TSS or 5'distance regions
   if [ $is_run_dmr_export == 1 ]; then
      dmr_analysis dmr_exportData \
            --input_mr_data_folder ${IN_MR_PATH} \
            --output_file_folder ${OUT_PATH}/out4dmr_in_deg_tss_5dist \
            --input_file_format 0 \
            --number_of_processes 10 --input_file ${OUT_PATH}'/uqdmr_regions_in_deg_tss_5dist'${in_data_str}'.bed' -wtStr '_Ctrl'
      echo "Export data of DMRs overlapping to TSS or 5distance - Done "
      echo ""

      # Export data for MRs that are not in TSS or enhancer regions
      dmr_analysis dmr_exportData  \
            --input_mr_data_folder ${IN_MR_PATH} \
            --output_file_folder ${OUT_PATH}/out4mr_not_in_tss_enhancer \
            --input_file_format 0 \
            --number_of_processes 10 --input_file ${OUT_PATH}'/mr_regions_not_in_enhancers'${in_data_str}'_tss.bed' -wtStr '_Ctrl'
      echo "Export data of MRs not in TSS or enhancers - Done "
      echo ""
   fi

   # Create background file list if it does not exist
   if ! [ -f $BACK_FILE ]; then
      echo $BACK_FILE " not exists and create one ! "
      if [ -e $FILE_FOLD ]; then
         ls  ./${FILE_FOLD}/chr*/data/*raw*.* > $BACK_FILE
         echo "Create " $BACK_FILE
      else
         echo "Cannot create background file because no data folder find! " $FILE_FOLD
      fi
   fi

Step 3: Running dds_analysis dTarget_methy_vs_express
_____________________________________________________
The third step involves running the `dds_analysis dTarget_methy_vs_express` command to predict putative target genes for DMRs based on their associations from either TSS or 5'distance regions. Here is the code:

.. code-block:: bash

   # Run dTarget_methy_vs_express for predicting target genes
   if [ $is_run_dtarget == 1 ]; then
      gene_mr_file=${OUT_PATH}'/uqGeneDmr_regions_in_deg_tss'${in_data_str}'.bed'
      gene_exp_file=${IN_DEG_PATH}'/Adrenal1vsAdrenal2_DEG_genes_zscores.tsv'
      in_mr_data_folder=${OUT_PATH}/out4dmr_in_deg_tss_5dist
      in_background_mr_file=$BACK_FILE
      number_of_samples=10

      # Test target gene and DMR associations from TSS regions
      dds_analysis dTarget_methy_vs_express -inGeneMRfile $gene_mr_file  -mrTAB \
            -inGeneEXPfile $gene_exp_file -expTAB \
            -inMRfolder $in_mr_data_folder -outName 'tss_region_' \
            -output_path $OUT_PATH -sampleName 'sample_name4replace.tsv' \
            -pathDepth 1 -inBackgroundList $in_background_mr_file -cutoff 0.05 -totalSamples $number_of_samples -numOfprocesses 10

      echo "Done with TSS target gene prediction"

      # Test target gene and DMR associations from 5'distance regions
      gene_mr_file=${OUT_PATH}'/uqGeneDmr_regions_in_deg_5dist'${in_data_str}'_overlap_enhancer.bed'
      dds_analysis dTarget_methy_vs_express -inGeneMRfile $gene_mr_file -mrTAB \
            -inGeneEXPfile $gene_exp_file -expTAB \
            -inMRfolder $in_mr_data_folder -outName 'distance_region_'  \
            -output_path $OUT_PATH -sampleName 'sample_name4replace.tsv' \
            -pathDepth 1 -inBackgroundList $in_background_mr_file -cutoff 0.01 -totalSamples $number_of_samples -numOfprocesses 10

      echo "Done with 5'distance target gene prediction"
   fi

Step 3: Plotting selected target gene and DMR associations
__________________________________________________________

.. code-block:: bash

    gene_exp_file=${IN_DEG_PATH}'/Adrenal1vsAdrenal2_DEG_genes_zscores.tsv'
    OUT_PATH='../../data/rat_data/out_data/'

    dds_analysis plot_mr_vs_exp -inGeneEXPfile ${gene_exp_file}  \
          -dpi 300 -inMRfolder ${OUT_PATH}/out4dmr_in_deg_tss_5dist \
          -sampleName sample_name4replace.tsv -expTAB -inGene 'Tab2' -inMR 'chr1:mr16' -wtStr '_Ctrl' -output_path ${OUT_PATH}

.. image:: Tab2_chr1_mr16.jpg
   :alt: chr1:mr16 vs tab2