plot_tss_enhancer_mrs
=====================


.. contents::
    :local:


This module plots the average DMR methylation levels at TSS and enhancer regions based on their predicted target genes from dTarget_methy_vs_express. It also exports differential gene and differential methylation information for each pair of putative targets and TSS or enhancers.

Optional arguments:

- `-h, --help`: Show the help message and exit.

Required arguments:

- `-exp_file DIFFERENTIAL_GENE_EXPRESSION_FILE, --differential_gene_expression_file`: A tab-delimited file with at least four columns (e.g., gene, group1_mean, group2_mean, rratio/log2folder change) containing differential gene expression information, which is exported by (bpb3 filterDEG4bpb3). If more than four columns are included in this file, the program will assume the first column is gene, and the last 3 columns are group1_mean, group2_mean, and rratio/log2folder change, respectively.

- `-genes SELECTED_GENES4PLOT, --selected_genes4plot SELECTED_GENES4PLOT`: A string of gene symbols that will be selected for plotting the average methylation levels in their predicted regulatory regions such as TSS and Enhancer. Gene symbols are separated by commas (,) such as gene1,gene2,gene3.

Optional arguments:

- `-dmr_file DIFFERENTIAL_METHYLATION_FILE, --differential_methylation_file`: A tab-delimited file with only five columns (e.g., chrom, start_pos, end_pos, infos, pval) containing differential methylation regions, which is exported by (dmr_analysis dmr_combine_multChrs4rank). For the column "infos," the data structure is chrom:mr_id:anything:hyper/hypo:anything, where the information in the first, second, and fourth strings separated by ":" will be used by the program for finding DMRs.

- `-tss_file TSS_TARGET_FILE, --tss_target_file TSS_TARGET_FILE`: A tab-delimited file with nine columns (index, gene, mr_id, rg_Tvalue, rg_Pvalue, enhancer, gene_type, utest_Pvalue, expectPvalue) containing predicted information of TSS target genes by using (dds_analysis dTarget_methy_vs_express).

- `-enc_file ENHANCER_TARGET_FILE, --enhancer_target_file ENHANCER_TARGET_FILE`: A tab-delimited file with nine columns (index, gene, mr_id, rg_Tvalue, rg_Pvalue, enhancer, gene_type, utest_Pvalue, expectPvalue) containing predicted information of enhancer target genes by using (dds_analysis dTarget_methy_vs_express).

- `-mr_folder MR_DATA_FOLDER, --mr_data_folder MR_DATA_FOLDER`: A file folder containing all DNA methylation data of methylated regions exported by dmr_analysis dmr_exportData.

- `-folder_name PREFIX_FOLDERNAME4CHROM, --prefix_foldername4chrom PREFIX_FOLDERNAME4CHROM`: A prefix folder name for exported methylation data of each chromosome under -mr_folder (e.g., prefix_foldername4chrom + chr1 means all methylation data for chr1).

- `-out_folder, --out_data_folder`: Output folder name for exported data from the current analysis. If this folder does not exist in the input data folder, it will be created. Default is "data".

- `-out_plots, --out_plot_folder`: Output folder name for all plots. Default is "plots".

- `-is_negative, --is_negative_target`: Select plot for negative, positive, or both responses of predicted DMR targets. Use 0, 1, or 2

for positive, negative, and all responses, respectively. Default is 1 for negative targets.

- `-gX, --gene_upstream_X`: Set the length of upstream to transcription start side (TSS) in the plot. Default is 500.

- `-gY, --gene_downstream_Y`: Set the length of downstream to transcription start side (TSS) in the plot. Default is 500.

- `-G, --gene_or_enhancer_length`: Set the gene or enhancer length in the plot. Default is 2000.

- `-w, --step_size_in_region`: Step size in a region for data smoothing. Default is 100.

- `-window, --window_size_in_combined_regions`: Window size of combined regions such as TSS and Enhancer. Default is 600. Increase this parameter if there is a gap between two combined regions (>500).

- `-sigma, --data_smooth_sigma`: Data smooth parameter sigma. Default is 50.

- `-flank, --flank_region`: Flank region to add on the two sides of a promoter or enhancer region. Default is 100.

- `-MRs_txt, --MRs_txt_description`: Text description of MRs. Default is "DMRs".

- `-dmr_compressed, --dmr_file_not_compressed`: Specify this option if the DMR file is not compressed. Default is False (dmr_file is compressed).

- `-flip_rr, --not_flip_sign_of_rratio`: Specify this option if you do not want to flip the sign of rratio. Default is False (sign of rratio is flipped).

Please replace the angle-bracketed placeholders (`<>`) with appropriate values when using the command.