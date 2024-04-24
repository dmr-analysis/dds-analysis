
dTarget_methy_vs_express
========================


.. contents::
    :local:


Predict long distance target genes for a specific region (e.g., mutation block) based on the coupling of methylation and gene expression across samples.

Optional Arguments:

- `-h, --help`: Show the help message and exit.

Required Arguments:

- `-inGeneMRfile IN_GENE_MR_FILE_FROM_DDS, --in_gene_mr_file_from_dds IN_GENE_MR_FILE_FROM_DDS`: Input file name for gene-MR relationships exported by dds_analysis check_block_gene_inTAD. The file format should be: gene_name; gene_type; block_id; new_mr_sites; patients; isTAD; enhancers; patient_id.

- `-inGeneEXPfile IN_GENE_EXPRESSION_FILE_FROM_BPB3, --in_gene_expression_file_from_bpb3`: Input file name for gene expression profiles exported by bpb3 differential_expression. The file format should be: gene_name; P-values, sample1, sample2, ...

- `-inMRfolder IN_METHYLATION_FOLDER_FROM_DMR, --in_methylation_folder_from_dmr`: Input file folder path for methylation data, extracted by dmr_analysis dmr_exportData based on predefined blocks or regions.

- `-inBackgroundList IN_MR_BACKGROUND_FILE_LIST, --in_mr_background_file_list`: A list of files to be used in background sampling, where each file name is listed on a separate line.

Optional Arguments with Default Values:

- `-cutoff, --expected_pval_cutoff`: P-value cutoff for the expected P-value after sampling. Default: 0.05. Set the desired cutoff value for the expected P-value.

- `-reg_cutoff, --regression_pval_cutoff`: P-value cutoff for the regression model when fitting expression data against methylation data. Default: 0.05. Set the desired cutoff value for the regression model P-value.

- `-totalSamples, --total_samples_in_random_selection`: Total number of samples to draw from random background samples. Default: 1000. Specify the desired number of samples to draw.

- `-numOfprocesses, --number_of_parallel_processes`: Number of parallel processes to be used in the calculation. Default: 10. Specify the number of parallel processes to be used.

- `-expStartCol, --gene_expression_start_column`: Start column number (e.g., 0, 1, 2, ...) for samples in the gene expression profile file. Default: 2 for bpb3 exported file. Specify the start column number for samples in the gene expression profile file.

- `-sampleName IN_SAMPLENAME_FOR_REPLACE, --in_sampleName_for_replace IN_SAMPLENAME_FOR_REPLACE`: Input file name for renaming the old sample names in the gene expression data. The first column is the old sample name, and the second column is the new name used for replacement. The default file name is sample_name4replace.tsv, which can be an empty file if no replacement is needed for sample names.

- `-pathDepth IN_DIFFEXP_PATHDEPTH4SAMPLENAME, --in_diffExp_pathDepth4sampleName`: File path depth for the column sample name. For example, if the column name is /path/folder/sample_name/..., then the depth of the sample name is 3, which is the number of "/" before the sample name. Default: 12. Specify the file path depth for the column sample name.

- `-expTAB, --is_tab_delimated_file4diffExp`: Whether the input differential expression gene file is in bpb3 differential_expression exported file format or a common tab-delimited file format (e.g., the first column is the gene name, and the rest of the columns are gene expression data in samples). Default: False. If this option is

 used, the input gene expression file will be treated as a common tab-delimited file.

- `-mrTAB, --is_tab_delimated_file4MR`: Whether the input MR regions files are in dds_analysis check_block_gene_inTAD exported format or a common tab-delimited file format (e.g., the first column is the gene_name, followed by the gene_type and new_mr_sites). Default: False. If this option is used, the input MR region file will be treated as a common tab-delimited file.

- `-no_dmr_format, --not_dmr_analysis_format`: The input MR data folder is in dmr_analysis export data format. Default: False. If this option is used, the data folder format is not in the dmr_analysis export format.

- `-outName OUTPUT_FILE_NAME, --output_file_name OUTPUT_FILE_NAME`: Prefix file name for the output data exported by this program. Default: output_result_of_*. Specify the prefix file name for the output data.

- `-mrColName MR_ID_COLUMN_NAME, --mr_id_column_name MR_ID_COLUMN_NAME`: Column name of MR IDs that contain strings of joined MR IDs (e.g., chr1:mr1~chr1:mr11). Default: new_mr_sites. If the MR IDs are in another column, please provide the new column name.

- `-output_path OUTPUT_FILES_PATH, --output_files_path OUTPUT_FILES_PATH`: Output files path for data exported by the program. Default: ./, which exports files in the local directory. Specify the output files path for the exported data.