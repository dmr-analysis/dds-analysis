
plot_mr_vs_exp
==============


.. contents::
    :local:



This command is used to plot DMR/MR methylation levels and gene expression for a pair of DMR and its target gene.

Optional arguments:

- `-h, --help`: Show the help message and exit.

Required arguments:

- `-inGeneEXPfile IN_GENE_EXPRESSION_FILE_FROM_BPB3, --in_gene_expression_file_from_bpb3`: Input file name for gene expression profiles exported by bpb3 differential_expression. The file format should be: gene_name; P-values, sample1, sample2, ...

- `-inMRfolder IN_METHYLATION_FOLDER_FROM_DMR, --in_methylation_folder_from_dmr`: Input file folder path for methylation data. The data should be extracted by dmr_analysis dmr_exportData based on predefined blocks or regions.

Optional arguments with default values:

- `-inGene, --in_gene_name`: Target gene name for an MR/DMR region. Default value is None. Specify the target gene name for the MR/DMR region.

- `-inMR, --in_mr_id`: MR/DMR ID of a target gene, such as chr1:mr1234. Default value is None. Specify the MR/DMR ID of the target gene.

- `-inFile, --in_gene_mr_file`: Input file that contains both DMR/MR ID and the target gene name. It should be a tab-delimited text file, where the column labels are mr_id and gene, respectively. The file format should be the same as the export file from dTarget_methy_vs_express.py. Default value is None. Specify the input file that contains both DMR/MR ID and the target gene name.

- `-expTAB, --is_tab_delimated_file4diffExp`: Specifies whether the input differential expression gene file is in the bpb3 differential_expression exported file format or a common tab-delimited file format. If set, the input gene expression file will be treated as a common tab-delimited file. Default value is False.

- `-expStartCol, --gene_expression_start_column`: Start column number (e.g., 0, 1, 2, ...) for samples in the gene expression profile file. Default value is 2 for bpb3 exported file. Specify the start column number for samples in the gene expression profile file.

- `-pathDepth IN_DIFFEXP_PATHDEPTH4SAMPLENAME, --in_diffExp_pathDepth4sampleName`: File path depth for the column sample name. For example, if the column name is /path/folder/sample_name/.., then the depth of the sample name is the number of "/" before the sample name. Default value is 12. Specify the file path depth for the column sample name.

- `-sampleName IN_SAMPLENAME_FOR_REPLACE, --in_sampleName_for_replace IN_SAMPLENAME_FOR_REPLACE`: Input file name for renaming the old sample names in the gene expression data. The first column (old_name) contains the old sample names, and the second column (new_name) contains the new names used for replacement. The default file name is sample_name4replace.tsv, which can be an empty file if no replacement is needed for sample names.

- `-dpi, --figure_resolution_dpi`: Exported figure resolution in dpi. Default dpi is 60. Set the exported figure resolution in dpi.

- `-wtStr, --wildType_fileString`: The first few character string that labels the file name as the Wild Type condition. For example, if the file name starts with gcb_meht1_* as the wild type control sample, then --wildType_fileString is gcb. This is the default setting in the program.

- `-output_path, --output_file_path`:

 Output path for exported figures. Default value is "./" (current directory). Specify the output path for exported figures.