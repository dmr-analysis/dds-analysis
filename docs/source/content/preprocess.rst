
preprocess
==========


.. contents::
    :local:


This module finds DEG in TSS and 5dist regions and preprocesses the data for dds_analysis.

Optional arguments:

- `-h, --help`: Show the help message and exit.

Required arguments:

- `-in_folder IN_FOLDER, --in_folder IN_FOLDER`: Input file folder for MRs overlapping with TSS and 5Distance regions that was done in dmr_analysis.

- `-in_string IN_STRING, --in_string IN_STRING`: String that will appear in the output file of DMR regions overlapped with TSS or 5dist.

- `-in_tss_file_mr IN_TSS_FILE_MR, --in_tss_file_mr IN_TSS_FILE_MR`: Path of a file that has DMR overlapped with TSS.

- `-in_dist_file IN_DIST_FILE, --in_dist_file IN_DIST_FILE`: Path of a file that has DMR overlapped with 5dist.

- `-in_deg_file IN_DEG_FILE, --in_deg_file IN_DEG_FILE`: Path of a file that has DEG and is exported using bpb3.

- `-out_folder OUT_FOLDER, --out_folder OUT_FOLDER`: Output folder path.

- `-gene_col_name GENE_COL_NAME, --gene_col_name GENE_COL_NAME`: Name of the gene column. See your DEG file and its header to find the column name.

- `-tss_file TSS_FILE, --tss_file TSS_FILE`: Path of a file that has TSS regions.

- `-full_mr_file FULL_MR_FILE, --full_mr_file FULL_MR_FILE`: File path of all ranked DMRs in the result folder of dmr_analysis.

- `-in_genome_file IN_GENOME_FILE, --in_genome_file IN_GENOME_FILE`: Genome file path which contains refFlat file and enhancer files.

