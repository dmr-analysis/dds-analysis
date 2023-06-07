
filterDEG4bpb3
==============


.. contents::
    :local:


This command is used to filter Differential Expressed Genes (DEG) by rratio based on the exported file from bpb3 differential_expression. It exports the filtered DEG with group mean and rratio.

Optional arguments:

- `-h, --help`: Show the help message and exit.

Required arguments for filtering:

- `--in_group1_str NAME`: Group1 column label.

- `--in_group2_str NAME`: Group2 column label.

- `--in_folder FILE_FOLDER`: Input file folder that contains the bpb3 exported differential expression gene file.

- `--in_file FILE_NAME`: Input file name of the bpb3 exported differential expression gene file.

Optional arguments for filtering (with default values):

- `--min_median_RPKM NUMBER`: Minimum median of RPKM value in each group. Default value is 1.0. If it is smaller than 0, this option will be switched off.

- `--rr_cutoff NUMBER`: Cutoff value for rratios to filter differential expression genes between two groups. Default value is 0.18. rr>0.14, 0.18, 0.22, 0.4, 0.67 represent 1.15, 1.2, 1.25, 1.5, 2.0 folder changes, respectively.

- `--sort_by_column_str NAME`: Column name of the dataframe column that will be sorted after filtering. Default value is differential_expressoin_T_test_p_value if the file is obtained from the bpb3 package.

- `--is_median`: Use mean or median of group values to calculate folder changes between groups. Default value is False, which uses the mean value of the group.

Please replace the angle-bracketed placeholders (`<>`) with appropriate values when using the command.