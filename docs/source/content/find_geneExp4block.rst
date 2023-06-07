find_geneExp4block
==================

.. contents::
    :local:

Find differential gene expression for mutation blocks.

Optional Arguments:

- `-h, --help`: Show this help message and exit.

Required Arguments:

- `-in_BGfolder IN_BLOCKS_GENOME_FOLDER, --in_blocks_genome_folder IN_BLOCKS_GENOME_FOLDER`: Path of blocks mapped to genomic region files, which are exported by map_block2genome.
- `-in_SBfile IN_SORTEDBLOCK_FILE, --in_sortedBlock_file IN_SORTEDBLOCK_FILE`: Bpb3 exported mutation block summary position file in BED format and sorted by chromosome and position.
- `-in_DGfile IN_DEG_FILE, --in_deg_file IN_DEG_FILE`: A differentially expressed gene file in tab delimited format, which is an export file from bpb3 differential_expression.
- `-out_folder OUT_FILE_FOLDER, --out_file_folder OUT_FILE_FOLDER`: Path of file folder to export results.

Optional Arguments with Default Values:

- `-in_FTlist, --in_feature_list`: A comma-separated string representing features of genomic regions that will be considered. Default=TSS,gene,TES,5dist,enhancers, where each feature is separated by a comma.
- `-in_cutoff IN_DMR_MINIMUM_CUTOFF, --in_dmr_minimum_cutoff IN_DMR_MINIMUM_CUTOFF`: Minimum cutoff values for selecting blocks or mrs to export. Default=None, which means no minimum cutoff for exporting data or there is no 9th column (p-value or other values can be used to filter by minimum cutoff) in the dataframe.