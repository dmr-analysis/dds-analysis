
map_block2chromSegment
======================


.. contents::
    :local:



Map mutation blocks to predicted chromatin segment regions.

Optional Arguments:

- `-h, --help`: Show this help message and exit

Required Arguments:

- `-in_BSpos IN_SORTEDBLOCK_FILE, --in_sortedBlock_file IN_SORTEDBLOCK_FILE`: Bpb3 exported mutation block summary position file in BED format and sorted by chromosome and position
- `-out_folder OUT_FILE_FOLDER, --out_file_folder OUT_FILE_FOLDER`: Export file folder
- `-in_chrSeg IN_CHROMSEGMENT_FOLDER, --in_chromSegment_folder IN_CHROMSEGMENT_FOLDER`: File folder path of chromatin segment files

Optional Arguments with Default Values:

- `-in_prefStr, --in_prefix_string4chromSegment_file`: Prefix string for chromatin segment file names. Default=combined*merged.bed.gz
- `-is_change, --change_name`: Whether to change the 4th column of name by reducing its length. Columns 1, 2, 3 are chr, start_pos, and end_pos, but column 4 is the id information. Default=False, does not change id information. If this option is used, the program will split the name by ":" and only take the first two elements for exporting a new id name.
- `-in_cutoff IN_DMR_MINIMUM_CUTOFF, --in_dmr_minimum_cutoff IN_DMR_MINIMUM_CUTOFF`: Minimum cutoff values for selecting blocks or mrs to export. Default=None, which means no minimum cutoff for exporting data or there is no 9th column (p-value or other values can be used to filter by minimum cutoff) in the dataframe.
- `-isMorB IS_MR_OR_BLOCKS, --is_MR_or_Blocks IS_MR_OR_BLOCKS`: Is input bed position files are MR (methylation regions chr#:mr#) or mutation blocks block_#_chr_start_end. 0 for MR, 1 for Blocks. Default=0 is MR.