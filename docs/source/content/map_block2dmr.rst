map_block2dmr
=================


.. contents::
    :local:


:map_block2dmr:

Map mutation blocks to differential methylated regions.

Optional Arguments:

- `-h, --help`: Show this help message and exit

Required Arguments:

- `-inBKfile IN_SORTEDBLOCK_FILE, --in_sortedBlock_file IN_SORTEDBLOCK_FILE`: Bpb3 exported mutation block summary position file in BED format and sorted by chromosome and position
- `-inDMRfile IN_DMR_FILE, --in_dmr_file IN_DMR_FILE`: DMRs of all chromosomes in BED format that ranked exported by dmr_analysis dmr_combine_multChrs4rank

Optional Arguments with Default Values:

- `-inFlank, --in_flank_region2block`: Add N bp flank regions on both sides of the blocks. Default=0, do not add flank regions.
- `-outFolder, --out_file_folder`: Path of output file folder. Default=out_blocks_dmr_single/
- `-min_cutoff, --dmr_min_cutoff`: Minimum cutoff value for selecting DMR regions from dmr_analysis exported file. Default=0.6.