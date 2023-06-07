filter_blocks
=================


.. contents::
    :local:


Filter mutation blocks by using DMR and/or DEG information.

Optional Arguments:

- `-h, --help`: Show the help message and exit.

Required Arguments:

- `-in_BKfile IN_COMBINED_DMRDEGBLOCK_FILE, --in_combined_DmrDegBlock_file IN_COMBINED_DMRDEGBLOCK_FILE`: A bpb3 mutation block file combined with DMR, DEG, and genomic region information that is exported by dds_analysis combine_dmr_deg2block. Here, we select blocks that satisfy condition 1: blocks associated with both DMR (not null) and DEG/enhancer (not null), and condition 2: blocks associated with either DMR (not null) or DEG/enhancer (not null). DMR and DEG/enhancer are already defined in previous analysis. If not, they are marked as "nan" in the input file.

Optional Arguments with Default Values:

- `--not_use_enhancer`: Whether to use/consider enhancer information during the filtering. The default value is True. If this option is used, enhancer information will not be used during the filtering process.