bpb3summary2bed_format
======================


.. contents::
    :local:

Usage: `dds_analysis bpb3summary2bed_format`

Convert bpb3 block summary file to a bed format file, where the fourth column is block ID, only works for human genome.

Optional arguments:

- `-h, --help`: Show this help message and exit

Required:

- `-in_block IN_BLOCK_SUMMARY_FILE, --in_block_summary_file IN_BLOCK_SUMMARY_FILE`: Input block summary file that is exported by bpb3 mussd analysis

Optional, with default values:

- `-out_path OUTPUT_FILE_PATH, --output_file_path OUTPUT_FILE_PATH`: Output file path for exported new bed files, default=./