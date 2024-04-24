
combine_dmr_deg2block
=====================


.. contents::
    :local:


Combine information from DMR, DEG, and mutation blocks.

Optional Arguments:

- `-h, --help`: Show the help message and exit.

Required Arguments:

- `-in_SBfile IN_SORTEDBLOCK_PATIENT_FILE, --in_sortedBlock_patient_file IN_SORTEDBLOCK_PATIENT_FILE`: Bpb3 exported mutation block summary position file in BED format and sorted by chromosome and position, where patient IDs are added by dds_analysis find_block_patientID.
- `-in_DMRfile IN_DMR_FILE, --in_dmr_file IN_DMR_FILE`: bpb3 block summary file with DMR information that is exported by dds_analysis map_block2dmr.
- `-in_DEGfiles IN_DEG_FOLDER_AND_FILE_SUFFIX, --in_deg_folder_and_file_suffix`: Path of a file folder that contains results exported by dds_analysis find_geneExp4block, where DEG information is added to blocks (e.g., file_path\/\*.tsv).

Optional Arguments with Default Values:

- `-outFold, --out_file_folder`: Path of output file folder.
