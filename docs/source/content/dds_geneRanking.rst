
dds_geneRanking
=================


.. contents::
    :local:


Select top-ranked genes and blocks from the final prediction.

Optional Arguments:

- `-h, --help`: Show the help message and exit.

Required Arguments:

- `-inUqGene IN_UNIQUE_GENE_FILE, --in_unique_gene_file IN_UNIQUE_GENE_FILE`: Input a file containing a list of unique genes exported from dds_analysis. This file should contain information about the unique genes identified in the analysis.

- `-inDEGFile IN_DEG_FILE, --in_DEG_file IN_DEG_FILE`: Input file containing a list of differentially expressed genes exported from BayesPI-BAR2. This file should contain information about the differentially expressed genes.

- `-inDMRFile IN_DMR_FILE, --in_DMR_file IN_DMR_FILE`: Input file containing differential methylation region information exported from dmr_analysis. This file should contain information about the differential methylation regions.

Optional Arguments with Default Values:

- `-inCutoff, --in_cutoff_pval4score`: Cutoff probability value for selecting top-ranked genes. Default: 0.5.

- `-TSS, --TSS_score`: Weight score for TSS (Transcription Start Site). Default: 4.0.

- `-gene, --gene_score`: Weight score for gene. Default: 2.0.

- `-TES, --TES_score`: Weight score for TES (Transcription End Site). Default: 2.0.

- `-dist, --dist_score`: Weight score for 5 distance regions. Default: 1.0.

- `-enhancer, --enhancer_score`: Weight score for enhancer. Default: 3.0.