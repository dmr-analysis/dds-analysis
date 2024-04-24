
chromSegment_test4blocks
========================


.. contents::
    :local:


Perform an enrichment test of mutation blocks or methylation regions associated with genes in the 7 chromatin segmentations of the human genome.

Optional Arguments:

- `-h, --help`: Show the help message and exit.

Required Arguments:

- `-inGeneFile IN_GENE_FILE, --in_gene_file IN_GENE_FILE`: Input file name for the gene-block association file, exported by dds_analysis dds_geneRanking. The columns in this file should include [gene_name, gene_type, block_id, new_mr_sites, patients, isTAD, enhancers, patient_id, geneType_meanScore, num_of_patients, diffExp_pval, dmr_pval, median_scores, mean_scores]. Only the gene_name and block_id columns will be used for this analysis.

- `-inChromFold IN_CHROMSEG_FOLDER, --in_chromSeg_folder IN_CHROMSEG_FOLDER`: Input file folder name that contains the 7 types of chromatin state classifications for all blocks. The files should be in BED format with the following columns: [chrom, start_pos, end_pos, chromatin_type, block_chrom, block_start_pos, block_end_pos, block_id]. The chromatin_type can be one of the following: E, T, R, PF, TSS, WE, CTCF.

Optional Arguments with Default Values:

- `-samples, --number_of_samples`: Number of randomly drawn samples. Default: 100. Specify the desired number of samples for the enrichment test.

- `-process, --number_of_processes`: Number of parallel processes to be used in the calculation. Default: 20. Specify the number of processes to be used for parallel computation.

- `-cutoff, --cutoff_of_absolute_log10_pvalue`: Cutoff value for filtering, which is the absolute value of log10(expected P-values). Default: 1.3, which is equivalent to a P-value of 0.05. Set the desired cutoff value for filtering.

- `-dpi, --figure_resolution_dpi`: Exported figure resolution in DPI (dots per inch). Default: 60. Specify the resolution (in DPI) for the exported figures.