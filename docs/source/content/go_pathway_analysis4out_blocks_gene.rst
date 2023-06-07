go_pathway_analysis4out_blocks_gene
===================================


.. contents::
    :local:


Perform GO pathway analysis for genes associated with mutation blocks.

Optional Arguments:

- `-h, --help`: Show the help message and exit.

Required Arguments:

- `-in_GEfile IN_GENE_FILE, --in_gene_file IN_GENE_FILE`: Input file containing unique gene information after combining DMR, DEG, TAD, and mutation block information. This file should contain information about the genes associated with the mutation blocks.

- `-in_GOfile IN_GO_FILE, --in_GO_file IN_GO_FILE`: Input file containing the DAVID export file. This file should contain the exported data from the DAVID (Database for Annotation, Visualization, and Integrated Discovery) tool for GO (Gene Ontology) pathway analysis.

Optional Arguments with Default Values:

- `-cutoff, --cutoff_Pvalue4GO`: Cutoff p-value for GO enrichment terms. Default: 0.05.