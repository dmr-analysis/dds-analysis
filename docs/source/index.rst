.. image:: oslo-universitetssykehus-ous-logo.png

.. image:: uio_segl_a.png

.. image:: helse-sor-ost.jpg

.. image:: ahhus.png

=================
DDS-Analysis Tool
=================
An integrated data analysis pipeline by considering both Differential Methylation Region (DMR) and Differentially Expressed Genes (DEG) in SNP analysis. This package depends on some functions from the bpb3 package!

**positional arguments**:
- **task**: Pipeline task to run

**optional arguments**:
- **-h, --help**: Show this help message and exit
Usage: dds_analysis <task> [<args>]

Tasks available for using:
--------------------------------

- **bpb3summary2bed_format**: Convert bpb3 block summary file to a bed format file
- **map_block2genome**: Map mutation block to genomic regions
- **map_block2chromSegment**: Map mutation block to chromatin Segment regions
- **map_block2dmr**: Map mutation block to differential methylated regions
- **find_geneExp4block**: Find differential expressed genes for mutation blocks
- **find_block_patieintID**: Find patient ID for mutation blocks
- **combine_dmr_deg2block**: Combine DMR, DEG, and mutation block information together
- **filter_blocks**: Filter mutation blocks by using DMR or/and DEG condition
- **collect_gene_names4blocks**: Collect unique gene names for mutation blocks with DMR and/or DEG
- **check_block_gene_inTAD**: Check whether block and gene are in the same TAD or TAD boundary
- **dds_geneRanking**: Select top-ranked genes from final prediction
- **go_pathway_analysis4out_blocks_gene**: GO pathway analysis of genes
- **find_enhancer_target_genes**: Find enhancer and its target genes overlapping with mutation blocks that associated with selected gene
- **chromSegment_test4blocks**: Enrichment test of mutation blocks or methylation regions that associated with genes in 7 chromatin segmentations of the human genome
- **dTarget_methy_vs_express**: Predict long-distance target gene for a specific region (e.g., mutation block) based on coupling of methylation and gene expression across samples
- **plot_mr_vs_exp**: Plot DMR/MR methylation level and Gene expression for a pair of DMR and its target gene
- **plot_tss_enhancer_mrs**: Plot the average methylation level of predicted DMRs at TSS and enhancer regions by the target genes predicted from dTarget_methy_vs_express
- **filterDEG4bpb3**: Filter Differential Expressed Genes (DEG) by rratio based on the exported file from bpb3 differential_expression then export it with group mean and rratio
- **preprocess**: This module first finds DEG in TSS, 5dist regions then preprocesses data for dds_analysis.


Contents:
---------
.. toctree::
   :maxdepth: 1

   content/installation
   content/bpb3summary2bed_format
   content/map_block2genome
   content/map_block2chromSegment
   content/map_block2dmr
   content/find_geneExp4block
   content/find_block_patientID
   content/combine_dmr_deg2block
   content/filter_blocks
   content/collect_gene_names4blocks
   content/check_block_gene_inTAD
   content/dds_geneRanking
   content/go_pathway_analysis4out_blocks_gene
   content/find_enhancer_target_genes
   content/chromSegment_test4blocks
   content/dTarget_methy_vs_express
   content/plot_mr_vs_exp
   content/plot_tss_enhancer_mrs
   content/filterDEG4bpb3
   content/preprocess
   content/demoFL
   content/demoHAP1
   content/demoRAT
   content/demoFLMR


