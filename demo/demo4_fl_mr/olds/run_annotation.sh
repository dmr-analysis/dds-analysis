#!/bin/bash
#this script is used to generate a gene position file without remove the TSS and TES part of regions .
dmr_analysis dmr_gene_annotation -i no -l 10 -xL 50000000 -X 10000 -Y 10000 \
        -rem no -hu yes -n no -M 10000 -N 1000000 \
	-r ../../data/final_demo_data/genome/hg19/gencode.v19.annotation_gene.refFlat \
	-g ../../data/final_demo_data/genome/hg19/hg19.chrom.sizes.clear.sorted \
	--folderOut ../../data/fl_mr/in_data/in_genome_regions/


