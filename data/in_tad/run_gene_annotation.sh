#!/bin/bash
#this script is used to generate a gene position file without remove the TSS and TES part of regions .
dmr_analysis dmr_gene_annotation -i no -l 10 -xL 50000000 -X 1 -Y 1 \
	-hu yes -n no -r ../genome/hg19/hg19.refFlat.txt -g ../genome/hg19/hg19.chrom.sizes.clear.sorted

