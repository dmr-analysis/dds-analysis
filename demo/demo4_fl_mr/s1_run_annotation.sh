#!/bin/bash

#generate predefined genomic regions files by using dmr_analysis based on genCode annontation 
OUT_FOLDER='../../data/fl_mr/out_data/in_genome_regions/'
IN_FOLDER='../../data/fl_mr/in_data/gencode/'

dmr_analysis dmr_gene_annotation -F ${OUT_FOLDER} -i no -l 10 \
        -rem no -xL 50000000 -X 1 -Y 1 -M 10000 -N 1000000 -hu yes -n no \
        -r ${IN_FOLDER}/gencode.v19.annotation_gene.refFlat \
        -g ${IN_FOLDER}/hg19.chrom.sizes.clear.sorted
echo gene_annotation-DON
#here enhancer file has to be added manually in list_region_files.txt



