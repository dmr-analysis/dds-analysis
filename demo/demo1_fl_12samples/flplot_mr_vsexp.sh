#!/bin/sh


echo ""
dds_analysis plot_mr_vs_exp -inGeneEXPfile '../../data/fl_12samples/out_data/demo3_fl_cohort_small/out/differentially_expressed_genes.txt'  \
        -dpi 300 -inMRfolder '../../data/fl_12samples/out_data/DMR_CpG_context/' \
	 -inGene 'BCL2' -inMR 'chr18:mr621' -wtStr 'counts' -output_path './olds/' -pathDepth 16
echo "Done with plot_vs_exp "
echo ""
