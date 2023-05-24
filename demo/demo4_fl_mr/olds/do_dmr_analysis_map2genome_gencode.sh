#!/bin/bash
# path of main folder to output data
OUT_FOLDER='../../data/fl_mr/out_data/'
IN_FOLDER='../../data/fl_mr/in_data/'

#make MRs to gencode annotation files
IN_GENOME_REGIONS=${IN_FOLDER}'/in_genome_regions/'
in_sorted_refFlat='gencode.v19.annotation_gene.ref_clean_sorted.bed'
logProb_cutoff=0.6


in_mr_result_folder=${IN_FOLDER}/'final_demo_data/fl_12samples/out_data/DMR_CpG_context/'
out_result_folder=${OUT_FOLDER}/'out_data_genecode/DMR_CpG_context2'
mr_IN_FILE='2_chroms_all_mr_data_range_dmrRanking'
out_folder4genome_map='out_map2genome'
out_file4genome_map=gcb_vs_tumor_DMR_hyper_hypo_mix_${logProb_cutoff}.csv


#b) map DMR to predefined genomic regions such as TSS, TES, gene et al.
dmr_analysis dmr_map2genome --in_sortedDMR_file ${in_mr_result_folder}/${mr_IN_FILE}.bed \
	--in_geneRegion_file ${IN_GENOME_REGIONS}/list_region_files.txt_gencode \
	--in_outFile_folder ${out_result_folder}/${out_folder4genome_map} \
	--in_refFlat_file ${IN_GENOME_REGIONS}/data/${in_sorted_refFlat}
echo dmr_map2genome - Done

#c) calculate percentage of DMR in annotated genomic regions
dmr_analysis dmr_cal2genome_percent --in_outFile_folder ${out_result_folder}/${out_folder4genome_map} \
        --in_outFile_name ${out_file4genome_map} --in_LogReg_proba ${logProb_cutoff} \
        --in_fileName_string $mr_IN_FILE  
echo dmr_cal2genome_percent - Done

#d) plot percentage of DMR in annotated genomic regions
dmr_analysis dmr_percent2plot --in_countFile_folder ${out_result_folder}/${out_folder4genome_map} \
        --in_countFile_name ${out_file4genome_map}
echo dmr_percent2plot - Done

