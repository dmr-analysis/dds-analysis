#!/bin/bash
bpb3 mussd --patient_mutations /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/bpb3/final_demo_data/demo3_fl_cohort_small/out/../patient_data/DO52689/icgc_mutations_on_18.tsv /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/bpb3/final_demo_data/demo3_fl_cohort_small/out/../patient_data/DO27859/icgc_mutations_on_18.tsv /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/bpb3/final_demo_data/demo3_fl_cohort_small/out/../patient_data/DO52656/icgc_mutations_on_18.tsv /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/bpb3/final_demo_data/demo3_fl_cohort_small/out/../patient_data/DO52680/icgc_mutations_on_18.tsv /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/bpb3/final_demo_data/demo3_fl_cohort_small/out/../patient_data/DO52677/icgc_mutations_on_18.tsv /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/bpb3/final_demo_data/demo3_fl_cohort_small/out/../patient_data/DO52694/icgc_mutations_on_18.tsv /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/bpb3/final_demo_data/demo3_fl_cohort_small/out/../patient_data/DO52669/icgc_mutations_on_18.tsv /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/bpb3/final_demo_data/demo3_fl_cohort_small/out/../patient_data/DO52685/icgc_mutations_on_18.tsv /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/bpb3/final_demo_data/demo3_fl_cohort_small/out/../patient_data/DO52653/icgc_mutations_on_18.tsv /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/bpb3/final_demo_data/demo3_fl_cohort_small/out/../patient_data/DO52676/icgc_mutations_on_18.tsv /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/bpb3/final_demo_data/demo3_fl_cohort_small/out/../patient_data/DO52681/icgc_mutations_on_18.tsv \
	 --result_folder /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/bpb3/final_demo_data/demo3_fl_cohort_small/out/mussd_blocks \
	 --genome /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/bpb3/final_demo_data/genome_data/hg19/human_g1k_v37_decoy.fasta \
	 --min_patients_in_block 3 \
	 --block_flank 25 \
	 --min_block_size 3 \
	 --cluster_distance 30 \
	 --block_distance 500\
	 --regions /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/bpb3/final_demo_data/demo3_fl_cohort_small/out/regions.bed