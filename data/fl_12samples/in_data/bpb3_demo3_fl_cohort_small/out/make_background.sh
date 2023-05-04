#!/bin/bash
        #SBATCH --job-name=bg_master
#SBATCH --time=200:00:00
#SBATCH --mem-per-cpu=1500M --partition=bigmem 
#SBATCH --cpus-per-task=1


bpb3 background_affinity_changes  \
	 --block_size 611 \
	 --result_folder /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/bpb3/final_demo_data/demo3_fl_cohort_small/out/background \
	 --genome /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/bpb3/final_demo_data/genome_data/hg19/human_g1k_v37_decoy.fasta \
	 --regions /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/bpb3/final_demo_data/demo3_fl_cohort_small/out/regions.bed \
	 --mutations_distribution 6 4 11 7 7 5 18 13 31 4 3 10 9 6 9 5 13 \
	 --max_rank 15 \
	 \
	 --pwm_folder /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/bpb3/final_demo_data/demo3_fl_cohort_small/out/pwm_for_background \
	 --chemical_potentials none -10 -13 -15 -18 -20 \
	 --iterations 1000 \
	 \
	 --integration mean_ddba \
	 --block_resample_iterations 10 \
	 --block_samples_to_take 8 \
	 --use_cores 10 \
	 --max_nodes 10 \
	 --reuse_output \
	 --background_mutations /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/bpb3/final_demo_data/demo3_fl_cohort_small/out/mussd_blocks/mutations_in_regions.tsv
