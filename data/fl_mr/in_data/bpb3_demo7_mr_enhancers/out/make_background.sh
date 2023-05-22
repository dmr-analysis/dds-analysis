#!/bin/bash
        #SBATCH --job-name=bg_master
#SBATCH --time=200:00:00
#SBATCH --mem-per-cpu=1500M --partition=bigmem 
#SBATCH --cpus-per-task=1


bpb3 background_affinity_changes  \
	 --block_size 768 \
	 --result_folder /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/to_github/bpb3_git/final_demo_data/demo7_mr_enhancers/out/background \
	 --genome /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/to_github/bpb3_git/final_demo_data/genome_data/hg19/human_g1k_v37_decoy.fasta \
	 --regions ../../final_demo_data/demo7_mr_enhancers/in/dmr_intersect_enhancer_regions_chr18noChr_5kb_flank.tsv \
	 --mutations_distribution 35 21 13 8 11 9 8 5 8 24 32 10 22 5 22 24 6 15 8 4 9 1 \
	 --max_rank 10 \
	 \
	 --pwm_folder /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/to_github/bpb3_git/final_demo_data/demo7_mr_enhancers/out/pwm_for_background \
	 --chemical_potentials none -10 -13 -15 -18 -20 \
	 --iterations 1000 \
	 --p_value_cutoff 0.1 \
	 \
	 --integration mean_ddba \
	 --block_resample_iterations 10 \
	 --block_samples_to_take 100 \
	 --use_cores 10 \
	 --max_nodes 10 \
	 --reuse_output \
	 --background_mutations /cluster/projects/nn4605k/junbaiw/rrpathology/local/FL_project2019/bayespi_bar_2_alireza/bayes_pi_bar_3.1/to_github/bpb3_git/final_demo_data/demo7_mr_enhancers/out/mussd_blocks/mutations_in_regions.tsv
