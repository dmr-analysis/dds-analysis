# 1. make genome region files based on genecode gtf annotation file
./run_annotation.sh

#2. do dds for genCode genome regions
./do_dds_analysis_fl_14samples_mr_enhancers_gencode.sh

#3. map MR to gencode annotation regions
./do_dmr_analysis_map2genome_gencode.sh

#4. prepare files before doing dmr vs deg analysis 
#./run_fl_mr_part1.sh

#5. export MR files and make background file list before doing mr vs deg prediction
#./run_fl_mr_part2.sh

#4. prerpocess before dTarget
./run_fl_preprocess.sh

#5. if shell background list creater is failed then use the python code
#to create the background file list
#make_background_files.py

#6. do dds_analysis of mr vs deg prediction
./run_fl_mr_part3.sh


