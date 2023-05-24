#!/bin/bash
#Before running following steps, it assumes taht DMRs are already predicted by dmr_abalysis
#
#the script is used to prepare files for dds_analysis

python find_DEG_in_tss_5dist_regions_fl_mr.py
echo "To find DMR regions that are overlapping with TSS or 5distance regions of DEG - Done"
 
