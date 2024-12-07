#this script is used to combine DMR DEG results to blocks
import pandas as pd
import glob
import os
import argparse
#from bpb3.script.script_high.other.common import check_folder
# from .script_high.functions_from_bpb3 import check_folder

#exec(open('combine_dmr_deg2block.py').read())
def my_parser(parser):
  required = parser.add_argument_group('Required')
  required.add_argument('-in_folder','--in_folder', required=True,
                  help='input file folder for MRs overlapping to TSS and 5Distance regions that was done in dmr_analysis')
  required.add_argument('-in_string','--in_string', required=True, help='string that will appear in output file of dmr region over lapped with tss or 5dist')
  required.add_argument('-in_tss_file_mr','--in_tss_file_mr',required=True,
            help='Path of a file that have dmr over lapped with tss ')
  required.add_argument('-in_dist_file', '--in_dist_file', required=True,
                        help='Path of a file that have dmr over lapped with 5 dist ')
  required.add_argument('-in_deg_file', '--in_deg_file', required=True,
                        help='Path of a file that have DEG and is exported using bpb3 ')
  required.add_argument('-out_folder', '--out_folder', required=True,
                        help='Path output folder ')
  required.add_argument('-gene_col_name', '--gene_col_name', required=True,
                        help='name of gene column, see your DEG file and its header to find column name.')

  required.add_argument('-tss_file', '--tss_file', required=True,
                        help='Path of a file that have tss region. ')
  required.add_argument('-full_mr_file', '--full_mr_file', required=True,
                        help='file path of all ranked DMRs in the result folder of dmr_analysis ')
  required.add_argument('-in_genome_file', '--in_genome_file', required=True,
                        help='genome file path which contains refFlat file and enhancer files ')


  return parser

def uq_of_column(in_df,col_name,sep_str='~'):
    '''Return a new df with a column with unique names/strings separsted by sep_str
    '''
    out_df=in_df.copy()
    out_df[col_name]=in_df[col_name].apply(lambda x: sep_str.join(list(set(x.split(sep_str)))))
    return out_df.copy()

def preprocess(full_mr_file, out_result_folder, in_genome_file, in_data_str, tss_file):
    ###preprocess file#######

    # this script is used to generate data for DMR vs. expressiondata


    # this script is used to prepare for files before running dds_Analysis dTarget_methy_vs_express



    # 4. to find mr not located in either enhacner or tss regions that will be used to build a list of background MRs
    # To build a list of MRs for randomly selected background MRs , we first need to remove MRs that are located in either TSS or enhancers regins


    out_file = os.path.join(out_result_folder, 'mr_regions_not_in_enhancers' + in_data_str + '.bed')
    os.system('bedtools intersect -a ' + full_mr_file + ' -b ' + in_genome_file + \
              ' -v > ' + out_file)
    out_file2mr_not_in_enhancer_tss = out_file.replace('.bed', '_tss.bed')
    os.system('bedtools intersect -a ' + out_file + ' -b ' + tss_file + \
              ' -v > ' + out_file2mr_not_in_enhancer_tss)
    os.system('rm -f ' + out_file)
    print('Export MR regions ddo not locate in TSS or Enhancer regions:')
    print('\t', out_file2mr_not_in_enhancer_tss, '\n')

    # 5. input DMR regions that overlapping with TSS and 5distance regions of DEGs that export by find_DEG_in_tss_5dist_regions_test.py"
    in_file = os.path.join(out_result_folder, 'dmr_regions_in_deg_tss_5dist' + in_data_str + '.bed')
    in_df = pd.read_csv(in_file, sep='\t', header=None)
    print('Read DMRs that are overlapping with either TSS or 5distance regions of DEGs : ')
    print(in_file, '\n')

    # 14a.
    # 6. to find unique dmr located in tss and 5dist regions of DEG
    # here assue columns 0,1,2,3 are chr, star, end, and name, and 5 is the genome type
    out_df = in_df.groupby([0, 1, 2, 3])[5].apply(';'.join).reset_index().copy()
    out_file = in_file.replace('dmr_', 'uqdmr_')
    out_df.to_csv(out_file, sep='\t', header=False, index=False)
    print('Export unique DMR regions that overlapping to TSS or 5distance regions of DEG: ')
    print('\t', out_file, '\n')

    # 7. to continue find unique gene dmrs in tss and 5dist regions of DEG
    # column 3 is DMR id, column 5 is the region/gene type such as TSS, 5dist et al
    in_df['mr_site'] = in_df[3].apply(lambda x: ':'.join(x.split(':')[0:2]))
    in_df['gene_type'] = in_df[5].apply(lambda x: x.split('||')[1])

    out_df2 = in_df.groupby([6, 'gene_type'])['mr_site'].apply('~'.join).reset_index().copy()
    out_df2.columns = ['gene_name', 'gene_type', 'new_mr_sites']
    out_file2 = in_file.replace('dmr_', 'uqGeneDmr_')

    #jbw may
    new_out_df2 =uq_of_column(out_df2,'new_mr_sites',sep_str='~')
    new_out_df2.to_csv(out_file2, sep='\t', index=False)
    print('Export unique genes in DMRs: ')

    # 8. to find tss overlapping with DMRs
    out_file3 = out_file2.replace('_5dist_', '_')
    cmd = 'grep TSS: ' + out_file2 + '  > ' + out_file3
    os.system(cmd)
    print('Export TSS overlapping with DMRs: ')
    # print('\t',out_file3)

    # 9. to add new column name
    out_df = pd.read_csv(out_file3, sep='\t', header=None)
    out_columns = ['gene_name', 'gene_type', 'new_mr_sites']
    out_df.columns = out_columns
    out_df.to_csv(out_file3, sep='\t', index=False)

    # 10. to find 5dist regions overlapping with enhancers
    out_file3 = in_file.replace('_tss_', '_')
    cmd = 'grep 5dist: ' + in_file + '  > ' + out_file3
    print('Export 5distance region overlapping with enhancers')
    os.system(cmd)

    out_file4 = out_file3.replace(in_data_str, in_data_str + '_overlap_enhancer')
    cmd2 = 'bedtools intersect -a ' + out_file3 + \
           ' -b ' +  in_genome_file + ' -wa -wb ' + \
           '  > ' + out_file4
    os.system(cmd2)
    print('Export 5dist overlapping with enhancers: ')
    # print(out_file4)

    # 11. reformato file  out_file4 for output
    in_file2 = out_file4
    in_df2 = pd.read_csv(in_file2, sep='\t', header=None)
    in_df2['mr_site'] = in_df2[3].apply(lambda x: ':'.join(x.split(':')[0:2]))
    in_df2['gene_type'] = in_df2[5].apply(lambda x: x.split('||')[1])
    out_df2 = in_df2.groupby([6, 'gene_type'])['mr_site'].apply('~'.join).reset_index().copy()
    out_df2.columns = out_columns
    out_file5 = in_file2.replace('dmr_', 'uqGeneDmr_')

    #jbw may 
    #remove duplicated mr in the same row
    new_out_df2 =uq_of_column(out_df2,'new_mr_sites',sep_str='~')
    new_out_df2.to_csv(out_file5, sep='\t', index=False)

    # here new_mr_sites was changed manually to block_id for doing chromSegment_test4blocks
    print('Convert file format to unique gene linked DMRs ..')
    print(out_file5, '\n')

    print('Remove temporary files: ')
    print('\n', out_file2, '\n', out_file3, '\n', out_file4, '\n')
    os.system('rm -f ' + out_file2)
    os.system('rm -f ' + out_file3)
    os.system('rm -f ' + out_file4)
    print('\n')


def find_DEG(in_data_str, in_tss_file, in_5dist_file, in_deg_file, out_result_folder, gene_col_name):

    # this script is used to find DEG genes in TSS and 5distnace regions, then export all of them with a bed format such
    # as chr, start, end, mr_ID, Pvalues et al.



    ###---START to Run --- #######
    # read all data such as DMR-in-TSS, DMR-in-5dist, and DEG files
    print('Read DMR file in TSS regions: ')
    print(in_tss_file, '\n')

    tss_df = pd.read_csv(in_tss_file, sep='\t', header=None)

    print('Read DMR file in 5distance regions: ')
    print(in_5dist_file, '\n')
    dist_df = pd.read_csv(in_5dist_file, sep='\t', header=None)

    print('Read DEG file: ')
    print(in_deg_file, '\n')
    deg_df = pd.read_csv(in_deg_file, sep='\t')
    # deg_df=out_df.copy()

    selected_col_genename = gene_col_name

    # combine DMR-in-TSS, DMR-in-5distance and DEG data to a single data frame
    combined_df = pd.concat([tss_df, dist_df]).copy()
    combined_df.reset_index(inplace=True, drop=True)
    combined_df.columns = ['chrom', 'start_pos', 'end_pos', 'gene_type', 'mr_chrom', 'mr_start_pos', 'mr_end_pos',
                           'mr_id', 'pval']

    # find overlappping between DMR-in-TSS or DMR-in-5istance and DEG , here 'D' means predicted DMR
    # or user can use other criteria to select DMRs from full list.
    dmr_combined_df = combined_df[combined_df.mr_id.apply(lambda x: ':D' in x)].copy()
    dmr_combined_df.reset_index(inplace=True, drop=True)

    # combine DMR with DEG
    deg_genes = set(deg_df[selected_col_genename].to_list())
    dmr_combined_df['gene_name'] = dmr_combined_df.gene_type.apply(lambda x: x.split('||')[2].split(':')[0])

    deg_dmr_combined_df = dmr_combined_df[dmr_combined_df.gene_name.apply(lambda x: x in deg_genes)].copy()
    deg_dmr_combined_df.reset_index(inplace=True, drop=True)

    # export DMR regions are overlapping with either TSS or 5dist of DEGs
    out_df = deg_dmr_combined_df[['mr_chrom', 'mr_start_pos', 'mr_end_pos', 'mr_id', 'pval', 'gene_type', 'gene_name']]
    out_file = os.path.join(out_result_folder, 'dmr_regions_in_deg_tss_5dist' + in_data_str + '.bed')
    out_df.to_csv(out_file, sep='\t', header=None, index=False)
    print('Export combined DMR-in-TSS or DMR-in-5distnace with DEG to file:')
    print('\t ', out_file)


def run(args):
    ###-- INPUT parameters which may need be changed in differen runs-- #####
    # input file folder for MRs overlapping to TSS and 5Distance regions that was done in dmr_analysis
    in_folder = args.in_folder
    in_data_str = args.in_string

    # file name of TSS regions that overlapping with DMRs that located in in_folder
    in_tss_file = args.in_tss_file_mr
    if not os.path.exists(in_tss_file):
        print(in_tss_file , " does not exists.")
        exit(1)

    # file name of 5'distance regions that overlapping with DMRs that located in in_folder
    # dist_file = '3_chroms_all_mr_data_range_dmrRanking_noGenes_5dist_Up1000000_Up5000removedShort_overlap1e-09.bed'

    in_5dist_file = args.in_dist_file
    if not os.path.exists(in_5dist_file):
        print(in_5dist_file , " does not exists.")
        exit(1)
    # DEG: differetial expression gene file in tab delimited format exported by bpb3
    # in_deg_file = '../../data/hap1_cell/in_data/final_demo_data/hap1_cell/in_data/DEG/HAP1_P1_vs_HAP1_KO1_differentially_expressed_genes_min1.1Fd_min1RPKM.txt'
    in_deg_file = args.in_deg_file
    if not os.path.exists(in_deg_file):
        print(in_deg_file , " does not exists.")
        exit(1)
    #column name of gene names
    gene_col_name = args.gene_col_name

    ###---OUTPUT parameters -- ######
    # output file path ?
    out_result_folder = args.out_folder
    if not os.path.exists(out_result_folder):
        os.makedirs(out_result_folder, exist_ok=True)

    find_DEG(in_data_str, in_tss_file, in_5dist_file, in_deg_file, out_result_folder, gene_col_name)




    #### preprocess part ####



    # file name of all DMRs in the result folder of dmr_analysis
    full_mr_file = args.full_mr_file
    print('Read ranked MRs from results of dmr_analysis :')
    print(full_mr_file, '\n')
    if not os.path.exists(full_mr_file):
        print(full_mr_file , " does not exists.")
        exit(1)

    # 3. genome file path which contains refFlat file and enhancer files
    in_genome_file = args.in_genome_file
    if not os.path.exists(in_genome_file):
        print(in_genome_file , " does not exists.")
        exit(1)

    # read TSS regions file from results of dmr_analysis
    tss_file = args.tss_file
    print('Read TSS region file: ')
    print(tss_file, '\n')
    if not os.path.exists(tss_file):
        print(tss_file , " does not exists.")
        exit(1)

    preprocess(full_mr_file, out_result_folder, in_genome_file, in_data_str, tss_file )


if __name__== '__main__':
 args=my_parser(argparse.ArgumentParser('python preprocess.py')).parse_args()
 run(args)








