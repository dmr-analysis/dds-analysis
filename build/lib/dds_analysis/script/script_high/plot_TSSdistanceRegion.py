# coding=utf-8
#this script is adapted from HMST-Seq-Analyzer (https://hmst-seq.github.io/hmst/)
#that was developed by Sindre Grønmyr and Junbai Wang
#
import argparse
import pandas as pd
import numpy as np
import scipy.interpolate as sc_i
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from scipy.ndimage import gaussian_filter1d
import datetime

def translate(value, orig_min, orig_max, new_min, new_max):
    """Mapping value in span [orig_min, orig_max] to span [new_min, new_max]
    """
    # Figure out how 'wide' each range is
    orig_span = int(orig_max) - int(orig_min)
    new_span = new_max - new_min

    # Convert the original range into a 0-1 range (float)
    value_scaled = float(value - orig_min) / float(orig_span)

    # Convert the 0-1 range into a value in the new range.
    return new_min + (value_scaled * new_span)

def get_methylation_type(regions):
    """Finds what is the methylation type for each of the files, returns it. If
    it's not the same for all files, return None"""
    all_5mC = []
    all_5hmC = []
    for region in regions:
        for f in regions[region]:
            all_5mC.append('5mC' in f)
            all_5hmC.append('5hmC' in f)
    if all(all_5mC):
        methylation_type = '5mC'
    elif all(all_5hmC):
        methylation_type = '5hmC'
    else:
        methylation_type = None
    return methylation_type

def same_length_regions(regions):
    """Checks if all TSS, geneBody and TES is provided for all samples to be plotted
    """
    length = []
    for region in regions:
        length.append(len(regions[region]))
    #added by jbw
    if sum(length)>0:
      out_term= all(x == length[0] for x in length)
    else:
      #all are empty 
      out_term=False
    return out_term 

def str_list_toArray(str):
    """Helper function that turns a string of an array into an array of numpy floats
    """
    str = str.replace('[', '').replace(']', '').replace(',', ' ')
    str_split = str.split()
    return np.array(str_split).astype(np.float)

def interpolater(df, start, end, step_size):
    """adds two columns to df: interpolated and x_bin. x_bins is positions in (start,end)
    interpolated is the methylation values at each x_bin
    df -- DataFrame
    start -- start position
    end -- end position
    step_size -- step size of arange"""
    df['interpolated'] = None
    df['x_bin'] = None
    for i, row in df.iterrows():
        tmp_pos = []

        for dpoint in row.meth_points_all:
            tmp_pos.append(np.floor(translate(dpoint, row.region_start, row.region_end, start, end)))
        #added wang
        x_bin = np.arange(int(min(tmp_pos)), int(max(tmp_pos)), step_size)
        if len(x_bin) < 1:
            x_bin = [int(tmp_pos[0])]
        df.at[i, 'x_bin'] = x_bin
        #added wang
        if len(row.powers)>2:
           df.at[i, 'interpolated'] = sc_i.interp1d(tmp_pos, row.powers, kind='nearest')(x_bin)
           #print x_bin, sc_i.interp1d(tmp_pos, row.powers, kind='nearest')(x_bin), row.powers, len(row.powers)
        else:
           #print row.powers
           df.at[i, 'interpolated'] = row.powers
           if len(row.powers) <=1:
              df.at[i,'x_bin']=[int(tmp_pos[0])]
              #print [int(tmp_pos[0])]
           elif len(row.powers)==2 :
              df.at[i,'x_bin']=np.asarray([int(tmp_pos[0]) ,int(tmp_pos[1])])
              #print np.asarray([int(tmp_pos[0]) ,int(tmp_pos[1])])
    return df

def find_average(df, s, e):
    """Returns dataframe with indexes between s and e, and the average at each site,
    where missing values is imputed with neighbors.
    s -- start
    e -- end"""
    #added wang
    end = pd.DataFrame(index=np.arange(s, e+1, 1), columns=['avg', 'count'])

    end['avg'] = np.nan
    end['count'] = np.nan
    for i, row in df.iterrows():
        for j, x in enumerate(row.x_bin):
            #if i==343 and j==1:
            #  print i,j , x
            #  print pd.isna(end.at[x,'avg'])

            if pd.isna(end.at[x, 'avg']):
                end.at[x, 'avg'] = row.interpolated[j]
                end.at[x, 'count'] = 1
            else:
                end.at[x, 'avg'] = end.at[x, 'avg'] + row.interpolated[j]
                end.at[x, 'count'] += 1


    end['orig_avg'] = (end['avg']/end['count'])
    return end['orig_avg']

def plot_combined(TSS_df, gene_df, TES_df, methylation_type, sigma, x, y, gene_len, step_size, window, MRs_txt, folder_out):
    """Plots the combined TSS-gene_body-TES figure."""
    combo_df = pd.DataFrame()
    plt.figure(figsize=(14,8))
    for sample in TSS_df.columns:
        frames = [TSS_df[sample], gene_df[sample], TES_df[sample]]
        result = pd.concat(frames)

        rolling = result.rolling(window=window, min_periods=1, center=True).mean()
        smooth = gaussian_filter1d(rolling, sigma=sigma)

        combo_df[sample] = smooth

        len_x = len(smooth)
        plt.plot(result.index, smooth, label=sample)

    combo_df.to_csv(folder_out+'/plotData/plotData_TSSgeneTES_smoothed_'+methylation_type+'_'+MRs_txt+'.csv')


    TSS_line = 0
    TES_line = y*2+gene_len
    plt.axvline(x=TSS_line, color='0', linewidth=1)
    plt.axvline(x=TES_line, color='0', linewidth=1)
    #plt.ylim(0, 3)
    plt.xlim(result.index[0], result.index[-1])
    plt.title(methylation_type+', '+MRs_txt)
    TSS_TES_ticks = [TSS_line, TES_line]
    Y_ticks = [y, y+gene_len]
    plt.xticks(list(np.linspace(-x, result.index[-1], 3))+TSS_TES_ticks+Y_ticks, [str(-x)+'bp', '50%', str(x)+'bp', 'TSS', 'TES', str(y)+'bp', str(-y)+'bp'])
    plt.ylabel('Average Methylation Level')
    plt.legend()
    plt.savefig(folder_out+'/'+methylation_type+'_TSSgeneTES_X'+str(x)+'_Y'+str(y)+'_G'+\
                str(gene_len)+'_binsize'+str(step_size)+'_sigma'+str(sigma)+'_'+MRs_txt+'_'+\
                str(datetime.datetime.now())[:10]+'.jpg',dpi=400)

#added jbw
def plot_combined2regions(TSS_df, gene_df, methylation_type, sigma, x, y, gene_len, step_size, window, MRs_txt, folder_out):
    """Plots the combined two regions such as promoter + 5Distance/Enhancer figure."""
    combo_df = pd.DataFrame()
    plt.figure(figsize=(14,8))
    for sample in TSS_df.columns:
        frames = [TSS_df[sample],gene_df[sample]]
        result = pd.concat(frames)

        rolling = result.rolling(window=window, min_periods=1, center=True).mean()
        smooth = gaussian_filter1d(rolling, sigma=sigma)

        combo_df[sample] = smooth

        len_x = len(smooth)
        plt.plot(result.index, smooth, label=sample)

    combo_df.to_csv(folder_out+'/plotData/plotData_TSS5dist_smoothed_'+methylation_type+'_'+MRs_txt+'.csv')

    #draw an vertical line to separete TSS and enhancer regions
    #jbw
    #TSS_line = x
    #print(x,y)
    TSS_line =y
    plt.axvline(x=TSS_line, color='0', linewidth=1)
    #set xlim location based on result position
    #print(result)
    plt.xlim(result.index[0], result.index[-1])
    plt.title(methylation_type+', '+MRs_txt)
    
    TSS_ticks = [TSS_line]
    #set promoter tick position
    Y_ticks = [0]
    #print(list(np.linspace(y, result.index[-1], 3))+Y_ticks)
    #print(['', 'Enhancer', '', 'Promoter'])
    plt.xticks(list(np.linspace(y, result.index[-1], 3))+Y_ticks, ['', 'Enhancer', '', 'Promoter'])
    plt.ylabel('Average Methylation Level')
    plt.legend()
    out_fig_file=folder_out+'/'+methylation_type+'_TSS5dist_X'+str(x)+'_Y'+str(y)+'_G'+\
                str(gene_len)+'_binsize'+str(step_size)+'_sigma'+str(sigma)+'_'+MRs_txt+'_'+\
                str(datetime.datetime.now())[:10]+'.jpg'

    plt.savefig(out_fig_file,dpi=400)
    print(out_fig_file)


def plot_alone(df, region, methylation_type, sigma, x, y, gene_len, step_size, window, MRs_txt, folder_out):
    """Plots TSS, gene body or TES samples by itself. If there for example exist no
    MRs in TSS, but there are MRs in TES and gene body, then the two will be plotted
    in their own figures."""
    if len(df) < 1:
        return

    if region in ('TSS', 'TES'):
        plt.figure(figsize=(8,8))
    else:
        plt.figure(figsize=(10,8))

    combo_df = pd.DataFrame()
    for sample in df.columns:
        rolling = df[sample].rolling(window=window, min_periods=1, center=True).mean()
        smooth = gaussian_filter1d(rolling, sigma=sigma)

        combo_df[sample] = smooth

        len_x = len(smooth)
        plt.plot(df.index, smooth, label=sample)

    combo_df.to_csv(folder_out+'/plotData/plotData_'+region+'_smoothed_'+methylation_type+'_'+MRs_txt+'.csv')

    if region == 'TSS':
        line = 0
        left_str = str(-x)+'bp'
        right_str = str(y)+'bp'
        middle_str = region
        plt.axvline(x=line, color='0', linewidth=1)
        mid_point = line
    elif region == 'TES':
        line = y*2+gene_len
        left_str = str(-y)+'bp'
        right_str = str(x)+'bp'
        middle_str = region
        plt.axvline(x=line, color='0', linewidth=1)
        mid_point = line
    else:
        left_str = region +  ' Start'
        right_str = region + ' End'
        middle_str = '50%'
        mid_point = (df.index[0]+df.index[-1])/2

    plt.xticks([df.index[0], mid_point, df.index[-1]], [left_str, middle_str, right_str])
    plt.xlim(df.index[0], df.index[-1])
    plt.title(region+', '+methylation_type+', '+MRs_txt)
    plt.ylabel('Average Methylation Level')
    plt.legend()
    out_fig_file=folder_out+'/'+methylation_type+'_'+region+'_X'+str(x)+'_Y'+str(y)+'_G'+\
                str(gene_len)+'_binsize'+str(step_size)+'_sigma'+str(sigma)+'_'+MRs_txt+'_'+\
                str(datetime.datetime.now())[:10]+'.jpg'
    plt.savefig(out_fig_file,dpi=400)
    print(out_fig_file)

def main(tss, tes, geneBody, x, y, gene_len, step_size, window, sigma, MRs_txt, folder_out,plot2regions=False):
    """Plots distribution of 5mc/5hmC levels in TSS, TES, and gene body. If want to plot only 2 joined regions then set plot2regions True (added jbw)"""
    regions = {}
    if tss is not None and len(tss)>0:
        regions['TSS'] = []
        for i in tss:
            regions['TSS'].append(i)
    if tes is not None and len(tes)>0:
        regions['TES'] = []
        for i in tes:
            regions['TES'].append(i)
    #add jbw
    if (geneBody is not None and len(geneBody)>0) and not plot2regions :
        regions['geneBody'] = []
        for i in geneBody:
            regions['geneBody'].append(i)
    elif (geneBody is not None and len(geneBody)>0) and plot2regions:
        regions['Enhancer']= []
        for i in geneBody:
            regions['Enhancer'].append(i)

    # Checking if all files given are all mC or hmC
    methylation_type = get_methylation_type(regions)
    if methylation_type is None:
        print('***** ERROR in plot_TSSgeneTES.py: All files must be of same methylation type, either 5mC or 5hmC. Exiting. *****')
        exit(1)


    conv = {'meth_points_all': str_list_toArray,
            'powers': str_list_toArray}

    TES_df = pd.DataFrame()
    TSS_df = pd.DataFrame()
    gene_df = pd.DataFrame()

    for region in regions:
        for file in regions[region]:
            print(file)
            sample = (file.split('/')[-1]).split('_')[0]
            tmp_df = pd.read_csv(file, index_col=False, header=0, converters=conv)
            if region == 'TSS':
                tmp_df = interpolater(tmp_df, -x, y, step_size)
                tmp_avg = find_average(tmp_df, -x, y)
                TSS_df[sample] = tmp_avg#[-x:y]
            elif region == 'geneBody' or region=='Enhancer':
                #added jbw
                tmp_df = interpolater(tmp_df, y, y+gene_len, step_size)
                tmp_avg = find_average(tmp_df, y, y+gene_len)
                gene_df[sample] = tmp_avg#[0:gene_len]
            elif region == 'TES':
                tmp_df = interpolater(tmp_df, y+gene_len, y+gene_len+(y+x), step_size)
                tmp_avg = find_average(tmp_df, y+gene_len, y+gene_len+(y+x)+1)
                TES_df[sample] = tmp_avg

    plt.rcParams.update({'font.size': 16})
    #sns.set(font_scale=1.5)
    #added jbw
    if plot2regions:
        #print(regions)
        if not same_length_regions(regions) or len(regions)<2:
           #print('plot_alone')
           plot_alone(TSS_df, 'TSS', methylation_type, sigma, x, y, gene_len, step_size, window, MRs_txt, folder_out)
           #plot_alone(TES_df, 'TES', methylation_type, sigma, x, y, gene_len, step_size, window, MRs_txt, folder_out)
           plot_alone(gene_df, 'Enhancer', methylation_type, sigma, x, y, gene_len, step_size, window, MRs_txt, folder_out)
        else:
           plot_combined2regions(TSS_df, gene_df, methylation_type, sigma, x, y, gene_len, step_size, window, MRs_txt, folder_out)
    elif not same_length_regions(regions):
        plot_alone(TSS_df, 'TSS', methylation_type, sigma, x, y, gene_len, step_size, window, MRs_txt, folder_out)
        plot_alone(TES_df, 'TES', methylation_type, sigma, x, y, gene_len, step_size, window, MRs_txt, folder_out)
        plot_alone(gene_df, 'Gene Body', methylation_type, sigma, x, y, gene_len, step_size, window, MRs_txt, folder_out)
    else:
        plot_combined(TSS_df, gene_df, TES_df, methylation_type, sigma, x, y, gene_len, step_size, window, MRs_txt, folder_out)
