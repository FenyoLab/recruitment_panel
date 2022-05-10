import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import seaborn as sns
import numpy as np
from scipy.stats import mannwhitneyu

root_dir = "/Users/snk218/Dropbox/mac_files/fenyolab/data_and_results/Rona_FRAP/final_results/Organized BER recruitment data only"
n_timepoints=242
subtract_t0=True
p_value_window = 50
p_cutoff=0.05

dir_list = next(os.walk(root_dir))[1]
num_proteins = len(dir_list)
make_plot=True

pvalue_matrix=[]
fig, axs = plt.subplots(nrows=7, ncols=5, figsize=(20,15),sharex=True) # room for 35 proteins (there are 31)
fig2, axs2 = plt.subplots(nrows=7, ncols=5, figsize=(20,15),sharex=True) # room for 35 proteins (there are 31)
row_i=0
col_i=0

data_for_hm_norm = pd.DataFrame()

for i,cur_dir in enumerate(dir_list):
    df_all=pd.read_excel(f"{root_dir}/{cur_dir}/SplitData.xlsx", sheet_name="SIALLDIFF", nrows=n_timepoints+1)
    df_ctrl=pd.read_excel(f"{root_dir}/{cur_dir}/SplitData.xlsx", sheet_name="SICTRLDIFF", nrows=n_timepoints+1)

    protein_name = re.split("_",cur_dir,)[1]
    print(protein_name)

    # PROCESSING INPUT
    # removed "Unnamed", "mean"
    df_all.drop(columns=df_all.filter(regex=r"^Unnamed").columns,inplace=True)
    df_ctrl.drop(columns=df_ctrl.filter(regex=r"^Unnamed").columns,inplace=True)
    if("mean" in df_all.columns):
        df_all.drop(columns=["mean",],inplace=True)
    if("mean" in df_ctrl.columns):
        df_ctrl.drop(columns=["mean",],inplace=True)

    # remove null row - there should be one null row as the last row in the df's
    # checking this - to make sure I am reading in all data each time correctly
    null_idx=df_all.index[df_all.isnull().all(1)]
    if(len(null_idx)!=1 or null_idx[0]!=n_timepoints):
        print(f"Error: SIALLDIFF has an unexpected number of rows in SplitData.xlsx for {cur_dir}")
    else:
        df_all.drop(null_idx,axis=0,inplace=True)
    null_idx=df_ctrl.index[df_ctrl.isnull().all(1)]
    if(len(null_idx)!=1 or null_idx[0]!=n_timepoints):
        print(f"Error: SICTRLDIFF has an unexpected number of rows in SplitData.xlsx for {cur_dir}")
    else:
        df_ctrl.drop(null_idx,axis=0,inplace=True)

    # Subtract initial starting value from all values so first time point intensity difference is 0
    if(subtract_t0):
        df_all[df_all.columns[df_all.columns != "Time (s)"]] = df_all[df_all.columns[df_all.columns != "Time (s)"]] - \
                                                               df_all[df_all.columns[df_all.columns != "Time (s)"]].iloc[0]
        df_ctrl[df_ctrl.columns[df_ctrl.columns != "Time (s)"]] = df_ctrl[df_ctrl.columns[df_ctrl.columns != "Time (s)"]] - \
                                                                  df_ctrl[df_ctrl.columns[df_ctrl.columns != "Time (s)"]].iloc[0]

    # MELT THE DATA FRAMES
    df_all=df_all.melt(id_vars=['Time (s)'], value_vars=df_all.columns[df_all.columns != "Time (s)"], var_name='cell', value_name='intensity-diff')
    df_all['cell']=df_all['cell'].astype('str') + '_0'
    df_all['condition']='siALL'

    df_ctrl = df_ctrl.melt(id_vars=['Time (s)'], value_vars=df_ctrl.columns[df_ctrl.columns != "Time (s)"], var_name='cell', value_name='intensity-diff')
    df_ctrl['cell'] = df_ctrl['cell'].astype('str') + '_1'
    df_ctrl['condition']='siCtrl'

    df_full=pd.concat([df_all,df_ctrl], ignore_index=True)

    # STATISTICS OVER TIME WINDOWS FROM START TO FINISH
    # Take the mean intensity of siALL/siCtrl for time windows of 50s
    # Use Mann Whitney to test for significance and get a p-value for each time window
    times_sig=[]
    times_sig_=[]
    max_time=int(np.max(df_all['Time (s)']))
    for start_t in range(25, 1150 + p_value_window, p_value_window):
        #print(start_t)
        end_t = start_t + p_value_window

        df_all_filt = df_all[(df_all['Time (s)'] >= start_t) & (df_all['Time (s)'] <= end_t)]
        df_all_data = df_all_filt.groupby('cell').mean()['intensity-diff']

        df_ctrl_filt = df_ctrl[(df_ctrl['Time (s)'] >= start_t) & (df_ctrl['Time (s)'] <= end_t)]
        df_ctrl_data = df_ctrl_filt.groupby('cell').mean()['intensity-diff']

        num_neg_siAll = len(df_all_data[df_all_data < 0])
        df_all_data[df_all_data<0]=0

        num_neg_siCtrl = len(df_ctrl_data[df_ctrl_data < 0])
        df_ctrl_data[df_ctrl_data<0]=0

        # Use Mann Whitney
        U1, p = mannwhitneyu(df_all_data, df_ctrl_data)
        # print(p)
        if(p < p_cutoff):
            sig=1
            times_sig_.append(start_t)
        else:
            sig=0
            if(len(times_sig_)>0):
                times_sig.append(times_sig_)
                times_sig_=[]

        pvalue_matrix.append([
            protein_name,
            start_t,
            min(end_t,max_time),
            num_neg_siAll,
            len(df_all_data),
            num_neg_siCtrl,
            len(df_ctrl_data),
            p,
            sig
        ])
    if(times_sig_!=[]):
        times_sig.append(times_sig_)

    p_times=[]
    if(len(times_sig) > 0):
        times_sig_len=[len(x) for x in times_sig]
        max_index = times_sig_len.index(max(times_sig_len))
        times_sig=times_sig[max_index]

        p_times.append(times_sig[0])
        p_times.append(min(times_sig[-1]+p_value_window,max_time))


    # PLOT
    colors_dict={'siALL':'orange','siCtrl':'blue'}
    # Plot mean with error shading of intensity-difference over time, for each protein
    if(make_plot):     # (have a flag here b/c plot takes a long time to make)
        sns.lineplot(data=df_full, x='Time (s)', y='intensity-diff', hue='condition', ax=axs[row_i][col_i], palette=colors_dict)
        axs[row_i][col_i].set_title(protein_name)
        axs[row_i][col_i].set_ylabel('Intensity')

        if(p_times != []):

            axs[row_i][col_i].axvline(p_times[0],ls='--',color='black')
            axs[row_i][col_i].axvline(p_times[1],ls='--',color='black')

        if (not (row_i == 0 and col_i == 0)):
            axs[row_i][col_i].get_legend().remove()

    # MEANS AND NORMALIZATION 0-1
    # ANOTHER PLOT - just showing norm'd means
    # SAVE DATA FOR HEATMAP
    # mean all cells for each time point
    # value at each time point = mean(siAll)-mean(siCtrl)
    df_mean = pd.DataFrame()
    df_mean['Time (s)'] = df_all['Time (s)'].unique()
    df_mean.index=df_mean['Time (s)']
    df_mean['siAll'] = df_all.groupby('Time (s)').mean()['intensity-diff']
    df_mean['siCtrl'] = df_ctrl.groupby('Time (s)').mean()['intensity-diff']

    # divide by max value over siAll and siCtrl so that data will be between 0 and 1
    max_val=np.max([df_mean['siAll'].max(),df_mean['siCtrl'].max()])
    min_val = np.min([df_mean['siAll'].min(), df_mean['siCtrl'].min()])

    #df_mean['siAll-norm'] = (df_mean['siAll']-min_val) / (max_val-min_val)
    #df_mean['siCtrl-norm'] = (df_mean['siCtrl']-min_val) / (max_val-min_val)

    df_mean['siAll-norm'] = df_mean['siAll'] / max_val
    df_mean['siCtrl-norm'] = df_mean['siCtrl'] / max_val

    # Plotting the normalized data:
    sns.lineplot(data=df_mean, x='Time (s)', y='siAll-norm', ax=axs2[row_i][col_i], label='siAll', color='orange')
    sns.lineplot(data=df_mean, x='Time (s)', y='siCtrl-norm', ax=axs2[row_i][col_i], label='siCtrl', color='blue')
    axs2[row_i][col_i].set_ylabel('Intensity (norm)')
    axs2[row_i][col_i].set_title(protein_name)
    if (not (row_i == 0 and col_i == 0)):
        axs2[row_i][col_i].get_legend().remove()

    data_for_hm_norm[protein_name]=df_mean['siAll-norm'] - df_mean['siCtrl-norm']

    if ((row_i + 1) % 7 == 0):
        row_i = 0
        col_i += 1
    else:
        row_i += 1

if(subtract_t0):
    suff='_subtract-t0'
else:
    suff=''

# SAVE STATS, HEATMAP DATA, PLOTS
pvalue_df = pd.DataFrame(pvalue_matrix, columns=['protein',
                                                 'start_t',
                                                 'end_t',
                                                 'siAll_neg',
                                                 'siAll_n',
                                                 'siCtrl_neg',
                                                 'siCtrl_n',
                                                 'MW p-value',
                                                 'sig'])
pvalue_df.to_csv(f"{root_dir}/all_pvalues{suff}-{p_cutoff}-{p_value_window}.csv")

data_for_hm_norm.to_csv(f"{root_dir}/all_for_heatmap{suff}.csv")

fig2.tight_layout()
fig2.savefig(f"{root_dir}/all_norm{suff}.png")

if(make_plot):
    fig.tight_layout()
    fig.savefig(f"{root_dir}/all_raw{suff}-{p_cutoff}-{p_value_window}.png")