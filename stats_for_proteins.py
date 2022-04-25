import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import seaborn as sns
import numpy as np
from scipy.stats import mannwhitneyu

root_dir = "/Users/snk218/Dropbox/mac_files/fenyolab/data_and_results/Rona_FRAP/final_results/Organized BER recruitment data only"
n_timepoints=242

dir_list = next(os.walk(root_dir))[1]
num_proteins = len(dir_list)
make_plot=False
start_t = 50
end_t = 400
pvalue_matrix=[]
fig, axs = plt.subplots(nrows=7, ncols=5, figsize=(20,15),sharex=True) # room for 35 proteins (there are 31)
row_i=0
col_i=0
data_for_hm = pd.DataFrame()
for i,cur_dir in enumerate(dir_list):
    #if(not os.path.exists(f"{root_dir}/{cur_dir}/SplitData.xlsx")):
    #    print(f"Error: did not locate SplitData.xlsx in {cur_dir}. Skipping...")
    #    continue

    df_all=pd.read_excel(f"{root_dir}/{cur_dir}/SplitData.xlsx", sheet_name="SIALLDIFF", nrows=n_timepoints+1)
    df_ctrl=pd.read_excel(f"{root_dir}/{cur_dir}/SplitData.xlsx", sheet_name="SICTRLDIFF", nrows=n_timepoints+1)

    protein_name = re.split("_",cur_dir,)[1]
    print(protein_name)

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

    df_all=df_all.melt(id_vars=['Time (s)'], value_vars=df_all.columns[df_all.columns != "Time (s)"], var_name='cell', value_name='intensity-diff')
    df_all['cell']=df_all['cell'].astype('str') + '_0'
    df_all['condition']='siALL'

    df_ctrl = df_ctrl.melt(id_vars=['Time (s)'], value_vars=df_ctrl.columns[df_ctrl.columns != "Time (s)"], var_name='cell', value_name='intensity-diff')
    df_ctrl['cell'] = df_ctrl['cell'].astype('str') + '_1'
    df_ctrl['condition']='siCtrl'

    df_full=pd.concat([df_all,df_ctrl], ignore_index=True)

    # Plot mean with error shading of intensity-difference over time, for each protein
    if(make_plot):
        sns.lineplot(data=df_full, x='Time (s)', y='intensity-diff', hue='condition', ax=axs[row_i][col_i])
        axs[row_i][col_i].set_title(protein_name)

        axs[row_i][col_i].axvline(50,ls='--',color='black')
        axs[row_i][col_i].axvline(400,ls='--',color='black')

        if (not (row_i == 0 and col_i == 0)):
            axs[row_i][col_i].get_legend().remove()

    # Heatmap
    # mean all cells for each time point
    # normalize between 0 and 1 (1 is the max Intensity of all values, both siAll and siCtrl)
    df_all_for_hm=df_all.groupby('Time (s)').mean()['intensity-diff']
    df_ctrl_for_hm=df_ctrl.groupby('Time (s)').mean()['intensity-diff']

    max_val=np.max([df_ctrl_for_hm.max(),df_all_for_hm.max()])
    df_all_for_hm=df_all_for_hm/max_val
    df_ctrl_for_hm=df_ctrl_for_hm/max_val

    data_for_hm[protein_name]=df_all_for_hm-df_ctrl_for_hm

    # Subtract the intensity-difference at time 0 from each intensity
    # Take the mean intensity of siALL from time=50s to time=400s
    # Take the mean intensity of siCtrl from time=50s to time=400s
    # Use ? to test for significance and get a p-value
    gap=50
    for start_t in range(25,625+gap,gap):
        end_t=start_t+gap

        df_all_filt=df_all[(df_all['Time (s)']>=start_t) & (df_all['Time (s)']<=end_t)]
        df_all_data=df_all_filt.groupby('cell').mean()['intensity-diff']

        df_ctrl_filt = df_ctrl[(df_ctrl['Time (s)'] >= start_t) & (df_ctrl['Time (s)'] <= end_t)]
        df_ctrl_data = df_ctrl_filt.groupby('cell').mean()['intensity-diff']

        orig_len=len(df_all_data)
        df_all_data=df_all_data[df_all_data>0]
        print(f"siAll: Negatives removed: {orig_len-len(df_all_data)} ({len(df_all_data)})")

        orig_len=len(df_ctrl_data)
        df_ctrl_data = df_ctrl_data[df_ctrl_data > 0]
        print(f"siCtrl: Negatives removed: {orig_len-len(df_ctrl_data)} ({len(df_ctrl_data)})")

        # Use Mann Whitney?
        U1, p = mannwhitneyu(df_all_data, df_ctrl_data)
        #print(p)

        pvalue_matrix.append([
            protein_name,
            start_t,
            end_t,
            p,
        ])

    if ((row_i + 1) % 7 == 0):
        row_i = 0
        col_i += 1
    else:
        row_i += 1

pvalue_df = pd.DataFrame(pvalue_matrix, columns=['protein','start_t','end_t','p-value'])
pvalue_df.to_csv(f"{root_dir}/all_pvalues.csv")

data_for_hm.to_csv(f"{root_dir}/all_for_heatmap2.csv")

if(make_plot):
    fig.tight_layout()
    fig.savefig(f"{root_dir}/all_raw.png")