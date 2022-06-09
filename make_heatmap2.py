import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import os
import re

# input each excel file and combine
root_dir = "/Users/sarahkeegan/Dropbox/mac_files/fenyolab/data_and_results/Rona_FRAP" + \
           "/final_results/final data"

#root_dir = "/Users/sarahkeegan/Dropbox/mac_files/fenyolab/data_and_results/Rona_FRAP" + \
#           "/final_results/normalized data export/Normalized Transform"

nrows=7
fig, axs = plt.subplots(nrows=nrows, ncols=5, figsize=(20,15), ) #sharex=True)
row_i=0
row_i_=0
col_i=0

metric="euclidean" #"correlation" #"euclidean" #"correlation"
method = 'ward' #'complete' #'average' #'ward'
core_proteins=['MSH6', 'MSH3', 'MSH2','MLH1']

# read in all files - combine
files_list=glob.glob(f"{root_dir}/normalized data export for heatmap/txt/Transform of*.txt")
print(len(files_list))
df_all=pd.DataFrame()
df_all_MSH=pd.DataFrame()

for file_ in files_list:
    protein_name=os.path.split(file_)[1]
    mo=re.match(r"(.+)_(.+)\.txt", protein_name)
    protein_name=mo.group(2)
    print(protein_name)

    df=pd.read_csv(file_, sep='\t')
    df = df.set_index('Time (s)', drop=False)

    # Plotting the normalized data:
    sns.lineplot(data=df, x='Time (s)', y='siCCND1/D2/D3', ax=axs[row_i][col_i], label='siAll', color='orange')
    sns.lineplot(data=df, x='Time (s)', y='siNT', ax=axs[row_i][col_i], label='siCtrl', color='blue')

    if (protein_name in ['MSH6', 'MSH3', 'MSH2']):
        sns.lineplot(data=df, x='Time (s)', y='siKIPs', ax=axs[row_i][col_i], label='siAll', color='red')

    axs[row_i][col_i].set_ylabel('Intensity (norm)')
    axs[row_i][col_i].set_title(protein_name)
    if (not (row_i == 0 and col_i == 0)):
        axs[row_i][col_i].get_legend().remove()

    if ((row_i + 1) % nrows == 0):
        row_i = 0
        col_i += 1
    else:
        row_i += 1

    df[protein_name] = df['siCCND1/D2/D3'] - df['siNT']
    df_all = pd.concat([df_all, df[[protein_name,]]], axis=1, ignore_index=False)

    if (protein_name in core_proteins):
        if(protein_name != 'MLH1'):
            #df[f"{protein_name}-siCCND1/D2/D3"] = df['siCCND1/D2/D3'] - df['siNT']
            #df_all_MSH = pd.concat([df_all_MSH, df[[f"{protein_name}-siCCND1/D2/D3", ]]], axis=1, ignore_index=False)
            df[f"{protein_name}-siKIPs"] = df['siKIPs'] - df['siNT']
            df_all_MSH = pd.concat([df_all_MSH, df[[f"{protein_name}-siKIPs", ]]], axis=1, ignore_index=False)

        # Read in the KO data
        files_list = glob.glob(f"{root_dir}/s phase/normalized/Transform of*_{protein_name}_*.csv")
        df = pd.read_csv(files_list[0], sep=',')
        df.rename(columns={'Hours':'Time (s)'}, inplace=True)
        df = df.set_index('Time (s)', drop=False)

        df[f"{protein_name}-KO-1"] = df['C11'] - df['WT']
        df[f"{protein_name}-KO-2"] = df['G11'] - df['WT']
        df_all_MSH = pd.concat([df_all_MSH, df[[f"{protein_name}-KO-1", f"{protein_name}-KO-2"]]],
                               axis=1, ignore_index=False)

        sns.lineplot(data=df, x='Time (s)', y='C11', ax=axs[row_i_][4], label='KO-1', color='green')
        sns.lineplot(data=df, x='Time (s)', y='G11', ax=axs[row_i_][4], label='KO-2', color='darkgreen')
        sns.lineplot(data=df, x='Time (s)', y='WT', ax=axs[row_i_][4], label='WT', color='black')

        axs[row_i_][4].set_ylabel('Intensity (norm)')
        axs[row_i_][4].set_title(f"{protein_name}-KO")
        axs[row_i_][4].set_xlim(0,600)

        if (row_i_ != 0):
            axs[row_i_][4].get_legend().remove()
        row_i_+=1

fig.tight_layout()
fig.savefig(f"{root_dir}/all_norm.pdf")
plt.clf()
plt.close(fig)



df_all['Time (s)']=df_all.index
df_all.index=range(len(df_all))

df_all_MSH['Time (s)']=df_all_MSH.index
df_all_MSH.index=range(len(df_all_MSH))
df_all_MSH=df_all_MSH[
    ['Time (s)',
    'MSH2-siKIPs',
    'MSH3-siKIPs',
    'MSH6-siKIPs',
    'MSH2-KO-1',
    'MSH2-KO-2',
    'MSH3-KO-1',
    'MSH3-KO-2',
    'MSH6-KO-1',
    'MSH6-KO-2',
    'MLH1-KO-1',
    'MLH1-KO-2'
     ]
]


# Make Heatmap
cutoff_in_sec=600 #1200
gap=3 #6
rolling=True
data_type=['siCCND123','MSH']
for i,hm_data in enumerate([df_all,df_all_MSH]):

    data_for_hm=hm_data
    data_for_hm.to_csv(f"{root_dir}/all_data_heatmap-{data_type[i]}.txt", sep='\t')

    data_for_hm=data_for_hm.iloc[1:,:]
    data_for_hm=data_for_hm[data_for_hm['Time (s)']<=cutoff_in_sec]
    time_in_min=np.round(data_for_hm['Time (s)']/60,2)
    data_for_hm.drop(columns=['Time (s)'],inplace=True)

    if(rolling):
        df = data_for_hm.rolling(gap).mean()
        select_idx=range(gap,len(df)+gap,gap)
        df_ = df.loc[select_idx]
        df_.index=time_in_min[select_idx]
    else:
        df_=data_for_hm
        df_.index=time_in_min

    if(i == 0):
        ret=sns.clustermap(df_.transpose(),
                           center=0,
                           col_cluster=False,
                           cmap=sns.diverging_palette(240, 10, n=81),
                           vmin=-0.4, vmax=0.4,
                           colors_ratio=0.01,
                           method=method,
                           metric=metric)

        hm_axes=ret.ax_heatmap

    else:
        grid_kws = {"width_ratios": (.05, 0.9), "wspace": .1,
                    "height_ratios": (.3, 0.7), "hspace": .1}
        f, ax_arr = plt.subplots(2,2, gridspec_kw=grid_kws)

        ret = sns.heatmap(df_.transpose(),
                          center=0,
                          cmap=sns.diverging_palette(240, 10, n=81),
                          #vmin=-0.3, vmax=0.3,
                          square=True,
                          ax=ax_arr[1][1],
                          cbar_ax=ax_arr[0][0],
                          cbar_kws={"orientation": "vertical"}
                          )
        ax_arr[0][1].axis('off')
        ax_arr[1][0].axis('off')
        hm_axes=ret.axes
        hm_axes.tick_params(axis='y', labelsize=6)

    xticks = hm_axes.get_xticks()
    hm_axes.set_xticks(xticks + 1,
                              np.arange(0.5, 10.5, 0.5))  # range(1,21,1)) # <--this one for full time frame!
    hm_axes.set_xlabel('Time (m)')


    if(i==0):
        plt.tight_layout()
    plt.savefig(f"{root_dir}/all_heatmap-v2-{data_type[i]}.pdf")
    plt.clf()
    plt.close()

    df_.index=df_.index*60
    df_.to_csv(f"{root_dir}/all_data_heatmap-ave-{data_type[i]}.txt", sep='\t')
