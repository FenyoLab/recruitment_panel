import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import os
import re

# input each excel file and combine
root_dir = "/Users/snk218/Dropbox/mac_files/fenyolab/data_and_results/Rona_FRAP/final_results/normalized data export/Normalized Transform"

# read in all files - combine
files_list=glob.glob(f"{root_dir}/Transform of*.txt")
df_all=pd.DataFrame()
for file_ in files_list:
    protein_name=os.path.split(file_)[1]
    mo=re.match(r"(.+)_(.+)\.txt", protein_name)
    protein_name=mo.group(2)
    print(protein_name)

    df=pd.read_csv(file_, sep='\t')

    df[protein_name]=df['siCCND1/D2/D3']-df['siNT']
    df = df[['Time (s)', protein_name]]
    df=df.set_index('Time (s)', drop=True)

    df_all=pd.concat([df_all,df], axis=1, ignore_index=False)
df_all['Time (s)']=df_all.index
df_all.index=range(len(df_all))

# Make Heatmap
cutoff_in_sec=600 #1200
gap=3 #6
rolling=True

data_for_hm=df_all
data_for_hm.to_csv(f"{root_dir}/all_data_heatmap.txt", sep='\t')

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


ret=sns.clustermap(df_.transpose(), center=0, col_cluster=False, cmap=sns.diverging_palette(240, 10, n=81),
                   vmin=-0.4, vmax=0.4,
                   colors_ratio=0.01,
                   method='ward')

xticks=ret.ax_heatmap.get_xticks()
ret.ax_heatmap.set_xticks(xticks+1,np.arange(0.5,10.5,0.5)) #range(1,21,1)) # <--this one for full time frame!
ret.ax_heatmap.set_xlabel('Time (m)')

plt.tight_layout()
plt.savefig(f"{root_dir}/all_heatmap-v2.png")

plt.savefig(f"{root_dir}/all_heatmap-v2.pdf")

df_.index=df_.index*60
df_.to_csv(f"{root_dir}/all_data_heatmap-ave.txt", sep='\t')
