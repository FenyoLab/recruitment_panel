import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

root_dir = "/Users/snk218/Dropbox/mac_files/fenyolab/data_and_results/Rona_FRAP/final_results/Organized BER recruitment data only"

subtract_t0=True
if(subtract_t0):
    suff='_subtract-t0'
else:
    suff=''
cutoff_in_sec=600 #1200
gap=3 #6
rolling=True

data_for_hm=pd.read_csv(f"{root_dir}/all_for_heatmap{suff}.csv")
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


ret=sns.clustermap(df_.transpose(), center=0, col_cluster=False, method='ward', vmin=-0.3, vmax=0.3, cmap=sns.diverging_palette(240, 10, n=81))

xticks=ret.ax_heatmap.get_xticks()
ret.ax_heatmap.set_xticks(xticks+1,np.arange(0.5,10.5,0.5)) #range(1,21,1)) # <--this one for full time frame!
ret.ax_heatmap.set_xlabel('Time (m)')

plt.tight_layout()
plt.savefig(f"{root_dir}/all_heatmap{suff}.png")
