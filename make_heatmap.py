import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

root_dir = "/Users/sarahkeegan/Dropbox/mac_files/fenyolab/data_and_results/Rona_FRAP/final_results/Organized BER recruitment data only"

data_for_hm=pd.read_csv(f"{root_dir}/all_for_heatmap_NORM.csv")
data_for_hm['Time (m)']=np.round(data_for_hm['Time (s)']/60,1)
data_for_hm.drop(columns=['Time (s)'],inplace=True)
data_for_hm.set_index('Time (m)',drop=True,inplace=True)

df = data_for_hm.rolling(6).mean()
df_ = df.iloc[::6, :]
df_ = df_.iloc[1:,:]

sns.clustermap(df_.transpose(), center=0, col_cluster=False, method='ward', vmin=-0.3, vmax=0.3, cmap=sns.diverging_palette(240, 10, n=36))
plt.tight_layout()
plt.savefig(f"{root_dir}/all_heatmap.png")
