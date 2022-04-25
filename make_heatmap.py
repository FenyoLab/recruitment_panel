import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import seaborn as sns
import numpy as np

root_dir = "/Users/snk218/Dropbox/mac_files/fenyolab/data_and_results/Rona_FRAP/final_results/Organized BER recruitment data only"

data_for_hm=pd.read_csv(f"{root_dir}/all_for_heatmap.csv")
data_for_hm.set_index('Time (s)',drop=True,inplace=True)


for protein in data_for_hm.columns:
    print(protein)
    print(np.percentile(data_for_hm[protein], [25,50,75]))

fig,ax=plt.subplots()
sns.clustermap(data_for_hm.transpose(), center=0, col_cluster=False, vmin=-50, vmax=50)
fig.tight_layout()
fig.savefig(f"{root_dir}/all_heatmap.png")
