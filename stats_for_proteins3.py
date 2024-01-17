import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import seaborn as sns
import numpy as np
from scipy.stats import mannwhitneyu
import getpass as gt

uid=gt.getuser()

root_dir = f"/Users/{uid}/Dropbox (NYU Langone Health)/mac_files/fenyolab/data_and_results/Rona_FRAP"
sub_dirs = ["/final_results/paper_revision_experiments/all_data",]

n_timepoints = 242
p_value_window = 60
p_cutoff = 0.01

dir_list = next(os.walk(root_dir+"/"+sub_dirs[0]))[1]
num_proteins = len(dir_list)
make_plot = True

sheet_names = ['siNT', 'siALL', 'siKIPs']
rescue_sheet_names = ['EV_siNT', 'EV_siUTR', 'WT_siUTR', 'HP_siUTR']

pvalue_matrix = []
if make_plot:
    ncols = 2
    fig, axs = plt.subplots(nrows=2, ncols=ncols, figsize=(20, 15),)
    row_i = 0
    row_i_ = 0
    col_i = 0


def read_and_process_sheet(data_sheets, sheet_names, sub_dir_index, cur_dir, n_tp):

    for sheet_name in sheet_names:
        data_sheets[sheet_name] = pd.read_excel(f"{root_dir}/{sub_dirs[sub_dir_index]}/{cur_dir}/SplitData.xlsx",
                                                sheet_name=sheet_name,
                                                nrows=n_tp)
    # PROCESSING INPUT
    df_full = pd.DataFrame()
    for i, sheet_name in enumerate(sheet_names):

        # removed "Unnamed", "mean"
        data_sheets[sheet_name].drop(columns=data_sheets[sheet_name].filter(regex=r"^Unnamed").columns,
                                     inplace=True)
        if ("mean" in data_sheets[sheet_name].columns):
            data_sheets[sheet_name].drop(columns=["mean", ], inplace=True)

        if ("time" in data_sheets[sheet_name].columns):
            # rename to "Time (s)"
            data_sheets[sheet_name].rename(columns={"time": "Time (s)"}, inplace=True)

        # Subtract initial starting value from all values so first time point intensity difference is 0
        cell_cols = data_sheets[sheet_name].columns[data_sheets[sheet_name].columns != "Time (s)"]
        data_sheets[sheet_name][cell_cols] = data_sheets[sheet_name][cell_cols] - \
                                             data_sheets[sheet_name][cell_cols].iloc[0]

        # MELT THE DATA FRAMES
        data_sheets[sheet_name] = data_sheets[sheet_name].melt(id_vars=['Time (s)'],
                                                               value_vars=cell_cols,
                                                               var_name='cell',
                                                               value_name='intensity-diff')
        data_sheets[sheet_name]['cell'] = data_sheets[sheet_name]['cell'].astype('str') + '_' + str(i)
        data_sheets[sheet_name]['condition'] = sheet_name

        df_full = pd.concat([df_full, data_sheets[sheet_name]], ignore_index=True)

    return df_full, data_sheets


for i, cur_dir in enumerate(dir_list):

    protein_name = cur_dir
    print(protein_name)

    if protein_name.endswith("res"):
        sheet_names_ = rescue_sheet_names
    else:
        sheet_names_ = sheet_names[:]

    data_sheets = {}
    df_full, data_sheets = read_and_process_sheet(data_sheets, sheet_names_, 0, cur_dir, n_timepoints)

    # STATISTICS OVER TIME WINDOWS FROM START TO FINISH ################################################################
    # Take the mean intensity of siALL/siCtrl for sliding windows of 60sec across the first 10min
    # Use Mann Whitney to test for significance and get a p-value for each time window
    times_sig = {'siALL-siNT': [], 'siKIPs-siNT': [], 'EV_siUTR-EV_siNT': [], 'WT_siUTR-EV_siNT': [],
                 'HP_siUTR-EV_siNT': [], 'HP_siUTR-EV_siUTR': [], }
    times_sig_ = {'siALL-siNT': [], 'siKIPs-siNT': [], 'EV_siUTR-EV_siNT': [], 'WT_siUTR-EV_siNT': [],
                  'HP_siUTR-EV_siNT': [], 'HP_siUTR-EV_siUTR': [], }
    p_times = {'siALL-siNT': [], 'siKIPs-siNT': [], 'EV_siUTR-EV_siNT': [], 'WT_siUTR-EV_siNT': [],
               'HP_siUTR-EV_siNT': [], 'HP_siUTR-EV_siUTR': [], }

    # ignore first 20 seconds, go up to 10 minutes
    # slide by every 5 seconds (the time step of the movies)
    start_time = 20  # 20 sec
    max_time = 10*60  # 10 min # int(np.max(df_all['Time (s)']))
    time_increment = 5  # 5 sec
    p_value_window = 60  # look for significance in average intensity over 1 minute intervals

    for start_t in range(start_time, max_time+time_increment, time_increment):
        #print(start_t)
        end_t = start_t + p_value_window
        if end_t > max_time:
            break

        pvalue_data = {}
        for sheet_name in data_sheets.keys():
            cur_df = data_sheets[sheet_name]
            cur_data = cur_df[(cur_df['Time (s)'] >= start_t) &
                              (cur_df['Time (s)'] <= end_t)].groupby('cell')['intensity-diff'].mean()

            num_neg = len(cur_data[cur_data < 0])

            cur_data[cur_data < 0] = 0
            pvalue_data[sheet_name] = cur_data

        # Use Mann Whitney - siCCNDs / siKIPs
        p = {'siALL-siNT': '', 'siKIPs-siNT': '', 'EV_siUTR-EV_siNT': '', 'WT_siUTR-EV_siNT': '',
             'HP_siUTR-EV_siNT': '', 'HP_siUTR-EV_siUTR': '', }
        for key in data_sheets.keys():
            if key == 'siALL' or key == 'siKIPs':
                # compare to siCtrl
                compare_keys = ['siNT', ]
            elif key == 'EV_siUTR' or key == 'WT_siUTR':
                compare_keys = ['EV_siNT', ]
            elif key == 'HP_siUTR':
                compare_keys = ['EV_siNT', 'EV_siUTR', ]
            else:
                continue

            for compare_key in compare_keys:
                key_ = f"{key}-{compare_key}"

                U1, p[f"{key_}"] = mannwhitneyu(pvalue_data[compare_key], pvalue_data[key])
                if p[key_] < p_cutoff:
                    times_sig_[key_].append(start_t)
                else:
                    if len(times_sig_[key_]) > 0:
                        times_sig[key_].append(times_sig_[key_])
                        times_sig_[key_] = []

        data_list = [
            protein_name,
            start_t,
            end_t
        ]

        for key in ['siALL-siNT', 'siKIPs-siNT', 'EV_siUTR-EV_siNT', 'WT_siUTR-EV_siNT',
                    'HP_siUTR-EV_siNT', 'HP_siUTR-EV_siUTR']:
            data_list.append(p[key])

        for key in ['siNT', 'siALL', 'siKIPs', 'EV_siNT', 'EV_siUTR', 'WT_siUTR', 'HP_siUTR']:
            if key in pvalue_data.keys():
                data_list.append(len(pvalue_data[key]))
            else:
                data_list.append('')

        pvalue_matrix.append(data_list)

    for key in times_sig.keys():

        if times_sig_[key]:
            times_sig[key].append(times_sig_[key])

        if len(times_sig[key]) > 1:
            print(f"*Note: multiple significant intervals for protein ({key})!*")

        if len(times_sig[key]) > 0:
            times_sig_len = [len(x) for x in times_sig[key]]
            max_index = times_sig_len.index(max(times_sig_len))
            times_sig[key] = times_sig[key][max_index]

            p_times[key].append(times_sig[key][0])
            p_times[key].append(times_sig[key][-1]+p_value_window)


    # PLOT #############################################################################################################
    colors_dict = {'siNT': 'black', 'siALL': 'red', 'siKIPs': 'blue',
                   'EV_siNT': 'red', 'EV_siUTR': 'blue', 'WT_siUTR': 'black', 'HP_siUTR': 'purple',
                   'siALL-siNT': 'red', 'siKIPs-siNT': 'blue',
                   'EV_siUTR-EV_siNT': 'blue', 'WT_siUTR-EV_siNT': 'black',
                   'HP_siUTR-EV_siNT': 'purple', 'HP_siUTR-EV_siUTR': 'green'}

    # Plot mean with error shading of intensity-difference over time, for each protein
    if make_plot:     # (have a flag here b/c plot takes a long time to make)
        sns.lineplot(data=df_full,
                     x='Time (s)',
                     y='intensity-diff',
                     hue='condition',
                     ax=axs[row_i][col_i],
                     palette=colors_dict)

        axs[row_i][col_i].set_title(protein_name)
        axs[row_i][col_i].set_ylabel('Intensity')

        for key in p_times.keys():
            if p_times[key]:
                axs[row_i][col_i].axvline(p_times[key][0], ls='--', color=colors_dict[key])
                axs[row_i][col_i].axvline(p_times[key][1], ls='--', color=colors_dict[key])

        if (row_i + 1) % 2 == 0:
            row_i = 0
            col_i += 1
        else:
            row_i += 1

# SAVE STATS, HEATMAP DATA, PLOTS
pvalue_df = pd.DataFrame(pvalue_matrix, columns=['protein',
                                                 'start_t',
                                                 'end_t',
                                                 'siAll-siNT_pvalue',
                                                 'siKIPs-siNT_pvalue',
                                                 'EV_siUTR-EV_siNT_pvalue',
                                                 'WT_siUTR-EV_siNT_pvalue',
                                                 'HP_siUTR-EV_siNT_pvalue',
                                                 'HP_siUTR-EV_siUTR_pvalue',
                                                 'siNT_n',
                                                 'siAll_n',
                                                 'siKIPs_n',
                                                 'EV_siNT_n',
                                                 'EV_siUTR_n',
                                                 'WT_siUTR_n',
                                                 'HP_siUTR_n'])

pvalue_df.to_csv(f"{root_dir}/all_pvalues-REV-{p_value_window}.txt", sep='\t')

if make_plot:
    fig.tight_layout()
    fig.savefig(f"{root_dir}/all_raw-REV-{p_cutoff}-{p_value_window}.pdf")


