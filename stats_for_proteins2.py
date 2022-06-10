import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import seaborn as sns
import numpy as np
from scipy.stats import mannwhitneyu
import getpass as gt

uid=gt.getuser()

root_dir = f"/Users/{uid}/Dropbox/mac_files/fenyolab/data_and_results/Rona_FRAP"
sub_dirs=["/final_results/Organized Recruitments data only/G1 recruitment siCCNDs",
          "/final_results/Organized Recruitments data only/S phase"]

n_timepoints=242
p_value_window=60
p_cutoff=0.01

dir_list = next(os.walk(root_dir+"/"+sub_dirs[0]))[1]
num_proteins = len(dir_list)
make_plot=True

sheet_names=['SICTRLDIFF','SIALLDIFF','SIKIPsDIFF']
sheet_names2=['AMBRA_WT_DIFF','AMBRA_KO_C11_DIFF','AMBRA_KO_G11_DIFF']
KO_Sphase_proteins=['MSH6', 'MSH3', 'MSH2', 'MLH1']
siKIPs_proteins=['MSH6', 'MSH3', 'MSH2']

pvalue_matrix=[]
if(make_plot):
    ncols=5
    fig, axs = plt.subplots(nrows=7, ncols=ncols, figsize=(20,15),) #sharex=True) # room for 35 proteins (there are 31)
    row_i=0
    row_i_=0
    col_i=0

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

for i,cur_dir in enumerate(dir_list):

    protein_name = cur_dir
    print(protein_name)

    # G1 phase, siCCNDs and siKIPs
    if(not protein_name in siKIPs_proteins):
        sheet_names_=sheet_names[:-1]
    else:
        sheet_names_=sheet_names[:]

    data_sheets = {}
    df_full, data_sheets=read_and_process_sheet(data_sheets, sheet_names_, 0, cur_dir, n_timepoints)

    # S-phase KO cell lines
    if(protein_name in KO_Sphase_proteins):
        df_full_KO_Sphase, data_sheets = read_and_process_sheet(data_sheets, sheet_names2, 1, cur_dir, int(n_timepoints/2+1))

    # STATISTICS OVER TIME WINDOWS FROM START TO FINISH ################################################################
    # Take the mean intensity of siALL/siCtrl for sliding windows of 60sec across the first 10min
    # Use Mann Whitney to test for significance and get a p-value for each time window
    times_sig = {'SIALLDIFF':[], 'SIKIPsDIFF':[], 'AMBRA_KO_C11_DIFF':[], 'AMBRA_KO_G11_DIFF':[]}
    times_sig_ = {'SIALLDIFF':[], 'SIKIPsDIFF':[], 'AMBRA_KO_C11_DIFF':[], 'AMBRA_KO_G11_DIFF':[]}
    p_times = {'SIALLDIFF':[], 'SIKIPsDIFF':[], 'AMBRA_KO_C11_DIFF':[], 'AMBRA_KO_G11_DIFF':[]}

    # ignore first 20 seconds, go up to 10 minutes
    # slide by every 5 seconds (the time step of the movies)
    start_time=20 # 20 sec
    max_time=10*60 # 10 min # int(np.max(df_all['Time (s)']))
    time_increment=5 # 5 sec
    p_value_window=60 # look for significance in average intensity over 1 minute intervals

    for start_t in range(start_time, max_time+time_increment, time_increment):
        #print(start_t)
        end_t = start_t + p_value_window
        if(end_t > max_time):
            break

        pvalue_data={}
        for i, sheet_name in enumerate(data_sheets.keys()):
            cur_df = data_sheets[sheet_name]
            cur_data = cur_df[(cur_df['Time (s)'] >= start_t) &
                              (cur_df['Time (s)'] <= end_t)].groupby('cell').mean()['intensity-diff']

            num_neg = len(cur_data[cur_data < 0])

            cur_data[cur_data<0]=0
            pvalue_data[sheet_name] = cur_data

        # Use Mann Whitney - siCCNDs / siKIPs
        p = {'SIALLDIFF': '', 'SIKIPsDIFF': '', 'AMBRA_KO_C11_DIFF':'', 'AMBRA_KO_G11_DIFF':''}
        for key in data_sheets.keys():
            if(key == 'SIALLDIFF' or key == 'SIKIPsDIFF'): # compare to siCtrl
                compare_key='SICTRLDIFF'
            elif(key == 'AMBRA_KO_C11_DIFF' or key == 'AMBRA_KO_G11_DIFF'):
                compare_key='AMBRA_WT_DIFF'
            else:
                continue

            U1, p[key] = mannwhitneyu(pvalue_data[compare_key], pvalue_data[key])
            if(p[key] < p_cutoff):
                times_sig_[key].append(start_t)
            else:
                if(len(times_sig_[key])>0):
                    times_sig[key].append(times_sig_[key])
                    times_sig_[key]=[]

        data_list=[
            protein_name,
            start_t,
            end_t,
            p['SIALLDIFF'],
            p['SIKIPsDIFF'],
            p['AMBRA_KO_C11_DIFF'],
            p['AMBRA_KO_G11_DIFF'],
            len(pvalue_data['SICTRLDIFF']),
            len(pvalue_data['SIALLDIFF'])]

        for key in ['SIKIPsDIFF','AMBRA_WT_DIFF','AMBRA_KO_C11_DIFF','AMBRA_KO_G11_DIFF']:
            if (key in pvalue_data.keys()):
                data_list.append(len(pvalue_data[key]))
            else:
                data_list.append('')

        pvalue_matrix.append(data_list)

    for key in times_sig.keys():

        if(times_sig_[key]!=[]):
            times_sig[key].append(times_sig_[key])

        if (len(times_sig[key]) > 1):
            print(f"*Note: multiple significant intervales for protein ({key})!*")

        if(len(times_sig[key]) > 0):
            times_sig_len=[len(x) for x in times_sig[key]]
            max_index = times_sig_len.index(max(times_sig_len))
            times_sig[key]=times_sig[key][max_index]

            p_times[key].append(times_sig[key][0])
            p_times[key].append(times_sig[key][-1]+p_value_window)


    # PLOT #############################################################################################################
    colors_dict={'SICTRLDIFF':'blue', 'SIALLDIFF':'orange','SIKIPsDIFF':'red',
                 'AMBRA_WT_DIFF':'blue','AMBRA_KO_C11_DIFF':'green','AMBRA_KO_G11_DIFF':'purple'}

    # Plot mean with error shading of intensity-difference over time, for each protein
    if(make_plot):     # (have a flag here b/c plot takes a long time to make)
        sns.lineplot(data=df_full,
                     x='Time (s)',
                     y='intensity-diff',
                     hue='condition',
                     ax=axs[row_i][col_i],
                     palette=colors_dict)

        axs[row_i][col_i].set_title(protein_name)
        axs[row_i][col_i].set_ylabel('Intensity')

        # Plot KO cell lines in the final column (rows 0-3)
        if (protein_name in KO_Sphase_proteins):
            sns.lineplot(data=df_full_KO_Sphase, x='Time (s)', y='intensity-diff',
                         hue='condition',
                         ax=axs[row_i_][ncols - 1],
                         palette=colors_dict)

            axs[row_i_][ncols - 1].set_ylabel('Intensity')
            axs[row_i_][ncols - 1].set_title(f"{protein_name}")
            axs[row_i_][ncols - 1].set_xlim(0, 600)

            if (not protein_name == 'MSH6'):
                axs[row_i_][ncols - 1].get_legend().remove()


        for key in p_times.keys():
            if(p_times[key] != []):
                if(key == 'SIALLDIFF' or key == 'SIKIPsDIFF'):
                    axs[row_i][col_i].axvline(p_times[key][0], ls='--', color=colors_dict[key])
                    axs[row_i][col_i].axvline(p_times[key][1], ls='--', color=colors_dict[key])
                else:
                    axs[row_i_][ncols - 1].axvline(p_times[key][0],ls='--',color=colors_dict[key])
                    axs[row_i_][ncols - 1].axvline(p_times[key][1],ls='--',color=colors_dict[key])

        if (not protein_name == 'MSH6'):
            axs[row_i][col_i].get_legend().remove()

        if ((row_i + 1) % 7 == 0):
            row_i = 0
            col_i += 1
        else:
            row_i += 1

        if (protein_name in KO_Sphase_proteins):
            row_i_ += 1

# SAVE STATS, HEATMAP DATA, PLOTS
pvalue_df = pd.DataFrame(pvalue_matrix, columns=['protein',
                                                 'start_t',
                                                 'end_t',
                                                 'siAll_pvalue',
                                                 'siKIPs_pvalue',
                                                 'Ambra_KO_C11_pvalue',
                                                 'Ambra_KO_G11_pvalue',
                                                 'siCtrl_n',
                                                 'siAll_n',
                                                 'siKIPs_n',
                                                 'Ambra_WT_n',
                                                 'Ambra_KO_C11_n',
                                                 'Ambra_KO_G11_n'])

pvalue_df.to_csv(f"{root_dir}/all_pvalues-{p_value_window}.txt", sep='\t')

if(make_plot):
    fig.tight_layout()
    fig.savefig(f"{root_dir}/all_raw-{p_cutoff}-{p_value_window}.pdf")


