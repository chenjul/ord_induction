import pandas as pd
import numpy as np
import scipy.stats as stats
import itertools
from sklearn.metrics import auc

# same folder
import sytox_scripts.supplementary as helper

######################## single coculture combination ########################
def calc_coculture_folddiff(co_df, add_df, fold_col='fold_dif'):
    '''
    Applies to a single combination/coculture.
    Note: returns transposed df
    '''
    fd = pd.concat([co_df, add_df]).T
    fd[fold_col] = fd.iloc[:,0]/fd.iloc[:,1]
    fd.columns = ['coculture', 'additive', 'fold_dif']

    return fd

def calc_mono_sum(left, right, mono_df, bug_col='bug_label', rfu_col='RFU', err_col='error'):
    '''
    Helper function in compiling the master, reformatted coculture data.
    Calculates the expected signal - sum of respective monocultures from a coculture.

    Inputs:
        left/right: (str) ID strains from coculture
        mono_df: imported .csv from pull_monocultures function
        bug_col: (str) column name in mono_df that stores left/right labels
        rfu_col: (str) column name in mono_df that stores monoculture summarized RFU
    '''
    subdf_left = mono_df[mono_df[bug_col].str.contains(left)].reset_index()
    subdf_right = mono_df[mono_df[bug_col].str.contains(right)].reset_index()

    left_rfu = subdf_left[rfu_col]
    right_rfu = subdf_right[rfu_col]
    left_err = subdf_left[err_col]
    right_err = subdf_right[err_col]
    mono_sum = left_rfu.add(right_rfu)

    return mono_sum, left_rfu, right_rfu, left_err, right_err

##################  reformat chip df to compile mono/co data ##################
def identify_mono_labels(df, media_key='MEDIA'):
    '''
    Identify the label order for monoculture combos.
    Set as global variable.
    '''
    if any(media_key in x for x in df.Label_right) == True:
        bug_label = 'Label_left'
        media_label = 'Label_right'
    else:
        bug_label = 'Label_right'
        media_label = 'Label_left'

    return bug_label, media_label

def pull_monocultures(df, df_err, save_dir, media_label, bug_label,
                      mono_key='BUG', media_key='MEDIA', save_desc=''):
    '''
    Returns a .csv with the summary statistic & respective error
    for all monocultures (bug/media) at all timepoints.

    Columns: timepoint, bug label, rfu, error
    Each set of combo data concatenated downwards.

    Inputs:
        df: imported .csv of mean or med
        df_err: imported .csv of std or sem
        mono_key: (str) common key in bacteria droplet inputs
        media_key: (str) common key in media droplet inputs
        media_ & bug_label: (str) ordered labels given the culture/media label format -- set globally
        save_dir: (str) destination path
        save_desc: (str) optional info appended to exported filenames
    '''

    ### relabel columns for pivot ###
    df = df.copy()
    df_err = df_err.copy()
    df.columns = ['Label_left', 'Label_right']+ ['RFU%d' %i for i in range(df.shape[1]-2)]
    df_err.columns = ['Label_left', 'Label_right']+['error%d' %i for i in range(df.shape[1]-2)]

    ### separate monoculture data from all, recognizing label order ###
    if media_label == 'Label_right':
        mono_df = helper.extract_combos(df, mono_key, media_key).reset_index(drop=True)
        mono_err = helper.extract_combos(df_err, mono_key, media_key).reset_index(drop=True)
    else:
        mono_df = helper.extract_combos(df, media_key, mono_key).reset_index(drop=True)
        mono_err = helper.extract_combos(df_err, media_key, mono_key).reset_index(drop=True)

    ### compile kinetic data ###
    pivot_rfu = pd.wide_to_long(mono_df, stubnames='RFU', j='tp',
                                i=[bug_label, media_label]).reset_index()
    pivot_err = pd.wide_to_long(mono_err, stubnames='error', j='tp',
                                i=[bug_label, media_label]).reset_index()
    full_mono_df = pivot_rfu.merge(pivot_err, on=['Label_left', 'Label_right', 'tp'])

    full_mono_df = full_mono_df.rename(columns={bug_label:'bug_label'})[['tp', 'bug_label', 'RFU', 'error']]
    full_mono_df['RFU_werr'] = full_mono_df.RFU.add(full_mono_df.error)

    full_mono_df.to_csv(save_dir+'summarized_monocultures'+save_desc+'.csv')

    return full_mono_df

def summarize_bs_dAUCs(df, left, right, media_label,
                       media_key='MEDIA'):
    ''' Returns summary data from bootstrapped distribution of combination. '''
    # coculture
    sub_aucs = helper.extract_combos(df, left, right).filter(like='bs').iloc[0]

    # monocultures
    if media_label == 'Label_right':
        sub_aucs_monoL = helper.extract_combos(df, left, media_key).filter(like='bs').iloc[0]
        sub_aucs_monoR = helper.extract_combos(df, right, media_key).filter(like='bs').iloc[0]
    else:
        sub_aucs_monoL = helper.extract_combos(df, media_key, left).filter(like='bs').iloc[0]
        sub_aucs_monoR = helper.extract_combos(df, media_key, right).filter(like='bs').iloc[0]

    co_AUC = np.mean(sub_aucs)
    L_AUC = np.mean(sub_aucs_monoL)
    R_AUC = np.mean(sub_aucs_monoR)

    subdf = pd.DataFrame({'co_AUC': sub_aucs, 'L_AUC': sub_aucs_monoL, 'R_AUC': sub_aucs_monoR})
    dAUC_df = helper.calc_dAUC_per_bs(subdf, 'co_AUC', 'L_AUC', 'R_AUC', 'dAUC_score')

    dAUC_df['Label_left'] = [left]*dAUC_df.shape[0]
    dAUC_df['Label_right'] = [right]*dAUC_df.shape[0]
    dAUC_mean = np.nanmean(dAUC_df.dAUC_score)
    dAUC_std = np.nanstd(dAUC_df.dAUC_score)
    dAUC_df = dAUC_df.reset_index()

    return dAUC_mean, dAUC_std, dAUC_df[['Label_left', 'Label_right', 'dAUC_score']], co_AUC, L_AUC, R_AUC

def pull_cocultures(df, df_err, mono_df, df_auc, chip_id, save_dir,
                    media_label, mono_key='BUG', save_desc=''):
    '''
    Returns a .csv with the summary statistic & respective error
    for all coculture combinations at all timepoints.

    Columns:
        timepoint, labels left&right, coculture rfu & error,
        calculated mono sum rfu & plotting_error,
        respective monoculture (left/right) rfu & error,

        AUCs for coculture, mono sum, left/right
        fraction (in decimal) of mono sum, left, or right AUC to co AUC
        calculated fold-change between actual (coculture) & expected (mono sum)

    Each set of combo data concatenated downwards.

    Inputs:
        df: imported .csv of mean or med
        df_err: imported .csv of std or sem
        df_auc: imported .csv of bootstrapped AUCs
        mono_df: imported .csv from pull_monocultures function
        mono_key: (str) common key in bacteria droplet inputs
        save_dir: (str) destination path
        save_desc: (str) optional info appended to exported filenames
    '''
    ### separate coculture data from all ###
    co_df = helper.extract_combos(df, mono_key, mono_key).reset_index()
    co_err = helper.extract_combos(df_err, mono_key, mono_key).reset_index()
    kinetic_df = co_df.filter(like='norm2')
    kinetic_err = co_err.filter(like='norm2')

    ### format compiled df ###
    chip_id_l, chip_id_s = [], []
    labels_left, labels_right, rfu, rfu_adjust, error = [], [], [], [], []
    left_rfu, left_adjust, left_err = [], [], []
    right_rfu, right_adjust, right_err = [], [], []
    mono_sum, sum_adjust = [], []
    peak_rfu_co, peak_rfu_left, peak_rfu_right = [], [], []
    peak_rfu_co_adj, peak_rfu_left_adj, peak_rfu_right_adj = [], [], []
    dAUC_scores_emp, dAUC_scores_emp_adj = [], []
    dAUC_est, dAUC_err, dAUC_adj = [], [], []
    q_25, q_75 = [], []
    co_auc_all, left_auc_all, right_auc_all = [], [], []

    num_combos = co_df.shape[0]
    num_tp = len(kinetic_df.columns)
    single_culture_tp = list(range(num_tp))
    all_culture_tps = single_culture_tp*num_combos # tp col values
    sum_err = [0]*num_combos*num_tp # for plotting purposes, no actual have error bars

    ### save separate df of bs dAUCs ###
    bs_dAUCs = []

    ### save dAUC for pval - non-kinetic ###
    labels_left_s, labels_right_s = [], []
    dAUC_scores_emp_s, dAUC_scores_emp_adj_s = [], []
    dAUC_s, dAUC_err_s, dAUC_adj_s = [], [], []
    q_25s, q_75s = [], []
    co_auc_s, left_auc_s, right_auc_s = [], [], []
    co_peak_s, left_peak_s, right_peak_s = [], [], []
    co_peak_s_adj, left_peak_s_adj, right_peak_s_adj = [], [], []

    ### per coculture ###
    # bug labels produced in the same order between df & df_error
    for n in co_df.index:
        left = co_df.Label_left[n]
        right = co_df.Label_right[n]
        chip_id_l.extend(itertools.repeat(chip_id, num_tp))
        chip_id_s.append(chip_id)
        labels_left.extend(itertools.repeat(left, num_tp))
        labels_right.extend(itertools.repeat(right, num_tp))
        labels_left_s.append(left)
        labels_right_s.append(right)

        # coculture
        rfu.extend(kinetic_df.iloc[n])
        error.extend(kinetic_err.iloc[n])

        # monoculture & mono sum RFU signal
        sum_and_mono = calc_mono_sum(left, right, mono_df)
        sum_and_mono2 = calc_mono_sum(left, right, mono_df, rfu_col='RFU_werr')

        mono_sum.extend(sum_and_mono[0])
        left_rfu.extend(sum_and_mono[1])
        right_rfu.extend(sum_and_mono[2])
        left_err.extend(sum_and_mono[3])
        right_err.extend(sum_and_mono[4])

        rfu_adjust.extend(kinetic_df.iloc[n] - kinetic_err.iloc[n])
        sum_adjust.extend(sum_and_mono2[0])
        left_adjust.extend(sum_and_mono2[1])
        right_adjust.extend(sum_and_mono2[2])

        # proxy for sustained, saturated RFU & RFU highpass filter
        peak_rfu_co.extend(itertools.repeat(kinetic_df.iloc[n].max(), num_tp))
        peak_rfu_left.extend(itertools.repeat(sum_and_mono[1].max(), num_tp))
        peak_rfu_right.extend(itertools.repeat(sum_and_mono[2].max(), num_tp))
        co_peak_s.append(kinetic_df.iloc[n].max())
        left_peak_s.append(sum_and_mono[1].max())
        right_peak_s.append(sum_and_mono[2].max())

        peak_rfu_co_adj.extend(itertools.repeat((kinetic_df.iloc[n] - kinetic_err.iloc[n]).max(), num_tp))
        peak_rfu_left_adj.extend(itertools.repeat(sum_and_mono2[1].max(), num_tp))
        peak_rfu_right_adj.extend(itertools.repeat(sum_and_mono2[2].max(), num_tp))
        co_peak_s_adj.append((kinetic_df.iloc[n] - kinetic_err.iloc[n]).max())
        left_peak_s_adj.append(sum_and_mono2[1].max())
        right_peak_s_adj.append(sum_and_mono2[2].max())

        # AUCs/dAUCs -- bootstrapped with estimated median dAUC (vs. empirical median)
        dAUCs = summarize_bs_dAUCs(df_auc, left, right, media_label=media_label)
        quants = np.nanquantile(dAUCs[2].dAUC_score, [0.25,0.75])
        dAUC_adj_score = dAUCs[0] - dAUCs[1]

        dAUC_est.extend(itertools.repeat(dAUCs[0], num_tp))
        dAUC_err.extend(itertools.repeat(dAUCs[1], num_tp))
        dAUC_adj.extend(itertools.repeat(dAUC_adj_score, num_tp))
        q_25.extend(itertools.repeat(quants[0], num_tp))
        q_75.extend(itertools.repeat(quants[1], num_tp))
        bs_dAUCs.append(dAUCs[2])
        co_auc_all.extend(itertools.repeat(dAUCs[3], num_tp))
        left_auc_all.extend(itertools.repeat(dAUCs[4], num_tp))
        right_auc_all.extend(itertools.repeat(dAUCs[5], num_tp))

        dAUC_s.append(dAUCs[0])
        dAUC_err_s.append(dAUCs[1])
        dAUC_adj_s.append(dAUC_adj_score)
        q_25s.append(quants[0])
        q_75s.append(quants[1])
        co_auc_s.append(dAUCs[3])
        left_auc_s.append(dAUCs[4])
        right_auc_s.append(dAUCs[5])

        # AUCs/dAUCs -- empirical
        co_AUC_emp = auc(single_culture_tp, kinetic_df.iloc[n])
        left_AUC_emp = auc(single_culture_tp, sum_and_mono[1])
        right_AUC_emp = auc(single_culture_tp, sum_and_mono[2])
        mono_sum_AUC_emp = left_AUC_emp + right_AUC_emp

        co_AUC_emp_adj = auc(single_culture_tp, (kinetic_df.iloc[n] - kinetic_err.iloc[n]))
        left_AUC_emp_adj = auc(single_culture_tp, sum_and_mono2[1])
        right_AUC_emp_adj = auc(single_culture_tp, sum_and_mono2[2])
        mono_sum_AUC_emp_adj = left_AUC_emp_adj + right_AUC_emp_adj

        dAUC_emp = helper.calc_dAUC(co_AUC_emp, mono_sum_AUC_emp)
        dAUC_emp_adj = helper.calc_dAUC(co_AUC_emp_adj, mono_sum_AUC_emp_adj)

        dAUC_scores_emp.extend(itertools.repeat(dAUC_emp, num_tp))
        dAUC_scores_emp_adj.extend(itertools.repeat(dAUC_emp_adj, num_tp))
        dAUC_scores_emp_s.append(dAUC_emp)
        dAUC_scores_emp_adj_s.append(dAUC_emp_adj)

    ### final dataframes ###
    bs_dAUC_df = pd.concat(bs_dAUCs)

    full_df = pd.DataFrame({'chip_ID': chip_id_l, 'tp': all_culture_tps,
                            'Label_left': labels_left, 'Label_right': labels_right,
                            'dAUC_score_emp': dAUC_scores_emp, 'dAUC_score_emp_adj': dAUC_scores_emp_adj,
                            'dAUC_score': dAUC_est, 'dAUC_score_err': dAUC_err, 'dAUC_score_adj': dAUC_adj,
                            'co_RFU': rfu, 'co_adj': rfu_adjust, 'co_error': error,
                            'mono_sum_RFU': mono_sum, 'mono_sum_RFU_adj': sum_adjust, 'mono_sum_error': sum_err,
                            'left_RFU': left_rfu, 'left_adj': left_adjust, 'left_error': left_err,
                            'right_RFU': right_rfu, 'right_adj': right_adjust, 'right_error': right_err,
                            'co_AUC': co_auc_all, 'left_AUC': left_auc_all, 'right_AUC': right_auc_all,
                            'co_peak_RFU': peak_rfu_co, 'left_peak_RFU': peak_rfu_left, 'right_peak_RFU': peak_rfu_right,
                            'co_peak_RFU_adj': peak_rfu_co_adj,
                            'left_peak_RFU_adj': peak_rfu_left_adj, 'right_peak_RFU_adj': peak_rfu_right_adj,
                            'dAUC_q_25': q_25, 'dAUC_q_75': q_75})

    ### calculate & append fold-change values ###
    full_df['fold_change'] = full_df['co_RFU']/full_df['mono_sum_RFU']
    full_df['RFU_diff'] = full_df['co_RFU'] - full_df['mono_sum_RFU']

    full_df['fold_change_adj'] = full_df['co_adj']/full_df['mono_sum_RFU_adj']
    full_df['RFU_diff_adj'] = full_df['co_adj'] - full_df['mono_sum_RFU_adj']

    # to help group by timepoint (DNA extraction)
    # "normalize" to max coculture RFU
    p_, p_adj = [], []

    for n,i in enumerate(full_df['RFU_diff']):
        if peak_rfu_co[n] != 0:
            p_.append(i/peak_rfu_co[n])
        else:
            p_.append(0)

    for n,i in enumerate(full_df['RFU_diff_adj']):
        if peak_rfu_co_adj[n] != 0:
            p_adj.append(i/peak_rfu_co_adj[n])
        else:
            p_adj.append(0)

    full_df['p_RFU_diff'] = p_
    full_df['p_RFU_diff_adj'] = p_adj

    non_kinetic = pd.DataFrame({'chip_ID': chip_id_s, 'Label_left': labels_left_s, 'Label_right': labels_right_s,
                                'dAUC_score_emp': dAUC_scores_emp_s, 'dAUC_score_emp_adj': dAUC_scores_emp_adj_s,
                                'dAUC_score': dAUC_s, 'dAUC_score_err': dAUC_err_s,'dAUC_score_adj': dAUC_adj_s,
                                'co_AUC': co_auc_s, 'left_AUC': left_auc_s, 'right_AUC': right_auc_s,
                                'co_peak_RFU': co_peak_s, 'left_peak_RFU': left_peak_s, 'right_peak_RFU': right_peak_s,
                                'co_peak_RFU_adj': co_peak_s_adj,
                                'left_peak_RFU_adj': left_peak_s_adj, 'right_peak_RFU_adj': right_peak_s_adj,
                                'dAUC_q_25': q_25s, 'dAUC_q_75': q_75s})

    # labels used downstream
    full_df['Left_simple'] = [helper.trim_label_name(n) for n in full_df.Label_left]
    full_df['Right_simple'] = [helper.trim_label_name(n) for n in full_df.Label_right]
    non_kinetic['Left_simple'] = [helper.trim_label_name(n) for n in non_kinetic.Label_left]
    non_kinetic['Right_simple'] = [helper.trim_label_name(n) for n in non_kinetic.Label_right]

    final_df = helper.concatenate_labels(full_df, col1='Label_left', col2='Label_right', combo_col='Labels_combo')
    final_df = helper.concatenate_labels(full_df, col1='Left_simple', col2='Right_simple', combo_col='Combo_simple')
    final_df = helper.concatenate_labels(full_df, col1='chip_ID', col2='Labels_combo', combo_col='full_ID')
    non_kinetic = helper.concatenate_labels(non_kinetic, col1='Label_left', col2='Label_right', combo_col='Labels_combo')
    non_kinetic = helper.concatenate_labels(non_kinetic, col1='Left_simple', col2='Right_simple', combo_col='Combo_simple')
    non_kinetic = helper.concatenate_labels(non_kinetic, col1='chip_ID', col2='Labels_combo', combo_col='full_ID')

    final_df.to_csv(save_dir+'summarized_cocultures_kinetic'+save_desc+'.csv')
    bs_dAUC_df.to_csv(save_dir+'bootstrapped_dAUC_scores'+save_desc+'.csv')
    non_kinetic.to_csv(save_dir+'summarized_cocultures_nonkinetic'+save_desc+'.csv')

    return final_df.head()
