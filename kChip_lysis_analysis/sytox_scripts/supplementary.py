import pandas as pd
from sklearn.metrics import auc

# found in same folder
import sytox_scripts.bootstrap_and_z as bsz

################################# calculations #################################
def calc_dAUC(co_AUC, exp_AUC):
    '''
    Normalize score to a -1 to 1 scale whereby
    sign assignment:
    - coculture signal > expected signal is positive
    - coculture signal > expected signal is negative
    scaled to the higher AUC value
    Expected = sum of monocultures
    '''
    if co_AUC > exp_AUC:
        dAUC = (co_AUC - exp_AUC)/co_AUC # +
    elif exp_AUC > co_AUC:
        dAUC = (co_AUC - exp_AUC)/exp_AUC # -
    elif co_AUC == exp_AUC:
        dAUC = 0 # in case both are zeros

    return dAUC

def calc_dAUC_per_bs(df, co_col, mono_col1, mono_col2, col_name):
    '''
    Each row represents a different iteration
    to calculate with for whole distribution. '''

    mono_sum = df[mono_col1].add(df[mono_col2])
    df[col_name] = [calc_dAUC(df[co_col][n],mono_sum[n]) for n,i in enumerate(df.index)]

    return df

def standardize_dAUC_scores(df, score_col='dAUC_score', err_col='dAUC_score_err'):

    constant = df[score_col].median()
    df['dAUC_by_SE'] = df[score_col]/df[err_col]

    return df

def get_null_pop(df):
    return df[df['Left_simple'] == df['Right_simple']].reset_index(drop=True)

def get_exp_pop(df):
    return df[~(df['Left_simple'] == df['Right_simple'])].reset_index(drop=True)

################################## plotting ##################################
def calc_sns_errorbars(df, col_name, col_err):
    '''
    Returns two series representing the lower and upper values to fill
    error area around a lineplot (seaborn).

    Inputs:
        df: a dataframe containing both the main and error value series
            seaborn calls from same df
        col_name: name of column containing main data series
        col_err: name of column containing respective error data series
    '''
    lower = df[col_name] - df[col_err]
    upper = df[col_name] + df[col_err]

    return lower, upper

########################### dataframe manipulation  ###########################
def extract_combos(df, left, right):
    '''
    Slices df for specific combo (via contains substring)
    '''
    return df[df.Label_left.str.contains(left) & df.Label_right.str.contains(right)]

def summarize_single_combo(df, left, right, save_dir, save_desc=''):
    '''
    Export .csv with all of a specific combos' wells &
    the summarized RFUs ('norm2') per timepoint (mean, std, median, sem).

    df: source dataframe with all droplets for all combos at all timepoints
    l/r: (str) left & right label to extract
    save_prefix: (str) beginning of save path for export
        e.g. save_path+'_media_media'
    '''
    subdf = extract_combos(df, left, right)
    tp_col = subdf.filter(like='norm2').columns
    df_sem = pd.DataFrame(bsz.boot_microwells(subdf.filter(like='norm2'))[1])
    df_sem.index = subdf.filter(like='norm2').columns
    summ_values = [subdf.mean(axis=0), subdf.std(axis=0), subdf.median(axis=0), df_sem]
    columns=['mean', 'std', 'median', 'sem']

    summ_df = pd.concat(summ_values, axis=1)
    summ_df.columns = columns

    # timepoints as column for easier downstream manipulation
    summ_df.T.to_csv(save_dir+'summarized'+save_desc+'.csv')

    return

def pull_overlap_df(source_df, ref_df, col_name1, col_name2):
    '''
    Returns a sliced source_df given matching values in specified df columns.
    '''
    return source_df[source_df[col_name1].isin(ref_df[col_name2])].reset_index(drop=True)

def remove_overlap_df(source_df, ref_df, col_name1, col_name2):
    '''
    Returns a sliced source_df removing rows with matching values in specified df columns.
    '''
    return source_df[~source_df[col_name1].isin(ref_df[col_name2])].reset_index(drop=True)

def filter_col_min(df, col, threshold):
    ''' Return a df filtered by a minimum threshold. '''
    filter_df = df[df[col] >= threshold]
    print(filter_df.shape)
    return filter_df.reset_index(drop=True)

def filter_col_equal(df, col, threshold):
    ''' Return a df filtered by an exact value. '''
    filter_df = df[df[col] == threshold]
    print(filter_df.shape)
    return filter_df.reset_index(drop=True)

def filter_col_max(df, col, threshold):
    ''' Return a df filtered by an exact value. '''
    filter_df = df[df[col] <= threshold]
    print(filter_df.shape)
    return filter_df.reset_index(drop=True)

############################## dataframe labels  ##############################
def replace_label(df, old, new):
    if old in df.values:
        return df.applymap(lambda x: new if x == old else x)
    else:
        return df

def replace_all_labels(df, label_map, from_col, to_col):
    for n,i in enumerate(label_map.index):
        df = replace_label(df, label_map[from_col][n], label_map[to_col][n])
    return df

def correct_labels(df, label_map, save_dir, og_col='original', temp_col='temp', new_col='corrected'):
    corrected = replace_all_labels(df, label_map, og_col, temp_col)
    corrected = replace_all_labels(df, label_map, temp_col, new_col)

    corrected.to_csv(save_dir+'distance_and_area_filtered.csv')

    return corrected

def concatenate_labels(df, col1='Label_left', col2='Label_right', combo_col='Labels_combo'):
    '''
    Returns a df which concatenates the left & right droplet labels.
    Used to easily pull information between dfs.'''

    separator_list = ['_']*len(df[col1])
    df[combo_col] = df[col1]+separator_list+df[col2]

    return df

def trim_label_name(label, split_by='_', multi=False, position1=1, position2=3):
    '''
    Returns a reformated barcode label name into a simplier, readable verison.
    e.g. 'BUG_12B09_M9' to '12B09' for cleaner plot labeling

    Inputs:
        label : (str) barcode name
        split_by : (str) separator e.g. '_'
        multi : (boolean) if keeping multiple slices (must be continuous)
        position1 & 2 : indices if taking multiple slices
    '''
    if multi == True:
        label = '_'.join(label.split(split_by)[position1:position2])
    else:
        label = label.split(split_by)[position1]
    return label

def get_unique_colabels(df):
    subdf = df.groupby(['Label_left', 'Label_right']).count().reset_index()
    return pd.unique(subdf[['Label_left', 'Label_right']])

#############################################################################
