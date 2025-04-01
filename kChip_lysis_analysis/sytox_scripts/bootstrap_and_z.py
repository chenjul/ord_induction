import pandas as pd
import numpy as np
import multiprocessing as mp
from pathos.multiprocessing import ProcessingPool as Pool
from sklearn.metrics import auc

# same folder
import sytox_scripts.supplementary as helper

######### bootstrapping #########
def boot(ar, ar_size, parallel=False):
    '''
    Returns the standard deviation of the distribution from bootstrapped values
    or the dataframe of the distribution if parallelized.
    '''
    # resample from label column
    if parallel:
        return pd.DataFrame(np.random.choice(ar, ar_size))
    else:
        return np.std(np.random.choice(ar, ar_size))

def boot_parallel(arglist):
    '''
    Input:
    arglist: 1D array (list) with values for bootstrapping

    * unpacks arglist and converts from list to pd.Series
    Required for ar.shape[0] in boot
    '''
    return boot(*arglist)

def boot_array(ar, ar_size_default=True, given_size='', bs_size=1000):
    '''
    Returns the standard error of the median (std of median distribution)
    Bootstraps the input array 1000 times (parallelized)
    Generates a distribution of medians from the 1000 bootstraps

    Input:
    ar : 1D array with values for bootstrapping
    Note [ar, i] makes parallel=True for boot
    '''
    cores = mp.cpu_count() # because we usually run this across four datasets
    pool = Pool(cores)

    if ar_size_default == False:
        ar_size = given_size
    else:
        ar_size = ar.shape[0]

    # process the DataFrame by mapping function to each df across the pool
    boot_ = pool.map(boot_parallel,([ar,ar_size,i] for i in range(1,bs_size+1)))

    # close down the pool and join
    pool.close()
    pool.join()
    pool.clear()

    booted = pd.concat(boot_, axis=1)
    return booted

def boot_microwells(df, ar_size_default=True, given_size='', bs_size=1000):
    '''
    df: subdf for specific combination, RFU/tp only (no Labels)
    - RFU/tp col only, labels removed
    - each row represents a microwell
    '''
    ar = df.to_numpy()

    if ar_size_default == False:
        booted = boot_np_array(ar, ar_size_default=False, given_size=given_size, bs_size=bs_size)
    else:
        booted = boot_np_array(ar, bs_size=bs_size)
    booted = booted.T

    # mean & SEMs of median distribution across tp for combo
    means = booted.mean().values
    errors = booted.std().values

    # all bs AUCs for combination
    num_tp = df.shape[1]
    tps = list(range(num_tp))
    bs_aucs = [auc(tps, booted.iloc[n]) for n,i in enumerate(booted.index)]

    return means, errors, bs_aucs

def boot_groupby(df, col1='Label_left', col2='Label_right', rfu_label='norm'):
    '''
    Returns the standard error of the median (from bootstrapping),
    for every timepoint for each specific combination.
    '''
    col_names = df.columns
    med_dict = {}
    sem_dict = {}
    auc_dict = {}

    for labels, group in df.groupby((col1, col2)):
        tps_rfu = group.filter(like=rfu_label).reset_index(drop=True)
        bs = boot_microwells(tps_rfu)
        med_dict[labels] = bs[0] # est median
        sem_dict[labels] = bs[1] # est error
        auc_dict[labels] = bs[2]

    df_med = pd.DataFrame(med_dict).T.reset_index()
    df_sem = pd.DataFrame(sem_dict).T.reset_index()
    df_med.columns = col_names
    df_sem.columns = col_names

    bs_col = ['Label_left', 'Label_right']+['bs%d' %i for i in range(1000)]
    df_auc = pd.DataFrame(auc_dict).T.reset_index()
    df_auc.columns = bs_col

    return df_med, df_sem, df_auc

######### Z' analysis #########
def calc_z(pos_m ,neg_m, pos_err, neg_err):
    '''
    Calculates the Z' between the positive & negative controls (combos).

    pos_m/neg_m: mean or median of positive/negative controls
    pos_err/neg_err: standard deviation or standard error of median for controls
    '''
    z = 1 - (3*(pos_err + neg_err) / (abs(pos_m - neg_m)))
    return z
