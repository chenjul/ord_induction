import argparse
import time
import pandas as pd
import numpy as np

import warnings
warnings.filterwarnings('ignore')

# assumes in same folder as other sytox_scripts
import sytox_scripts.bootstrap_and_z as bsz
import sytox_scripts.supplementary as helper

#############################################################################
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csvfile",
                        help="path to .csv file with all droplets, often blanked",
                        type=str)
    parser.add_argument("--save_dir",
                        help="save path for exported summary .csvs",
                        type=str)
    parser.add_argument("--chip_id",
                        help="initial part to file name",
                        type=str)
    parser.add_argument("--save_desc",
                        help="optional description to file names",
                        default="",
                        type=str)
    parser.add_argument("--mono_key",
                        help="substring to find culture droplets",
                        default="BUG",
                        type=str)
    parser.add_argument("--media_key",
                        help="substring to find media droplets",
                        default="MEDIA",
                        type=str)
    args = parser.parse_args()
    return args
#############################################################################
def main(csvfile, save_dir, chip_id, save_desc, mono_key, media_key):
    '''
    For every combination, generate its summary statistics
    (mean, median, std, SEM) and save the bootstrapped AUCs.
    '''
    # import main csv
    df=pd.read_csv(csvfile, index_col=0)
    path = save_dir+chip_id

    # summarize data
    subdf_mean = df.groupby(['Label_left', 'Label_right']).mean().reset_index()
    subdf_std = df.groupby(['Label_left', 'Label_right']).std().reset_index()
    subdf_med = df.groupby(['Label_left', 'Label_right']).median().reset_index()
    subdf_boot = bsz.boot_groupby(df)

    # export summarized & calculated values
    subdf_mean.to_csv(path+'summarized_combos'+save_desc+'_mean.csv')
    subdf_std.to_csv(path+'summarized_combos'+save_desc+'_std.csv')
    subdf_med.to_csv(path+'summarized_combos'+save_desc+'_med.csv')
    subdf_boot[0].to_csv(path+'summarized_combos'+save_desc+'_estmed.csv')
    subdf_boot[1].to_csv(path+'summarized_combos'+save_desc+'_sem.csv')
    subdf_boot[2].to_csv(path+'summarized_combos'+save_desc+'_bs_aucs.csv')

if __name__=="__main__":

    start_time = time.time()
    args = get_args().__dict__

    main(**args)
    print("--- %s seconds ---" % (time.time() - start_time))
