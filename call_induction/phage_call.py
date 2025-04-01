import argparse
import time
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import glob as glob
import statistics as stats

import warnings
warnings.filterwarnings('ignore')

# same folder
import phage_call_functions as call

#############################################################################
def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--date",
                        help='YYYYMMDD',
                        type=str)
    parser.add_argument("--source_dir",
                        help='directory to samtools depth file',
                        type=str)
    parser.add_argument("--save_dir",
                        help="save path for exported summary .csvs",
                        type=str)
    parser.add_argument("--clr",
                        help="plotting colour",
                        default="#DAC6D4",
                        type=str)

    args = parser.parse_args()
    return args
#############################################################################
def main(date, source_dir, save_dir, clr):

    ''' For every library (sorted.cov.txt file),
        call putative prophage induction regions. '''

    call.call_all(folder=source_dir, save_dir=save_dir, date=date, clr=clr)

if __name__=="__main__":

    start_time = time.time()
    args = get_args().__dict__

    main(**args)
    print("--- %s seconds ---" % (time.time() - start_time))
