import numpy as np
import matplotlib.pyplot as plt
import glob as glob
import pandas as pd
import statistics as stats

######################## import & filter ########################
def import_samcov(file_path):
    ''' from samtools depth, per bp '''

    df = pd.read_csv(file_path, sep='\t', header=0)
    df.columns = ['Cx', 'Pos', 'Cov']

    return df

def subtract_bg(bp_array):
    ''' subtract median depth across entire library from all positions '''

    bg = bp_array.median().Cov
    print('background coverage: %s'%bg)
    bp_array.Cov = pd.Series(bp_array.Cov - bg)

    return bp_array, bg

def export_filtered_contigs_sam(df, contigs, file_name, date, save_dir):
    ''' save only the bp depth for passing contigs
        for downstream phage calling '''

    subdf = df[df.Cx.isin(contigs)]
    subdf.to_csv(save_dir+'%s_filtered_phage_bpcov_%s.csv' %(date, file_name))

    return

def import_filter_samcov(file_path, date, save_dir, depth_filter=50000, ext='.sorted.cov.txt'):
    ''' calculates the sum depth across all bp per contig (background-subtracted)
        remove contigs < depth sum
        whereby we assume min phage len 10k with 10X cov for induction, loosened to 50k sum
        from contigs that pass, export bp depth

        also saves summary of filtering steps '''

    file_name = file_path.split('/')[-1].split(ext)[0]

    df = import_samcov(file_path)
    df_, bg = subtract_bg(df)
    df_sum = df_.groupby(['Cx']).sum().reset_index()

    filtered = df_sum[df_sum.Cov >= depth_filter]
    filtered.to_csv(save_dir+'%s_filtered_depth_sum_%s_%s.csv' %(date, file_name, depth_filter))

    bg_info = pd.DataFrame({'library':file_name, 'bgCov': bg}, index=[0])
    bg_info.to_csv(save_dir+'%s_background_cov_%s.csv'%(date, file_name))

    export_filtered_contigs_sam(df, filtered.Cx, file_name, date, save_dir)

    return df

######################## enriched cov regions ########################
def check_enrichment(array, zero_check, zero_check_p, median_check, cov):
    ''' helper: identifies coordinates for putative phage induction region '''

    coord = 'hold'

    for n,i in enumerate(array):
        # stop at non-zero
        if i == 0:
            pass
        else:
            # is it a random read?
            if array[n:n+zero_check].count(0) >= zero_check*zero_check_p:
                pass
            else:
                # is the coverage reasonable?
                if stats.median(array[n:n+median_check]) >= cov:
                    coord = n
                    break
                else:
                    pass

    return coord

def identify_enriched_region(array, zero_check=50, zero_check_p=0.2,
                             median_check=500, cov=10, min_len=8000,
                             zero_ct_threshold=0.2, check_multi_len=100000):
    '''
    array: (list)

    Here, scan until non-zeros.
    Check proportion of zeros in next 50 bp, must be < 20% zeros.
    Then check what is the median of the next 500 bp.
    If median >= 10X, accept as start coordinate.

    Repeat from other end.
    Check that the length is minimum 10 kb --> 8 kb to be loose.
    Check that the region does not contain more than 20% zeroes.
    Notes which contigs to check for more than 1 phage in the region manually.
    Calculate the median coverage of the region between.

    '''
    # scan for coordinates of enriched region
    start = check_enrichment(array, zero_check, zero_check_p, median_check, cov)
    stop_ = check_enrichment(list(reversed(array)), zero_check, zero_check_p, median_check, cov)

    if start == 'hold' or stop_ == 'hold':
        start = 0
        stop_ = len(array)

    stop = len(array) - stop_
    phage_len = stop - start + 1
    phage_region = array[start:stop+1]

    # how sparse is the region & what is the median coverage
    zero_ct = list(phage_region).count(0)
    if len(phage_region) == 0:
        zero_portion = 0
        med_cov = 0
    else:
        zero_portion = zero_ct/len(phage_region)
        med_cov = stats.median(phage_region)

    # classify if phage
    if (phage_len >= min_len) and (med_cov >= cov) and (zero_portion < zero_ct_threshold):
        phage_check = True
    elif (phage_len >= check_multi_len):
        phage_check = 'CHECK' # see if contig coverage has multiple phages & recover manually
    else:
        phage_check = False

    start = start + 1 # coordinates start at 1

    return start, stop, phage_len, med_cov, zero_ct, zero_portion, phage_check

def add_enrich_summary(sum_df, bp_df, file_name, date, save_dir,
                       zero_check=50, zero_check_p=0.2, median_check=500,
                       cov=10, min_len=8000):
    ''' produces a .csv with the details of all passing contigs '''
    starts, stops, lens, covs, zero_cts, zero_p, checks, = [],[],[],[],[],[],[]

    for n,i in enumerate(sum_df.Cx):
        subdf = list(bp_df[bp_df.Cx.str.contains(i)].Cov)
        summary = identify_enriched_region(subdf, zero_check, zero_check_p,
                                           median_check, cov, min_len=8000)

        starts.append(summary[0])
        stops.append(summary[1])
        lens.append(summary[2])
        covs.append(summary[3])
        zero_cts.append(summary[4])
        zero_p.append(summary[5])
        checks.append(summary[6])

    sum_df['Start'] = starts
    sum_df['Stop'] = stops
    sum_df['Len'] = lens
    sum_df['MedCov'] = covs
    sum_df['ZeroCount'] = zero_cts
    sum_df['PortionZero'] = zero_p
    sum_df['IsPhage'] = checks

    sum_df.to_csv(save_dir+'%s_filtered_depth_sum_coords_%s.csv' %(date,file_name))

    return sum_df

######################## plotting ########################
def pull_phage_region(start, stop, bp_df, contig):
    ''' slices a df with cov of putative phage region '''

    subdf = bp_df[bp_df.Cx == contig]
    phage = subdf[(subdf.Pos >= start) & (subdf.Pos <= stop)]

    return phage

def plot_enriched_region(start, stop, bp_df, contig, file_name, date, save_dir,
                         clr, w, h, ylog_on, ylim_on, ylim_h, font_size):

    ''' Plots depth along entire region per bp '''

    phage = pull_phage_region(start, stop, bp_df, contig)

    fig,ax = plt.subplots(figsize=(w,h))
    ax.plot(phage.Pos, phage.Cov, color=clr)

    if ylog_on:
        ax.set_yscale('log')

    if ylim_on:
        ax.set_ylim(1, ylim_h)

    plt.title(file_name+': '+contig, size=font_size)
    plt.savefig(save_dir+'%s_phage_depth_%s_%s.png' %(date, file_name, contig), dpi=300)

    return phage

def plot_all_phages(sum_df, bp_df, file_name, date, save_dir, clr, w, h, ylog_on, ylim_on, ylim_h, font_size):

    ''' From a library, plot all called putative phages '''

    for n,i in enumerate(sum_df.Cx):
        phage = sum_df[sum_df.Cx == i]
        start = phage.Start.values[0]
        stop = phage.Stop.values[0]

        plot_enriched_region(start, stop, bp_df, i, file_name, date, save_dir,
                             clr, w, h, ylog_on, ylim_on, ylim_h, font_size)

    return

######################## combine steps ########################
def combine_calling(file_path, date, save_dir, clr, w, h, ylog_on, ylim_on, ylim_h,
                    font_size, ext='.sorted.cov.txt'):

    ''' NOTE: Hard-coded import files.
        All phage hit calling steps. '''

    name = file_path.split('/')[-1].split(ext)[0]
    print('Processing %s ...' %name)
    depth = import_filter_samcov(file_path, date, save_dir)
    filt = pd.read_csv(save_dir+'%s_filtered_depth_sum_%s_50000.csv' %(date, name), index_col=0)
    bpcov = pd.read_csv(save_dir+'%s_filtered_phage_bpcov_%s.csv' %(date, name), index_col=0)

    sum_df = add_enrich_summary(filt, bpcov, name, date, save_dir)
    phages = sum_df[sum_df.IsPhage == True]
    phages.to_csv(save_dir+'%s_called_phages_%s.csv'%(date, name))

    plot_all_phages(phages, bpcov, name, date, save_dir, clr,
                    w, h, ylog_on, ylim_on, ylim_h, font_size)

    return


def call_all(folder, date, save_dir, clr, w=8, h=2, ylog_on=False,
             ylim_on=False, ylim_h=5000, font_size=14, ext='.sorted.cov.txt'):

    ''' Run phage hit calling for all libraries '''

    all_depths = glob.glob(folder+'*%s'%ext)

    for n,i in enumerate(all_depths):
        print(i)
        combine_calling(i, date, save_dir, clr, w, h, ylog_on, ylim_on, ylim_h, font_size)

    return
