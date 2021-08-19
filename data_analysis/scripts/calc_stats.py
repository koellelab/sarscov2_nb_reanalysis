import argparse
import pandas as pd
import numpy as np
import glob
import itertools
import os
import scipy.stats


def process_vcf(vcf_path, filter_allow):
    sample_name = \
        'CoV_' + vcf_path.split('/')[-1].split('_')[1]
    with open(vcf_path, 'r') as fp:
        for line in fp:
            if line[0:2] != '##':
                header=line.split('\n')[0].split('\t')
                break
    vcf = \
        pd.read_csv(vcf_path, sep='\t', header=None, comment='#')
    vcf.columns = header
    vcf['SAMPLE'] = sample_name
    vcf['FILTER_PASS'] = \
        vcf['FILTER'].isin(filter_allow)
    # todo make this better
    vcf['AF'] = \
        vcf.loc[:,'INFO'].str.split(';').apply(lambda k: 
            [i.split('=')[1] for i in k if i.split('=')[0]=='AF'][0]) 
    vcf['DP4'] = \
        vcf.loc[:,'INFO'].str.split(';').apply(lambda k: 
            [sum([int(j) for j in i.split('=')[1].split(',')]) 
            for i in k if i.split('=')[0]=='DP4'][0])
    vcf['AF'] = vcf.loc[:,'AF'].astype(float)
    vcf['TYPE'] = \
        np.where((vcf['REF'].str.len()==1) & 
            (vcf['REF'].str.len()==1), 'SNP','INDEL')
    vcf = vcf[['SAMPLE', '#CHROM', 
        'POS','REF', 'ALT', 'AF', 'DP4', 
        'FILTER_PASS', 'TYPE']]
    return(sample_name, vcf)


def get_transmission_pairs(dat):
    import ast
    dat['donor_pairs'] = \
        dat['donor_pairs'].apply(lambda k: ast.literal_eval(k))
    dat['recip_pairs'] = \
        dat['recip_pairs'].apply(lambda k: ast.literal_eval(k))
    # key is transmission pair ID, value is donor 
    donor_dict = {}
    # key is transmission pair ID, value is recipient
    recip_dict = {}
    for idx, row in dat.iterrows():
        for item in row['donor_pairs']:
            donor_dict[item] = row['sample_name']
        for item in row['recip_pairs']:
            recip_dict[item] = row['sample_name']
    if set(donor_dict.keys()) != set(recip_dict.keys()):
        raise Exception('incomplete transmission pair data in metadata')
    pair_dict = {}
    for key,value in donor_dict.items():
        pair_dict[key] = (value, recip_dict[key])
    return(pair_dict)


def filter_pair_dat(pair_dat, only_shared=False, min_af=0.00):
    filter_pair_dat = pair_dat.copy()
    # exclude any that failed filter
    filter_pair_dat = \
        filter_pair_dat[(filter_pair_dat['FILTER_PASS_x'] != False) & 
            (filter_pair_dat['FILTER_PASS_y'] != False)]
    # set anything to 0/1 that is above/below a certain threshold
    filter_pair_dat.loc[filter_pair_dat['AF_x'] < min_af, 'AF_x'] = 0
    filter_pair_dat.loc[filter_pair_dat['AF_y'] < min_af, 'AF_y'] = 0
    filter_pair_dat.loc[filter_pair_dat['AF_x'] > 1-min_af, 'AF_x'] = 1
    filter_pair_dat.loc[filter_pair_dat['AF_y'] > 1-min_af, 'AF_y'] = 1
    # only shared variants
    if only_shared==True:
        # when min_af is > 0 we want to include min_af
        # to be consistent with above
        # when it is 0 though we need to make sure not to count that 
        # as shared
        filter_pair_dat = \
            filter_pair_dat[(filter_pair_dat['AF_x'] >= min_af) & 
                (filter_pair_dat['AF_x'] > 0) &
                (filter_pair_dat['AF_y'] >= min_af) & 
                (filter_pair_dat['AF_y'] > 0)]
    return(filter_pair_dat)


def get_pair_dat(dat_1, dat_2):
    pair_dat = dat_1.merge(dat_2, 
        left_on=['POS', 'REF', 'ALT'], 
        right_on=['POS', 'REF', 'ALT'], 
        how='outer')
    pair_dat.loc[:,['AF_x', 'AF_y']] = \
        pair_dat.loc[:,['AF_x', 'AF_y']].fillna(0)
    pair_dat.loc[:,['AF_x', 'AF_y']] = \
        pair_dat.loc[:,['AF_x', 'AF_y']].astype(float)
    return(pair_dat)


def prop_shared_test(vcf_dat=None, transmission_pairs=None, random_pairs=None, min_af=0.0):
    def calc_prop_transmitted_binned(samples, bins, min_af=0.0):
            # each bin includes the data at that bin value and
            # up until the next
            donor_dat = vcf_dat[samples[0]]
            donor_dat = donor_dat[donor_dat['FILTER_PASS'] == True]
            donor_af_binned = np.histogram(donor_dat['AF'], bins=bins)[0]
            samples_dat = get_pair_dat(vcf_dat[samples[0]], 
                vcf_dat[samples[1]])
            samples_dat = \
                filter_pair_dat(samples_dat, 
                    min_af=min_af, only_shared=True)
            transmitted_af_binned = np.histogram(samples_dat['AF_x'], bins=bins)[0]
            with np.errstate(divide='ignore', invalid='ignore'):
                prop_transmitted_binned = transmitted_af_binned/donor_af_binned
            prop_transmitted_binned = np.nan_to_num(prop_transmitted_binned)
            return(prop_transmitted_binned)
    bins = np.array([0, 0.01, 0.02, 0.03, 0.06])
    prop_transmitted_pairs = np.zeros((len(transmission_pairs), bins.size-1)) 
    for pair_idx, pair in enumerate(transmission_pairs.values()):
        # calculate the proportion of transmitted variants in each 
        # bin for this pair
        prop_transmitted_pairs[pair_idx,:] = \
            calc_prop_transmitted_binned(pair, bins, min_af=min_af)
    prop_transmitted_non_pairs = np.zeros((len(transmission_pairs), bins.size-1))
    for non_pair_idx, non_pair in enumerate(random_pairs.values()):
        prop_transmitted_non_pairs[non_pair_idx,:] = \
            calc_prop_transmitted_binned(non_pair, bins, min_af=min_af)
    p_vals = {}
    for bin_idx in [1,2,3]:
        bin_p_value = scipy.stats.kstest(prop_transmitted_pairs[:, bin_idx],
            prop_transmitted_non_pairs[:, bin_idx])
        p_vals[bins[bin_idx]] = bin_p_value
        print(f'pvalue comparing prop of shared variants between pairs and non pairs for the bin \
with the lower limit {bins[bin_idx]} is {bin_p_value[1]}')
    return(p_vals)


def nb_change_t_test(nb_dat=None, realign_vcf_dat=None, max_donor_af_cutoff=0.06):
    to_test = []
    for pair in nb_dat['pair']:
        donor_dat = realign_vcf_dat[pair[0] + '_' + pair[0]]
        # only alleles that passed filtering
        donor_dat = donor_dat[donor_dat['FILTER_PASS'] == True]
        if donor_dat['AF'].max() > max_donor_af_cutoff:
            to_test.append(pair)
    test_dat = nb_dat[nb_dat['pair'].isin(to_test)]
    t_test = scipy.stats.ttest_rel(test_dat['1%BottleneckSize_mle'], 
        test_dat['3%BottleneckSize_mle'])
    print(f'pvalue comparing 1%BottleneckSize_mle and 3%BottleneckSize_mle for all samples with max donor AF > {max_donor_af_cutoff} is \
    {t_test[1]}')
    return(t_test)


def prop_vars_shared(realign_vcf_dat=None, transmission_pairs=None, nb_dict=None, donor_af_range=[0, 1.0], min_recipient_af=0.0, min_nb=0.0):
    n_donor_vars = 0
    n_shared_vars = 0
    for pair_id, pair in transmission_pairs.items():
        if nb_dict[pair] >= min_nb:
            pair_dat = \
                get_pair_dat(realign_vcf_dat[pair[0] + '_' + pair[0]], 
                    realign_vcf_dat[pair[1] + '_' + pair[0]])
            # gets only variants which pass filter in both samples
            # and are present in at least one of the samples 
            # at at least the min recipient frequency
            pair_dat = \
                filter_pair_dat(pair_dat, 
                    min_af=min(min_recipient_af, donor_af_range[0]), 
                    only_shared=False)
            # only ones present in the the donor in the donor_af_range
            pair_dat = \
                pair_dat[(pair_dat['AF_x'] <= donor_af_range[1]) & (pair_dat['AF_x'] >= donor_af_range[0])]
            n_donor_vars += pair_dat.shape[0]
            n_shared_vars += pair_dat[(pair_dat['AF_y'] >= min_recipient_af)].shape[0]
    print(f'of the {n_donor_vars} variants present in donor samples at >= {donor_af_range[0]} and <= {donor_af_range[1]}, \
only {n_shared_vars} ({n_shared_vars/n_donor_vars}) are present in recipients at >= {min_recipient_af}')
    return(n_donor_vars, n_shared_vars, n_donor_vars/n_shared_vars)



def get_consensus_diffs(pair_dat, threshold=0.5):
    return(pair_dat[((pair_dat['AF_x'] > threshold) & 
            (pair_dat['AF_y'] < threshold)) | 
        ((pair_dat['AF_x'] < threshold) & 
            (pair_dat['AF_y'] > threshold))])
    


def get_pair_dat(dat_1, dat_2):
    pair_dat = dat_1.merge(dat_2, 
        left_on=['POS', 'REF', 'ALT'], 
        right_on=['POS', 'REF', 'ALT'], 
        how='outer')
    pair_dat.loc[:,['AF_x', 'AF_y']] = \
        pair_dat.loc[:,['AF_x', 'AF_y']].fillna(0)
    pair_dat.loc[:,['AF_x', 'AF_y']] = \
        pair_dat.loc[:,['AF_x', 'AF_y']].astype(float)
    return(pair_dat)


def get_snps(pairs, vcf_dict):
    out = pd.DataFrame(pairs).T
    out.columns = ['donor', 'recipient']
    pair_dat = \
        [get_pair_dat(vcf_dict[i[0]+'_'+i[0]], vcf_dict[i[1]+'_'+i[0]]) for i in pairs.values()]
    diff_dat = [get_consensus_diffs(i) for i in pair_dat]
    n_diffs, diffs = [[i.shape[0] for i in diff_dat], [list(i['REF']+i['POS'].astype(str)+i['ALT']) for i in diff_dat]]
    out['n_diffs'] = n_diffs
    out['diffs'] = diffs
    print(out.to_string())
    return(out)    


def hpd(dat, qs):
    width = qs[2] - qs[0]
    # sorts from smallest to largest
    dat = sorted(dat)
    # length of data
    n = len(dat)
    # number of values we are keeping
    # thus, the index of the begining of 
    # the HPD interval must be <= this value
    # this gives us the tail of the distribution
    interval_idx_inc = int(np.floor(width * n))
    # number of values we are excluding
    # thus, possible number of HPD intervals
    # this gives us the head of the distribution
    n_intervals = n - interval_idx_inc
    # the width of each possible interval
    # for each possible head and tail value, 
    # what is the difference between them
    interval_width = [a_i - b_i for a_i, b_i in 
                      zip(dat[interval_idx_inc:], 
                          dat[:n_intervals])]
    # find the shortest interval
    min_idx = interval_width.index(min(interval_width))
    hpd_interval = (dat[min_idx], dat[min_idx+interval_idx_inc])
    dat_hpd = [item for item in dat if (item >= hpd_interval[0]) & (item <= hpd_interval[1])]
    dat_mid = np.quantile(dat_hpd, qs[1])
    return((hpd_interval[0], dat_mid, hpd_interval[1]))


def calc_corr(dat, reps=1000, alpha=0.05):
    from scipy.stats import spearmanr
    rng = np.random.default_rng(10)
    corrs = []
    for rep in range(reps):
        sampled_dat = dat.iloc[rng.integers(dat.shape[0], size=dat.shape[0]),:]
        corrs.append(spearmanr(sampled_dat.iloc[:,0], sampled_dat.iloc[:,1])[0])
    corr_ci = hpd(corrs, [alpha/2, 0.5, 1-alpha/2])
    return(corr_ci)


def calc_shared_corr(realign_vcf_dat, transmission_pairs, nb_dict, min_nb=0, min_af=0.01, simulated_dat)
    from numpy import cov
    from scipy.stats import spearmanr
    # for simulated data get shared variants
    simulated_shared = simulated_data[(simulated_data[0] > min_af) & (simulated_data[1] > min_af)]
    # now get the correlation between them
    simulated_corr_ci = calc_corr(simulated_shared)
    # now do the observed data
    all_pair_dat = pd.DataFrame()
    for pair_id, pair in transmission_pairs.items():
        if nb_dict[pair] >= min_nb:
            pair_dat = \
                get_pair_dat(realign_vcf_dat[pair[0] + '_' + pair[0]], 
                    realign_vcf_dat[pair[1] + '_' + pair[0]])
            pair_dat = \
                filter_pair_dat(pair_dat, 
                    min_af=min_af, only_shared=True)
            all_pair_dat = pd.concat([all_pair_dat, pair_dat[['AF_x', 'AF_y']]])
    observed_corr_ci = calc_corr(all_pair_dat)


def run():
parser = argparse.ArgumentParser()
# input files
# todo outdir currently doesn't get used by anything
parser.add_argument('--metadata', 
    default='data/abe255_Data_file_format.csv')
parser.add_argument('--vcfDir', 
    default='data/seq/*/*_filter_norm.vcf',
    help='directory with all vcf files')
parser.add_argument('--vcfRealignDir', 
    default='data/seq_realign/*/*_lofreq_filter_norm.vcf',
    help='directory with realigned vcf files')
parser.add_argument('--nbData', 
    default='data/Popa_new_Nb_estimates.csv',
    help='file with reanalyzed nb data')
parser.add_argument('--filterAllow',
    default=['PASS', 'min_af_0.010000'],
    nargs='+',
    help='which vcf filter strings to allow')
parser.add_argument('--simulatedData', 
    default='data/stochastic_sims_6percent.csv',
    help='file with simulated nb data')
args = parser.parse_args()
# reads in metadata
metadata = pd.read_csv(args.metadata, sep=',')
transmission_pairs = get_transmission_pairs(metadata)
# reads transmission pair VCF data (variants relative to wuhan) into a dictionary
# process results
vcf_dat = {}
for path in glob.glob(args.vcfDir):
    sample, sample_dat = \
        process_vcf(path, args.filterAllow)
    vcf_dat[sample] = sample_dat

# reads transmission pair VCF data (variants relative to donor) into a dictionary
# process results
realign_vcf_dat = {}
for path in glob.glob(args.vcfRealignDir):
    sample_ref=  '_'.join(['_'.join(i.split('_')[0:2]) for i in path.split('/')[-2].split('-')])
    sample, sample_dat = \
        process_vcf(path, args.filterAllow)
    realign_vcf_dat[sample_ref] = sample_dat


    # for each donor/recipient pair, get their SNP difference
    dist_dat = get_snps(transmission_pairs, realign_vcf_dat)
    dist_dat.to_csv('.'.join(args.metadata.split('.')[:-1]) + '_dists.tsv', sep='\t', header=None, index=None)
    #### NB STATS ####
    nb_dat = pd.read_csv(args.nbData)
    nb_dat['donor'] = \
        'CoV_' + nb_dat['donor'].astype(str).str.zfill(3)
    nb_dat['recipient'] = \
        'CoV_' + nb_dat['recipient'].astype(str).str.zfill(3)
    nb_dat['pair'] = list(zip(nb_dat['donor'], nb_dat['recipient']))
    nb_dict = {i['pair']: i['1%BottleneckSize_mle'] for idx, i in nb_dat.iterrows()}
    # t_test on botleneck estimate change between 1% and 3%
    t_test_out = nb_change_t_test(nb_dat=nb_dat, realign_vcf_dat=realign_vcf_dat)
    # KS test on whether known transmission pairs share more variants
    if not os.path.exists('data/random_pairs.csv'):
        raise Exception('random pairs have not been generated, run plot_figures.py to generate random pairs')
    random_pairs = pd.read_csv('data/random_pairs.csv', header=None)
    random_pairs = {idx: (i[0], i[1]) for idx, i in random_pairs.iterrows()}
    prop_shared_test_out = prop_shared_test(vcf_dat=vcf_dat, 
                transmission_pairs=transmission_pairs, 
                random_pairs=random_pairs, min_af=0.01)
    # calculates proportion of variants within the donor_af_range
    # that are present at > min_recipient_af
    prop_shared_out = \
        prop_vars_shared(realign_vcf_dat=realign_vcf_dat, transmission_pairs=transmission_pairs, 
            nb_dict=nb_dict, donor_af_range=[0.02, 0.06], min_recipient_af=0.01, min_nb=1000)
    # calcualtes correlations between shared variant frequencies
    shared_var_corr = calc_shared_corr(realign_vcf_dat, pd.read_csv(args.simulatedData, sep=',', header=None))

if __name__ == "__main__":
    run()

