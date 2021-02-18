import argparse
import pandas as pd
import numpy as np
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import string
import itertools
import os


def plot_style(grey='#333333'):
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = 'Helvetica'
    mpl.rcParams['font.weight'] = 'light'
    mpl.rcParams['text.color'] = grey
    mpl.rcParams['axes.labelcolor'] = grey
    mpl.rcParams['xtick.color'] = grey
    mpl.rcParams['ytick.color'] = grey
    # Font sizes
    mpl.rcParams['figure.titlesize'] = 16
    mpl.rcParams['axes.titlesize'] = 16
    mpl.rcParams['axes.labelsize'] = 16
    mpl.rcParams['xtick.labelsize'] = 16
    mpl.rcParams['ytick.labelsize'] = 16
    # Border colors
    mpl.rcParams['axes.edgecolor'] = grey


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

     
def get_clonal_diffs(pair_dat, threshold=0.05):
    return(pair_dat[((pair_dat['AF_x'] < threshold) & 
            (pair_dat['AF_y'] > 1-threshold)) | 
        ((pair_dat['AF_x'] > 1-threshold) & 
            (pair_dat['AF_y'] < threshold))])


def plot_nb_estimates(ax, nb_dat=None, realign_vcf_dat=None):
    max_donor_af = {}
    for pair in nb_dat['pair']:
        donor_dat = realign_vcf_dat[pair[0] + '_' + pair[0]]
        # only alleles that passed filtering
        donor_dat = donor_dat[donor_dat['FILTER_PASS'] == True]
        # max allele frequency
        max_donor_af[(pair[0], pair[1])] = donor_dat['AF'].max()
    nb_dat['max_donor_af'] = nb_dat['pair'].map(max_donor_af)
    jitter = [-0.005, 0, 0.005]
    cols = ['#333333', '#D08770', '#5E81AC']
    for row_idx,row in nb_dat.iterrows():
        for nb_idx,nb in enumerate(['1%BottleneckSize', '3%BottleneckSize', '6%BottleneckSize']):
            # CI
            _ = ax.plot([row['max_donor_af']-jitter[nb_idx]]*2,
                [row[f'{nb}_lower95'], row[f'{nb}_upper95']], 
                lw=1, color=cols[nb_idx], zorder=2)
            # MLE
            _ = ax.scatter(row['max_donor_af']-jitter[nb_idx],
                row[f'{nb}_mle'],
                color=cols[nb_idx],
                edgecolor='#333333', 
                alpha=0.85, zorder=3)
    ax.text(0.75,0.925,
            '1% cutoff', 
            color='#333333', size=15,
            transform=ax.transAxes)
    ax.text(0.75,0.875,
            '3% cutoff', 
            color='#D08770', size=15,
            transform=ax.transAxes)
    ax.text(0.75,0.825,
            '6% cutoff', 
            color='#5E81AC', size=15,
            transform=ax.transAxes)
    ax.set_xlabel('max. donor AF')
    ax.set_ylabel('bottleneck estimate')
    ax.set_yscale('log')
    ax.set_aspect(1./ax.get_data_ratio())
    return(ax)


def plot_tv(ax, tv_vcf_dat=None, pair=None, ct_dat=None, lod=0.01, min_af=0.0):  
    tv_plot_dat = \
        get_pair_dat(tv_vcf_dat[pair[0] + '_' + 
                pair[0]], 
            tv_vcf_dat[pair[1] + '_' +
                pair[0]])
    # we want all variants not just those that are shared
    # remove any < a certain threshold
    tv_plot_dat = \
        filter_pair_dat(tv_plot_dat, only_shared=False, 
            min_af=min_af)
    donor_name = \
                f'{pair[0]} (CT value: {ct_dat[pair[0]]})'
    recipient_name = \
       f'{pair[1]} (CT value: {ct_dat[pair[1]]})'
    axins = inset_axes(ax, height="60%", width="50%", borderpad=3, loc=4)
    ax.scatter(tv_plot_dat['AF_x'], 
        tv_plot_dat['AF_y'],
        alpha=0.65,
        color='#333333', edgecolor='#333333', zorder=3)
    axins.scatter(tv_plot_dat['AF_x'], 
        tv_plot_dat['AF_y'],
        alpha=0.65,
            color='#333333', edgecolor='#333333', zorder=3)
    for i_ax in [ax, axins]:
        i_ax.grid(True,  ls='--')
        i_ax.set_axisbelow(True)
        i_ax.axvline(x=lod, ls='--', color='#BF616A', zorder=1)
        i_ax.axhline(y=lod, ls='--', color='#BF616A', zorder=1)
        i_ax.axvline(x=1-lod, ls='--', color='#BF616A', zorder=1)
        i_ax.axhline(y=1-lod, ls='--', color='#BF616A', zorder=1)
        i_ax.set_aspect('equal')
    ax.set_xlabel(donor_name)
    ax.set_ylabel(recipient_name)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(-0.05, 1.05)
    big_ticks = [0,0.25,0.5,0.75,1.0]
    ax.set_xticks(big_ticks)
    ax.set_yticks(big_ticks)
    axins.set_ylim(-0.0025, 0.0625)
    axins.set_xlim(-0.0025, 0.0625)
    axins.tick_params(axis="both", labelsize=12)
    small_ticks=[0.00, 0.02, 0.04, 0.06]
    axins.set_xticks(small_ticks)
    axins.set_yticks(small_ticks)
    return(ax)


def plot_snv_density_indiv(ax, density_vcf_dat=None, transmission_pairs=None, kde_range=[0.00, 1.0], title=None):
    from scipy import stats
    # iterate through transmission pairs
    # get pair dat, making share variant didn't fail filtering in either sample
    # get only shared variants
    # determine if max AF in recipient is > threshold and append 
    # shared variant freqs in donor to appropriate list
    kde_range.sort()        # operates in place
    x_s = np.linspace(0.0,kde_range[1], 1000)
    #cols = ['#81a1c1', '#a3be8c', '#A78CBE']
    cols = ['#5E81AC', '#D08770','#333333']
    for pair_id, pair in transmission_pairs.items():
        pair_dat = \
            get_pair_dat(density_vcf_dat[pair[0] + '_' + pair[0]], 
                density_vcf_dat[pair[1] + '_' + pair[0]])
        filtered_pair_dat = \
            filter_pair_dat(pair_dat, min_af=kde_range[0], only_shared=True)
        kde_dat = filtered_pair_dat['AF_x']
        kde_dat = kde_dat[(kde_dat >= kde_range[0]) & (kde_dat <= kde_range[1])]
        pair_kde = stats.gaussian_kde(kde_dat)
        # we need to find the frequency of the maximum clonal variant in the recipient
        # all pair variants, not just shared
        filtered_all_pair_dat = \
            filter_pair_dat(pair_dat, min_af=0.00, only_shared=False)
        only_recipient = \
            filtered_all_pair_dat[(filtered_all_pair_dat['AF_x'] < 0.0001)]
        max_clonal_af = only_recipient['AF_y'].max()
        if max_clonal_af > 1-kde_range[1]:
            ax.plot(x_s, pair_kde(x_s), color=cols[2], alpha=1.0, zorder=3)
        elif max_clonal_af > kde_range[1]:
            ax.plot(x_s, pair_kde(x_s), color=cols[1], alpha=0.5, zorder=2)
        else:
            ax.plot(x_s, pair_kde(x_s), color=cols[0], alpha=0.5, zorder=2)
    t0 = ax.text(0.57,0.925,
        'max. de novo AF', 
        color='#333333', size=15,
        transform=ax.transAxes)
    t1 = ax.text(0.57,0.875,
        '$>$ ' + str(1-kde_range[1]), 
        color=cols[2], size=15,
        transform=ax.transAxes)
    t2 = ax.text(0.57,0.825,
        '$>$ ' + str(kde_range[1]) + ', $\leq$ ' + str(1-kde_range[1]), 
        color=cols[1], size=15,
        transform=ax.transAxes)
    t3 = ax.text(0.57,0.775,
        '$\leq$ ' + str(kde_range[1]), 
        color=cols[0], size=15,
        transform=ax.transAxes)
    _ = \
        [t.set_bbox(dict(facecolor='white', 
            alpha=0.5, pad=0, edgecolor='none')) for t in [t0, t1,t2,t3]]
    #ax.set_xticks([0.00, 0.01, 0.02, 0.03, 0.04, 0.05])
    ax.set_xlabel('shared AF in donor')
    ax.set_xlim(-0.005, kde_range[1] + 0.005)
    ax.set_xticks([0,0.02, 0.04, 0.06]) # hard coded
    ax.set_ylabel('density')
    ax.set_yticks([])
    ax.set_aspect(1./ax.get_data_ratio())
    return(ax)


def plot_multi_tv(ax, tv_vcf_dat=None, transmission_pairs=None, nb_dict=None, simulated_data=None, lod=0.0, min_af=0.0, min_nb=0.0):
    for pair_id, pair in transmission_pairs.items():
        if nb_dict[pair] >= min_nb:
            pair_dat = \
                get_pair_dat(tv_vcf_dat[pair[0] + '_' + pair[0]], 
                    tv_vcf_dat[pair[1] + '_' + pair[0]])
            pair_dat = \
                filter_pair_dat(pair_dat, 
                    min_af=min_af, only_shared=False)
            ax.scatter(pair_dat['AF_x'], 
                pair_dat['AF_y'],
                alpha=0.65,
                color='#333333', edgecolor='#333333', zorder=3)
    #axins = inset_axes(ax, height="50%", width="40%", borderpad=3, loc=4)
    ax.scatter(simulated_data[0],
        simulated_data[1], 
        alpha=0.65, 
        color='#D08770', edgecolor='#333333', zorder=2)
    ax.set_xlabel('donor AF')
    ax.set_ylabel('recipient AF')
    ax.grid(True,  ls='--')
    ax.set_axisbelow(True)
    ax.axvline(x=lod, ls='--', color='#BF616A', zorder=1)
    ax.axhline(y=lod, ls='--', color='#BF616A', zorder=1)
    ax.axvline(x=1-lod, ls='--', color='#BF616A', zorder=1)
    ax.axhline(y=1-lod, ls='--', color='#BF616A', zorder=1)
    ax.set_aspect('equal')
    ax.set_ylim(-0.0025, 0.0625)
    ax.set_xlim(-0.0025, 0.0625)
    t1 = ax.text(0.075,0.925,
            f'observed ($N_b \geq$ {min_nb})', 
            color='#333333', size=15,
            transform=ax.transAxes)
    t2 = ax.text(0.075,0.875,
            'simulated ($N_b$ = 1000)', 
            color='#D08770', size=15,
            transform=ax.transAxes)
    _ = \
        [t.set_bbox(dict(facecolor='white', 
            alpha=0.5, pad=0, edgecolor='none')) for t in [t1,t2]]
    small_ticks=[0.00, 0.02, 0.04, 0.06]
    ax.set_xticks(small_ticks)
    ax.set_yticks(small_ticks)
    return(ax)


def plot_tv_all(vcf_dat=None, transmission_pairs=None, ct_dat=None, lod=0.01, min_af=0.0):
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages('figures/figure_s1.pdf') as pdf:
        for pair_idx, pair in transmission_pairs.items():
            tv_plot_dat = \
                get_pair_dat(vcf_dat[pair[0] + '_' + 
                        pair[0]], 
                    vcf_dat[pair[1] + '_' +
                        pair[0]])
            # we want all variants not just those that are shared
            # remove any < a certain threshold
            tv_plot_dat = \
                filter_pair_dat(tv_plot_dat, only_shared=False, 
                    min_af=min_af)
            donor_name = \
                f'{pair[0]} (CT value: {ct_dat[pair[0]]})'
            recipient_name = \
               f'{pair[1]} (CT value: {ct_dat[pair[1]]})'
            fig, axs = plt.subplots(1,2, figsize=(10, 4.8), 
                constrained_layout=True)
            for ax in axs:
                ax.scatter(tv_plot_dat['AF_x'], 
                    tv_plot_dat['AF_y'],
                    alpha=0.65,
                    color='#333333', edgecolor='#333333', zorder=3)
                ax.grid(True,  ls='--')
                ax.set_axisbelow(True)
                ax.axvline(x=lod, ls='--', color='#BF616A', zorder=1)
                ax.axhline(y=lod, ls='--', color='#BF616A', zorder=1)
                ax.axvline(x=1-lod, ls='--', color='#BF616A', zorder=1)
                ax.axhline(y=1-lod, ls='--', color='#BF616A', zorder=1)
                ax.set_aspect('equal')
                ax.set_xlabel(donor_name)
                ax.set_ylabel(recipient_name)
            axs[0].set_ylim(-0.05, 1.05)
            axs[0].set_xlim(-0.05, 1.05)
            big_ticks = [0,0.25,0.5,0.75,1.0]
            axs[0].set_xticks(big_ticks)
            axs[0].set_yticks(big_ticks)
            axs[1].set_ylim(-0.0025, 0.0625)
            axs[1].set_xlim(-0.0025, 0.0625)
            small_ticks=[0.00, 0.02, 0.04, 0.06]
            axs[1].set_xticks(small_ticks)
            axs[1].set_yticks(small_ticks)
            label = string.ascii_lowercase[pair_idx//26] + \
                string.ascii_lowercase[pair_idx - pair_idx//26 * 26]
            axs[0].text(-0.25, 1.05, label,
                size=20, weight='bold',  transform=axs[0].transAxes)
            pdf.savefig(fig)
            plt.close()


def generate_random_pairs(transmission_pairs=None, family_dat=None):
    success=False
    # doing this in a while loop to make sure we successfuly
    # get a full set of random pairings
    while success == False:
        all_recipients = \
            set([val[1] for val in transmission_pairs.values()])
        non_pairs = []
        for pair_idx, pair in enumerate(transmission_pairs.values()):
            # get all recipients that this donor transmitted too
            donor_recipients = \
                set([val[1] for val in transmission_pairs.values() if val[0] == pair[0]])
            # removes recipients which have already been sampled from this donor
            # removes recipients which have already been sampled
            random_recipients = all_recipients - \
                set(pair) - \
                donor_recipients
            # removes any random recipients in the same family as the donor or recipient
            random_recipients = \
                [i for i in random_recipients if family_dat[i] != family_dat[pair[0]] and 
                    family_dat[i] != family_dat[pair[1]]]
            if len(random_recipients) == 0:
                break
            if pair[0] in random_recipients:
                raise Exception('donor in random recipients')
            if pair[1] in random_recipients:
                raise Exception('real recipient in random recipients')
            random_recipient = np.random.choice(random_recipients, size=1)[0]
            non_pairs.append([pair[0], random_recipient])
            all_recipients = all_recipients - set([random_recipient])
            print(all_recipients)
        success = True
        if len(non_pairs) != 39:
            raise Exception('non pairs list is either too long or too short')
    pd.DataFrame(non_pairs).to_csv('data/random_pairs.csv', header=None, index=None)


def plot_prob_transmission(ax, vcf_dat=None, transmission_pairs=None, random_pairs=None, family_dat=None, min_af=0.0):
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
    pairs_mean = np.mean(prop_transmitted_pairs, axis=0)
    non_pairs_mean = np.mean(prop_transmitted_non_pairs, axis=0)
    x_vals = np.arange(1,4)
    x_val_expand = 1.5
    for x_idx, x_val in enumerate(x_vals):
        _ = ax.plot([x_val*x_val_expand-0.2-0.125, x_val*1.5-0.2+0.125],
            [pairs_mean[x_val], pairs_mean[x_val]],
            color='#5E81AC', lw=3, zorder=2)
        _ = ax.plot([x_val*x_val_expand+0.2-0.125, x_val*1.5+0.2+0.125],
            [non_pairs_mean[x_val], non_pairs_mean[x_val]],
            color='#333333', lw=3, zorder=2)
    for row_idx, row in enumerate(prop_transmitted_pairs):
        _ = ax.scatter(x_vals*x_val_expand - 0.2 + np.random.normal(scale=0.05), 
            row[x_vals], 
            color='#5E81AC', edgecolor='#333333',
            alpha=0.65, zorder=3)
    for row_idx, row in enumerate(prop_transmitted_non_pairs):
        _ = ax.scatter(x_vals*x_val_expand + 0.2 + np.random.normal(scale=0.05), 
            row[x_vals], 
            color='#333333', edgecolor='#333333',
            alpha=0.65, zorder=3)
    _ = ax.set_xticks(x_vals*x_val_expand)
    _ = ax.set_xticklabels(['1-2%', '2-3%', '3-6%'])
    t1 = ax.text(0.45,0.875,
                'transmission pairs', 
                color='#5E81AC', size=15,
                transform=ax.transAxes)
    t2 = ax.text(0.45,0.825,
            'non transmission pairs', 
            color='#333333', size=15,
            transform=ax.transAxes)
    _ = \
        [t.set_bbox(dict(facecolor='white', 
            alpha=0.5, pad=0, edgecolor='none')) for t in [t1,t2]]
    ax.set_xlabel('donor AF')
    ax.set_ylabel('prob. of transmission')
    ax.set_aspect(1./ax.get_data_ratio())
    return(ax)


def plot_shared_vars(ax, vcf_dat=None, transmission_pairs=None, n_plot=5, min_af=0.01):
    # get all samples involved in transmission pairs
    transmission_samples = \
        list(set(itertools.chain(*transmission_pairs.values())))
    subclonal_dat = pd.DataFrame()
    for sample in transmission_samples:
        sample_dat = vcf_dat[sample]
        sample_dat = sample_dat[sample_dat['FILTER_PASS'] == True]
        sample_dat = \
            sample_dat[(sample_dat['AF'] <= 1-min_af) & 
                (sample_dat['AF'] >= min_af)]
        sample_dat['label'] = \
            sample_dat['REF'] + \
            sample_dat['POS'].astype(str) + \
            sample_dat['ALT']
        subclonal_dat = \
            pd.concat([subclonal_dat, sample_dat])
    subclonal_dat_counted = \
        subclonal_dat.groupby('label').size().\
            sort_values(ascending=False).reset_index()
    to_plot = \
        subclonal_dat[subclonal_dat['label'].isin(subclonal_dat_counted['label'][0:n_plot])]
    x_val_dict = \
        {idx: i['label'] for idx, i in subclonal_dat_counted.iterrows()}
    rev_x_val_dict = \
        {i['label']: idx for idx, i in subclonal_dat_counted.iterrows()}
    count_dict = \
        {idx: i[0] for idx, i in subclonal_dat_counted.iterrows()}
    for idx, row in to_plot.iterrows():
        x_val = rev_x_val_dict[row['label']] + np.random.normal(scale=0.075)
        _ = ax.scatter(x_val, row['AF'], color='#333333', alpha=0.65)
    ax.set_xticks(subclonal_dat_counted.index[0:n_plot])
    x_ticklabels = ax.get_xticks()
    ax.set_xticklabels([f'{x_val_dict[i]}\nn={count_dict[i]}' for i in x_ticklabels], size=7, rotation=45)
    ax.axhline(min_af, ls='--', color='#BF616A')
    ax.set_xlabel('variant')
    ax.set_ylabel('frequency')
    ax.set_ylim(0,)
    ax.set_aspect(1./ax.get_data_ratio())
    return(ax)


def plot_shared_vars_all(vcf_dat=None, transmission_pairs=None, min_af=0.0, min_abundance=3):
    from collections import Counter
    # get all samples involved in transmission pairs
    transmission_samples = \
        list(set(itertools.chain(*transmission_pairs.values())))
    subclonal_dat = pd.DataFrame()
    for sample in transmission_samples:
        sample_dat = vcf_dat[sample]
        sample_dat = sample_dat[sample_dat['FILTER_PASS'] == True]
        sample_dat = \
            sample_dat[(sample_dat['AF'] <= 1-min_af) & 
                (sample_dat['AF'] >= min_af)]
        sample_dat['label'] = \
            sample_dat['REF'] + \
            sample_dat['POS'].astype(str) + \
            sample_dat['ALT']
        subclonal_dat = \
            pd.concat([subclonal_dat, sample_dat])
    # want to filter so that we only keep rows 
    # where the variant is in at least min_abundance samples
    subclonal_dat_counted = \
        subclonal_dat.groupby('label').size().\
            sort_values(ascending=False).reset_index()
    vars_to_plot = \
        subclonal_dat_counted[subclonal_dat_counted[0] >= min_abundance]
    to_plot = \
        subclonal_dat[subclonal_dat['label'].isin(vars_to_plot['label'].unique())]
    x_val_dict = \
        {idx: i['label'] for idx, i in subclonal_dat_counted.iterrows()}
    rev_x_val_dict = \
        {i['label']: idx for idx, i in subclonal_dat_counted.iterrows()}
    count_dict = \
        {idx: i[0] for idx, i in subclonal_dat_counted.iterrows()}
    fig, ax = plt.subplots(figsize=(50, 6), constrained_layout=True)
    for idx, row in to_plot.iterrows():
        x_val = rev_x_val_dict[row['label']] + np.random.normal(scale=0.075)
        _ = ax.scatter(x_val, row['AF'], color='#333333', edgecolor='#333333', alpha=0.65)
    ax.set_xticks(vars_to_plot.index)
    x_ticklabels = ax.get_xticks()
    ax.set_xticklabels([f'{x_val_dict[i]}\nn={count_dict[i]}' for i in x_ticklabels], size=4, rotation=90)
    ax.set_xlim(-1, vars_to_plot.index[-1]+1)
    ax.axhline(min_af, ls='--', color='#BF616A')
    ax.set_xlabel('variant', size=24)
    ax.set_ylabel('frequency', size=24)
    ax.set_yscale('log')
    fig.savefig('figures/figure_s2.pdf')  
    plt.close()  


def plot_var_abundance(vcf_dat=None, transmission_pairs=None, min_af=0.0, min_abundance=3):
    transmission_samples = \
        list(set(itertools.chain(*transmission_pairs.values())))
    subclonal_dat = pd.DataFrame()
    for sample in transmission_samples:
        sample_dat = vcf_dat[sample]
        sample_dat = sample_dat[sample_dat['FILTER_PASS'] == True]
        sample_dat = \
            sample_dat[(sample_dat['AF'] <= 1-min_af) & 
                (sample_dat['AF'] >= min_af)]
        sample_dat['label'] = \
            sample_dat['REF'] + \
            sample_dat['POS'].astype(str) + \
            sample_dat['ALT']
        subclonal_dat = \
            pd.concat([subclonal_dat, sample_dat])
    # want to filter so that we only keep rows 
    # where the variant is in at least min_abundance samples
    subclonal_dat_counted = \
        subclonal_dat.groupby('label').size().\
            sort_values(ascending=False).reset_index()
    vars_to_plot = \
        subclonal_dat_counted[subclonal_dat_counted[0] >= min_abundance]
    fig, ax = plt.subplots(figsize=(6.4*2, 4.8), constrained_layout=True)
    ax.scatter(vars_to_plot.index, 
        vars_to_plot[0], 
        color='#333333',
        edgecolor='#333333',
        alpha=0.65)
    ax.set_xlabel('variant')
    ax.set_xticks([])
    ax.set_xlim(-4, vars_to_plot.index[-1]+4)
    ax.set_ylim(0,)
    ax.set_ylabel('abundance')
    fig.savefig('figures/figure_s3.pdf')



def figure_1(nb_dat=None, realign_vcf_dat=None, vcf_dat=None, transmission_pairs=None, 
    family_dat=None, ct_dat=None, simulated_data=None, highlight_pair=None):
    fig, axs = \
        plt.subplots(2, 3, figsize=(15, 10), 
            constrained_layout=True)
    axs[0,0] = \
        plot_nb_estimates(axs[0,0], nb_dat=nb_dat, realign_vcf_dat=realign_vcf_dat)
    axs[0,1] = \
        plot_tv(axs[0,1], tv_vcf_dat=realign_vcf_dat, pair=highlight_pair, 
            ct_dat=ct_dat, lod=0.01)
    axs[0,2] = \
        plot_snv_density_indiv(axs[0,2], density_vcf_dat=realign_vcf_dat, transmission_pairs=transmission_pairs, 
        kde_range=[0.01, 0.06])
    nb_dict = {i['pair']: i['1%BottleneckSize_mle'] for idx, i in nb_dat.iterrows()}
    axs[1,0] = \
        plot_multi_tv(axs[1,0], tv_vcf_dat=realign_vcf_dat, transmission_pairs=transmission_pairs, nb_dict=nb_dict,
            simulated_data=simulated_data, lod=0.01, min_nb=1000)
    # if we haven't generated the random pairs yet, do that
    if not os.path.exists('data/random_pairs.csv'):
        generate_random_pairs(transmission_pairs=transmission_pairs, family_dat=family_dat)
    random_pairs = pd.read_csv('data/random_pairs.csv', header=None)
    random_pairs = {idx: (i[0], i[1]) for idx, i in random_pairs.iterrows()}
    axs[1,1] = \
        plot_prob_transmission(axs[1,1], vcf_dat=vcf_dat, 
            transmission_pairs=transmission_pairs, 
            random_pairs=random_pairs, 
            family_dat=family_dat, min_af=0.01)
    axs[1,2] = \
        plot_shared_vars(axs[1,2], vcf_dat=vcf_dat, transmission_pairs=transmission_pairs, n_plot=12, min_af=0.01)
    for ax_idx, ax in enumerate(axs.flat):
        ax.text(-0.175, 1.0, string.ascii_lowercase[ax_idx], transform=ax.transAxes, 
                size=20, weight='bold')
    fig.savefig('figures/figure_1.pdf')
    plt.close()


def figure_2(nb_dat=None):
    cols = ['#333333', '#D08770', '#5E81AC']
    jitter = [-0.25, 0, 0.25]
    fig, ax = plt.subplots(figsize=(6.4*2, 4.8), constrained_layout=True)
    for row_idx, row in nb_dat.iterrows():
        for nb_idx,nb in enumerate(['1%BottleneckSize', '3%BottleneckSize', '6%BottleneckSize']):
            # CI
            _ = ax.plot([row_idx-jitter[nb_idx]]*2,
                [row[f'{nb}_lower95'], row[f'{nb}_upper95']], 
                lw=1, color=cols[nb_idx], zorder=2)
            # MLE
            _ = ax.scatter(row_idx-jitter[nb_idx],
                row[f'{nb}_mle'],
                color=cols[nb_idx],
                edgecolor='#333333', 
                alpha=0.85, zorder=3)
    ax.set_xticks(nb_dat.index)
    x_ticklabels = ax.get_xticks()
    ax.set_xticklabels([f'{nb_dat.iloc[i,:]["donor"]} -> {nb_dat.iloc[i,:]["recipient"]}' \
        for i in x_ticklabels], rotation=45, size=8, ha='right')
    t1 = ax.text(0.91,0.925,
                '1% cutoff', 
                color='#333333', size=15,
                transform=ax.transAxes)
    t2 = ax.text(0.91,0.86,
            '3% cutoff', 
            color='#D08770', size=15,
            transform=ax.transAxes)
    t3 = ax.text(0.91,0.795,
            '6% cutoff', 
            color='#5E81AC', size=15,
            transform=ax.transAxes)
    _ = \
        [t.set_bbox(dict(facecolor='white', 
            alpha=0.5, pad=0, edgecolor='none')) for t in [t1,t2, t3]]
    ax.set_xlabel('transmission pair')
    ax.set_ylabel('bottleneck esitmate')
    fig.savefig('figures/figure2.pdf')


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
        default='data/seq_realign/*/*_filter_norm.vcf',
        help='directory with realigned vcf files')
    parser.add_argument('--simulatedData', 
        default='data/stochastic_sims_6percent.csv',
        help='file with simulated nb data')
    parser.add_argument('--nbData', 
        default='data/Popa_new_Nb_estimates.csv',
        help='file with reanalyzed nb data')
    parser.add_argument('--highlightPair',
        default=['CoV_162', 'CoV_161'],
        nargs=2, 
        help='transmission pair to highlight')
    parser.add_argument('--filterAllow',
        default=['PASS', 'min_af_0.010000'],
        nargs='+',
        help='which vcf filter strings to allow')
    args = parser.parse_args()
    plot_style()
    # reads in metadata
    metadata = pd.read_csv(args.metadata, sep=',')
    transmission_pairs = get_transmission_pairs(metadata)
    family_dat = {i['sample_name']: i['family'] for idx, i in metadata.iterrows()}
    ct_dat = {i['sample_name']: i['ct'] for idx, i in metadata.iterrows()}
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
    nb_dat = pd.read_csv(args.nbData)
    nb_dat['donor'] = \
        'CoV_' + nb_dat['donor'].astype(str).str.zfill(3)
    nb_dat['recipient'] = \
        'CoV_' + nb_dat['recipient'].astype(str).str.zfill(3)
    nb_dat['pair'] = list(zip(nb_dat['donor'], nb_dat['recipient']))
    simulated_data = pd.read_csv(args.simulatedData, sep=',', header=None) 
    # //// Figure 1 ////
    figure_1(nb_dat=nb_dat, realign_vcf_dat=realign_vcf_dat, vcf_dat=vcf_dat, 
        ct_dat=ct_dat, transmission_pairs=transmission_pairs, family_dat=family_dat, 
        simulated_data=simulated_data, highlight_pair=args.highlightPair)
    # //// Figure 2 ////
    figure_2(nb_dat=nb_dat)
    # //// Figure S1 ////
    plot_tv_all(vcf_dat=realign_vcf_dat, transmission_pairs=transmission_pairs, ct_dat=ct_dat)
    # //// Figure S2 ////
    plot_shared_vars_all(vcf_dat=vcf_dat, transmission_pairs=transmission_pairs, min_af=0.01, min_abundance=3)
    # //// Figure S3 ////
    plot_var_abundance(vcf_dat=vcf_dat, transmission_pairs=transmission_pairs, min_af=0.01)


if __name__ == "__main__":
    run()
