import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages



def process_allele_freqs(pairs_path=None, donor_col='donor', recipient_col='recipient',
    allele_freqs_path=None, min_allele_freq=0.0, snp_only=True):
    transmission_pairs = pd.read_csv(pairs_path, sep=',')
    pairs = \
     [(row[donor_col], row[recipient_col]) for idx, row in transmission_pairs.iterrows()]
    allele_freqs = pd.read_csv(allele_freqs_path, sep='\t', low_memory=False)
    # subset just to allele freq data
    allele_freqs = allele_freqs[allele_freqs['Type'] == 'AF']
    # subset just to samples we have data for (those on SRA)
    allele_freqs = allele_freqs[~allele_freqs['Sample'].isnull()]
    # subsets to just those alleles > min_allele_freq
    allele_freqs['value'] = allele_freqs['value'].astype(float)
    allele_freqs = allele_freqs[allele_freqs['value'] >= min_allele_freq]
    # subsets to just SNPs
    if snp_only:
        allele_freqs = allele_freqs[(allele_freqs['REF'].str.len() == 1) & (allele_freqs['ALT'].str.len()==1)]
    tv_plot_data = pd.DataFrame()
    for idx, pair in enumerate(pairs):
        pair_allele_freqs = allele_freqs[allele_freqs['Sample'].isin(pair)][['POS', 'REF', 'ALT', 'Sample', 'value']]
        pair_allele_freqs = \
            pair_allele_freqs.pivot(index=['POS', 'REF', 'ALT'], columns='Sample', values='value').reset_index()
        pair_allele_freqs = \
            pair_allele_freqs.rename(columns={pair[0]: 'fromFreq', pair[1]: 'toFreq'})
        pair_allele_freqs['from'] = pair[0]
        pair_allele_freqs['to'] = pair[1]
        pair_allele_freqs['pairIdx'] = idx
        tv_plot_data = pd.concat([tv_plot_data, pair_allele_freqs])
    tv_plot_data.to_csv('tv_plot_data.csv', index=None)
    return(tv_plot_data)



def make_tv_plots(dat=None):
    pairs = set(dat['pairIdx'])
    with PdfPages('tv_plots.pdf') as pdf:
        for pair_idx in pairs:
            pair_dat = dat[dat['pairIdx'] == pair_idx]
            fig, axs = plt.subplots(1,2, figsize=(6.4*2, 4.8), constrained_layout=True)
            for ax in axs:
                ax.scatter(pair_dat['fromFreq'], pair_dat['toFreq'])
                ax.set_xlabel(pair_dat.iloc[0]['from'])
                ax.set_ylabel(pair_dat.iloc[1]['to'])
            axs[0].set_xlim(-0.05, 1.05)
            axs[0].set_ylim(-0.05, 1.05)
            axs[0].set_xticks([0.00, 0.25, 0.50, 0.75, 1.00])
            axs[0].set_yticks([0.00, 0.25, 0.50, 0.75, 1.00])
            axs[1].set_xlim(-0.005, 0.055)
            axs[1].set_ylim(-0.005, 0.055)
            axs[1].set_xticks([0.00, 0.01, 0.02, 0.03, 0.04, 0.05])
            axs[1].set_yticks([0.00, 0.01, 0.02, 0.03, 0.04, 0.05])
            fig.savefig('test.pdf')


def run():
parser = argparse.ArgumentParser()
# input files
# todo outdir currently doesn't get used by anything
parser.add_argument('--transmission_pairs', 
 default='data/abe2555_data_file_S4.csv',
 help='file with information on transmission pairs')
parser.add_argument('--donor', 
 default='fromSample',
 help='column name with donor sample name')
parser.add_argument('--recipient',
 default='toSample',
 help='column name with recipient sample name')
parser.add_argument('--allele_freqs',
 help='file with allele frequencies (assumes long format data with columns: POS, Sample, value, Type', 
 default='data/SARSCoV_allSamples_split_AF001_long_format.tab')
parser.add_argument('--min_allele_freq',
 help='minum allele frequency for a variant to be considered', 
 type=float,
 default=0.0)
parser.add_argument('--snp_only',
 help='plot just SNPs or SNPs and indelx', 
 type=bool,
 default=True)
args = parser.parse_args()
tv_plot_data = process_allele_freqs(
    pairs_path=args.transmission_pairs, donor_col=args.donor, recipient_col=args.recipient,
    allele_freqs_path=args.allele_freqs, min_allele_freq=args.min_allele_freq, snp_only=args.snp_only)
make_tv_plots(tv_plot_data)
# todo add seps as arguments



