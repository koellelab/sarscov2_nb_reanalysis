import argparse
import pandas as pd
from collections import defaultdict


def run():
    # formats the multiple metadata files from Popa et al. into a 
    # single file with one line per sample
    parser = argparse.ArgumentParser()
    # input files
    # todo outdir currently doesn't get used by anything
    parser.add_argument('--metadata', 
        default='data/abe2555_Data_file_S1_SampleInformation.csv',
        help='csv file with general metadata (Data File S1 Sample Information from Popa et al. ')
    parser.add_argument('--metadataLongitudinal', 
        default='data/abe2555_Data_file_S5.csv',
        help='csv file with longitudinal metadata (Data File S5 Sample Information from Popa et al. ')
    parser.add_argument('--metadataTransmission', 
        default='data/abe2555_Data_file_S4.csv',
        help='csv file with transmission pair metadata (Data File S4 Sample Information from Popa et al. ')
    args = parser.parse_args()
    metadata = pd.read_csv(args.metadata, sep=',', header=1)
    metadata = metadata[['Sample.Name', 'Sample.Date', 'Ct.Value']]
    long_metadata = pd.read_csv(args.metadataLongitudinal, sep=',')[['Sample', 'PatientIndex']]
    metadata =\
        metadata.merge(long_metadata, left_on='Sample.Name', 
            right_on='Sample', how='left')[['Sample.Name', 'Sample.Date', 'Ct.Value', 'PatientIndex']]
    transmission_metadata = pd.read_csv(args.metadataTransmission, sep=',')
    transmission_metadata['fromFamily_id'] = \
        transmission_metadata['fromFamily'].str.split('Family ', 
            expand=True)[1].str.split('/', expand=True)[0] 
    transmission_metadata['toFamily_id'] = \
        transmission_metadata['toFamily'].str.split('Family ', 
            expand=True)[1].str.split('/', expand=True)[0]     #
    # keys = sample names
    # values = indices of the transmission pairs for which
    #          this sample is the recippient
    donor_dict = defaultdict(list)
    recip_dict = defaultdict(list)
    family_dict = {}
    for i, row in transmission_metadata.iterrows():
        donor_dict[row['fromSample']].append(i)
        recip_dict[row['toSample']].append(i)
        family_dict[row['fromSample']] = row['fromFamily_id']
        family_dict[row['toSample']] = row['toFamily_id']
    metadata['donor_pairs'] = metadata['Sample.Name'].map(donor_dict)
    metadata['recip_pairs'] = metadata['Sample.Name'].map(recip_dict)
    metadata['family'] = metadata['Sample.Name'].map(family_dict)
    metadata['Sample.Name'] = \
        'CoV_' + metadata['Sample.Name'].str.replace('CeMM0', '').str.replace('CeMM', '')
    metadata.columns = \
        [['sample_name', 'date', 'ct', 'patient_index', 'donor_pairs', 'recip_pairs', 'family']]
    metadata.to_csv('data/abe255_Data_file_format.csv', sep=',', index=None)


if __name__ == "__main__":
    run()

