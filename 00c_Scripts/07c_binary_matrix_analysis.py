#!/usr/bin/env python

''' Analyse the binary matrix of samples and gene clusters.

Analyzes all combinations of the genomes/samples (columns) for shared
gene content vs unique gene content.

Careful this gets crazy with larger collections of genomes and samples.

Alternately, you can provide a list of combinations to analyze rather
than iterating through all combinations. For this file, make a plain
txt file with one combination to test per line. Separate genome/sample
names on each line with a single space. ex file:

LabelA, Name1 Name3 Name5 Name10
LabelB, Name4 Name6 Name7 Name8 Name12
LabelC, Name2 Name9 Name11

This script outputs tsv files.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Nov 2021
License :: GNU GPLv3
Copyright 2021 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, itertools, os
from collections import defaultdict
import pandas as pd


def parse_key_file(keyfile):

    key = {}

    with open(keyfile, 'r') as file:
        for line in file:
            X = line.rstrip().split(', ')
            cluster = X[0]
            geneName = X[1]
            key[cluster] = geneName

    return key


def parse_combination_list(lcom):

    # Parses the combinations file into a dictionary
    # initialize dict of {label: combination}
    # each combination is a list of column names to select from the df
    combination_dict = {}

    with open(lcom, 'r') as file:
        for line in file:
            X = line.rstrip().split(', ')
            label = X[0]
            combination = X[1:]
            combination_dict[label] = combination

    return combination_dict


def generate_combination_list(df):

    # initialize dict of {combination_number: combination}
    # each combination is a list of column names to select from the df
    combination_dict = {}
    # Get the list of column names from the df.
    name_list = df.columns.tolist()
    # iterate through names and create combinations of all lengths
    for i in range(2, len(name_list)+1):
        combinations = list(itertools.combinations(name_list, i))
        for n, comb in enumerate(combinations):
            combination_dict[f'combination_{n}'] = comb

    return combination_dict


def write_genes(key, clusters, out):

    with open(out, 'w') as file:
        for clust in clusters:
            file.write(f'{clust}, {key[clust]}\n')

    return True


def analyze_single(df, key, onam, outdir, core):

    # Calculate core vs accessory for all samples
    name_list = df.columns.tolist()
    n = len(name_list)
    specific = 1/n
    df['Sum'] = df.sum(axis=1)
    df['Fraction'] = df['Sum']/n
    df_core = df.loc[df['Fraction'] >= core]
    df_accs = df.loc[(df['Fraction'] < core) & (df['Fraction'] > specific)]
    df_spec = df.loc[df['Fraction'] == specific]

    # Calculations
    total_total = len(df)
    core_total = len(df_core)
    accs_total = len(df_accs)
    spec_total = len(df_spec)
    core_fract = core_total / total_total
    accs_fract = accs_total / total_total
    spec_fract = spec_total / total_total

    # get accessory gene counts for each sample
    accs_sum = df_accs.sum()
    accs_ndx = accs_sum.index.tolist()
    accs_sms = accs_sum.tolist()
    accs_lst = [f'{ndx}:{accs_sms[i]}' for i, ndx in enumerate(accs_ndx)]
    accs_brkd = ', '.join(accs_lst)
    
    # get specific gene counts for each sample
    spec_sum = df_spec.sum()
    spec_ndx = spec_sum.index.tolist()
    spec_sms = spec_sum.tolist()
    spec_lst = [f'{ndx}:{spec_sms[i]}' for i, ndx in enumerate(spec_ndx)]
    spec_brkd = ', '.join(spec_lst)

    # Print summary results
    print(f'Total Genes:\t{total_total}')
    print(f'Core Total:\t{core_total}')
    print(f'Core Fraction:\t{core_fract}')
    print(f'Accessory Total:\t{accs_total}')
    print(f'Accessory Fraction:\t{accs_fract}')
    print(f'Specific Total:\t{spec_total}')
    print(f'Specific Fraction:\t{spec_fract}')
    print(f'Accessory Breakdown:\t{accs_brkd}')
    print(f'Specific Breakdown:\t{spec_brkd}')

    ''' I don't think I need anymore. They were diagnostic
    # Write sub-dataframes for each combination
    df_core.to_csv(f'{outdir}/{onam}_ALL_core.tsv', sep='\t')
    df_accs.to_csv(f'{outdir}/{onam}_ALL_accs.tsv', sep='\t')
    df_spec.to_csv(f'{outdir}/{onam}_ALL_spec.tsv', sep='\t')
    '''

    # Write gene cluster names for each category.
    _ = write_genes(key, list(df_core.index.values), f'{onam}_ALL_core.csv')
    _ = write_genes(key, list(df_accs.index.values), f'{onam}_ALL_accs.csv')
    _ = write_genes(key, list(df_spec.index.values), f'{onam}_ALL_spec.csv')


def analyze_multiple(df, key, onam, outdir, combination_dict, core):

    # initialize dict for summary data
    sumdata = {
            'Combination': [],
            'Total_Genes': [],
            'Core_Genes': [],
            'Core_Fraction': [],
            'Accessory_Genes': [],
            'Accessory_Fraction': [],
            'Specific_Genes': [],
            'Specific_Fraction': [],
            'Samples/Genomes': [],
            'Accessory-Breakdown (Sample:GeneCount)': [],
            'Specific-Breakdown (Sample:GeneCount)': [],
            }

    # Calculate core vs accessory for each combination in list
    for i, (label, combi) in enumerate(combination_dict.items()):
        print(f'Computing combination: {label} ...')
        # Select columns for current combination
        df2 = df[combi].copy()
        # Get length of current combination
        n = len(combi)
        # Set fraction for specific genes
        specific = 1/n
        # Sum each row
        df2['Sum'] = df2.sum(axis=1)
        # Drop rows equal to zero as those genes are not present in the
        # current sample combination
        df2 = df2.loc[df2['Sum'] != 0]
        # Fraction of total for each row
        df2['Fraction'] = df2['Sum']/n
        # Select core genes
        df_core = df2.loc[df2['Fraction'] >= core]
        # Select accessory genes
        df_accs = df2.loc[
                            (df2['Fraction'] < core) & 
                            (df2['Fraction'] > specific)
                            ]
        # Select specific genes
        df_spec = df2.loc[df2['Fraction'] == specific]

        # Calculations
        total_total = len(df2)
        core_total = len(df_core)
        accs_total = len(df_accs)
        spec_total = len(df_spec)
        core_fract = round(core_total / total_total, 4)
        accs_fract = round(accs_total / total_total, 4)
        spec_fract = round(spec_total / total_total, 4)

        # get accessory gene counts for each sample
        accs_sum = df_accs.drop(['Sum', 'Fraction'], axis=1).sum()
        accs_ndx = accs_sum.index.tolist()
        accs_sms = accs_sum.tolist()
        accs_lst = [f'{ndx}:{accs_sms[i]}' for i, ndx in enumerate(accs_ndx)]
        accs_brkd = ', '.join(accs_lst)
        
        # get specific gene counts for each sample
        spec_sum = df_spec.drop(['Sum', 'Fraction'], axis=1).sum()
        spec_ndx = spec_sum.index.tolist()
        spec_sms = spec_sum.tolist()
        spec_lst = [f'{ndx}:{spec_sms[i]}' for i, ndx in enumerate(spec_ndx)]
        spec_brkd = ', '.join(spec_lst)

        # Store data for summary table
        sumdata['Combination'].append(label)
        sumdata['Total_Genes'].append(total_total)
        sumdata['Core_Genes'].append(core_total)
        sumdata['Core_Fraction'].append(core_fract)
        sumdata['Accessory_Genes'].append(accs_total)
        sumdata['Accessory_Fraction'].append(accs_fract)
        sumdata['Specific_Genes'].append(spec_total)
        sumdata['Specific_Fraction'].append(spec_fract)
        sumdata['Samples/Genomes'].append(', '.join(combi))
        sumdata['Accessory-Breakdown (Sample:GeneCount)'].append(accs_brkd)
        sumdata['Specific-Breakdown (Sample:GeneCount)'].append(spec_brkd)

        # set outfile prefix
        outpre = f'{outdir}/{onam}_{i+1:03}_{label}'
        ''' I don't think I need anymore. They were diagnostic
        # Write sub-dataframes for each combination
        df_core.to_csv(f'{outpre}_core.tsv', sep='\t')
        df_accs.to_csv(f'{outpre}_accs.tsv', sep='\t')
        df_spec.to_csv(f'{outpre}_spec.tsv', sep='\t')
        '''
        # Write gene cluster names for each category.
        _ = write_genes(key, list(df_core.index.values), f'{outpre}_core.csv')
        _ = write_genes(key, list(df_accs.index.values), f'{outpre}_accs.csv')
        _ = write_genes(key, list(df_spec.index.values), f'{outpre}_spec.csv')

    # Convert summary data to df and write to tsv
    df_summary = pd.DataFrame(sumdata)
    df_summary.to_csv(f'{outdir}/{onam}_000_Summary.tsv', sep='\t', index=False)


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-b', '--binary_matrix_file',
        help='Please specify binary matrix file for input!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-k', '--GeneCluster_RepGeneName_key',
        help='Please specify the cluster-key.tsv file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-d', '--output_directory',
        help='What directory would you like to write files to?',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-p', '--output_file_prefix',
        help='What would you like to call the output files?',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-l', '--list_of_combinations',
        help='OPTIONAL: List of combinations to analyze.',
        metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-a', '--all_samples_only',
        help='OPTIONAL: Ignore sub-combinations and only compare all samples.',
        action="store_true",
        default=False,
        required=False
        )
    parser.add_argument(
        '-c', '--set_core_threshold',
        help='OPTIONAL: Set ratio to consider gene core (default = 1.0).',
        metavar='',
        type=float,
        required=False,
        default=1.0
        )
    args=vars(parser.parse_args())

    # Run this scripts main function
    print('\n\nRunning Script...\n\n')

    # Store the input variables
    bmat = args['binary_matrix_file']
    keyfile = args['GeneCluster_RepGeneName_key']
    odir = args['output_directory']
    onam = args['output_file_prefix']
    lcom = args['list_of_combinations']
    core = args['set_core_threshold']
    alls = args['all_samples_only']

    # Create output directory if it does not exist
    os.makedirs(odir, exist_ok=True)
    
    # Read in the binary matrix
    df = pd.read_csv(bmat, sep='\t', index_col=0, header=0)
    # parse the key file
    key = parse_key_file(keyfile)
    # generate/analyse sample combinations
    if lcom:
        combination_dict = parse_combination_list(lcom)
        _ = analyze_multiple(df, key, onam, odir, combination_dict, core)
    elif alls:
        _ = analyze_single(df, key, onam, odir, core)
    else:
        combination_dict = generate_combination_list(df)
        _ = analyze_multiple(df, key, onam, odir, combination_dict, core)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')


if __name__ == "__main__":
    main()