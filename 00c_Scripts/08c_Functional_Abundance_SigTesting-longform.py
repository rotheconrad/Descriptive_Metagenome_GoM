#!/usr/bin/env python

''' Gene abundance correlation testing

Tests for correlation of shared gene abundance between samples.


Inputs:

    - Longform data table with columns: Gene, Sample, Annotation, Abundance
    - metadata mapping file with Sample column plus additional metadata cols
    - list of samples to test.
    - list of genes to test.

Outputs:

    - Scatter plots with correlation and p values for each gene.
    - Summary tsv file with columns: Gene, correlation, p value

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: February 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import pandas as pd
from scipy.stats import pearsonr
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(style="whitegrid")


# Set the default sans-serif font to Helvetica
# matplotlib.rcParams['font.sans-serif'] = "Helvetica"
# Set default font to sans-serif
# matplotlib.rcParams['font.family'] = "sans-serif"
# If you don't have Helvetica and need it see this:
# https://fowlerlab.org/2019/01/03/changing-the-sans-serif-font-to-helvetica/
# If you don't need Helvetica font just delete the above code.



def parse_list_file(list_file):
    # reads in a single column file to a list. Returns list
    flist = []
    with open(list_file, 'r') as file:
        for line in file:
            flist.append(line.rstrip())

    return flist


def test_genes_in_data(longform, gene_list):

    temp = {}
    genes_two = []

    with open(longform, 'r') as file:
        for line in file:
            X = line.split('\t')
            gene = X[0]
            temp[gene] = ''

    count = 0
    for gene in gene_list:
        if gene in temp:
            genes_two.append(gene)
        else:
            count += 1
            print('\t\t', gene)

    return genes_two, count


def add_metadata(longform, metadata, sample_list, gene_list):
    # read in files and create a longform dataform

    # Read the data matrix, stack it and reset the index
    m1 = pd.read_csv(longform, sep='\t')

    # Get list of genes to test
    genes_one = parse_list_file(gene_list)
    # test all genes in list are in the data
    genes, count = test_genes_in_data(longform, genes_one)
    # retreive annotations from longform data for genes in list
    annot_dict = pd.Series(m1.Annotation.values, index=m1.Gene)
    annots = set([annot_dict[gene] for gene in genes])

    # read the metadata and match to m1 using sample names
    m2 = pd.read_csv(metadata, sep='\t')
    for col in m2.columns:
        d1 = pd.Series(m2[col].values, index=m2.Sample).to_dict()
        d2 = []

        for S in m1.Sample.values:
            d2.append(d1[S])

        m1[col] = d2

    # Get list of samples to test
    samples = parse_list_file(sample_list)
    # subset data for samples to test
    print('\n\t\tNumber of samples in file:', len(samples))
    print('\t\tNumber of samples in data before subsampling:', len(m1))
    m1 = m1[m1['Sample'].isin(samples)]
    print('\t\tData rows after subsampling for sample:', len(m1))

    print('\n\t\tNumber of genes in file:', len(genes_one))
    print('\t\tNumber of genes from file NOT in longform data:', count)
    print('\t\tNumber of genes in file:', len(genes), '\n\n')

    return m1, annots


def clean_outfile_gene_name(string):

    sub_dict = {
                ',': '',
                ' ': '_',
                '(': '_',
                ')': '_',
                '[': '_',
                ']': '_',
                '+': '',
                '/': '_',
                '-': '_',
                '\'': '',
                ':': '_',
                '"': '',
                }

    new_string = ''

    for s in string:
        if s in sub_dict:
            new_string += sub_dict[s]
        else:
            new_string += s

    if "____" in new_string: new_string = new_string.replace("____", "_")

    if "___" in new_string: new_string = new_string.replace("___", "_")

    if "__" in new_string: new_string = new_string.replace("__", "_")

    if new_string.startswith('_'): new_string = new_string[1:]

    if new_string.endswith('_'): new_string = new_string[:-1]

    return new_string


def correlate_and_plot(df, outpre, expVar, resVar, annots):

    single_out = {'Gene': [], 'Pearson r': [], 'p value': []}
    sum_out = {'Gene': [], 'Pearson r': [], 'p value': []}


    # iterate through all annotations, compute correlation and plot
    for gene in annots:
        # Select only rows for current annotations
        single_df = df[df["Annotation"] == gene]
        sum_df = single_df.groupby(expVar).sum()
        
        # First we look all genes with the assigned annotation
        single_corr = pearsonr(single_df[expVar], single_df[resVar])
        single_r = single_corr[0]
        single_p = single_corr[1]

        single_out['Gene'].append(gene)
        single_out['Pearson r'].append(single_r)
        single_out['p value'].append(single_p)

        #print(f'{gene}\t{r}\t{p}')

        # Next we look at just the sum of expVar categories
        sum_corr = pearsonr(sum_df.index, sum_df[resVar])
        sum_r = sum_corr[0]
        sum_p = sum_corr[1]

        sum_out['Gene'].append(gene)
        sum_out['Pearson r'].append(sum_r)
        sum_out['p value'].append(sum_p)

        #print(f'{gene}\t{single_r}\t{single_p}\t{sum_r}\t{sum_p}')

        # Build plot if significant correlation
        #if abs(sum_r) >= 0.6 and sum_p <= 0.05:
        if abs(sum_r) >= 0.8 and sum_p <= 0.05:
        #if abs(sum_r) >= 0.2 and sum_p <= 0.2:
        #if 'transposase' in gene:
        #if sum_df[resVar].max() - sum_df[resVar].min() >= 1.0:
            g = sns.jointplot(
                            data=single_df,
                            x=expVar,
                            y=resVar,
                            kind="reg",
                            color="#2166ac",
                            label="Single Genes"
                            )
            sns.regplot(
                            data=sum_df,
                            x=sum_df.index,
                            y=resVar,
                            ax=g.ax_joint,
                            color="#b2182b",
                            label="Summed Genes"
                            )
            g.fig.suptitle(gene)
            stats_line = (
                f"Single Genes:\n"
                f"Pearson r: {single_r:.4f}\np value: {single_p:.4f}\n"
                f"Summed Genes:\n"
                f"Pearson r: {sum_r:.4f}\np value {sum_p:.4f}")
            g.ax_joint.text(
                0.55, .85, stats_line,
                verticalalignment='top', horizontalalignment='right',
                transform=g.ax_joint.transAxes
                )
            g.ax_joint.legend(loc='upper center')
            plt.tight_layout()
            outname = clean_outfile_gene_name(gene)
            print(f'\t\t{outname}')
            #if '/' in gene: gene = gene.replace('/', '-')

            plt.savefig(f'{outpre}_HistKDE_{outname}.pdf')
            plt.close()


    single_dfout = pd.DataFrame(single_out)
    single_dfout.to_csv(f'{outpre}_corrs_singles.tsv', sep='\t', index=False)
    sum_dfout = pd.DataFrame(sum_out)
    sum_dfout.to_csv(f'{outpre}_corrs_summed.tsv', sep='\t', index=False)

    return True
    

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-d1', '--longform_data_tsv',
        help='Longform data table with binary, count, or continuous values!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-d2', '--metadata_mapping_file',
        help='File containing sample names and metadata to test!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-d3', '--sample_test_list',
        help='List of samples to test.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-d4', '--rep_gene_test_list',
        help='List of cluster representative genes to test.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--prefix_for_output_files',
        help='Prefix to use for the output files!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-e', '--explanatory_variable_name',
        help='Column name for the explanatory variable to test',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-r', '--response_variable_name',
        help='Column name for the explanatory variable(s) to test',
        metavar='',
        type=str,
        required=True
        )

    args=vars(parser.parse_args())

    # Define input variable
    ldata = args['longform_data_tsv']
    mdata = args['metadata_mapping_file']
    sdata = args['sample_test_list']
    gdata = args['rep_gene_test_list']
    outpre = args['prefix_for_output_files']
    expVar = args['explanatory_variable_name']
    resVar = args['response_variable_name']

    # Do what you came here to do:
    print('\n\nRunning Script ...\n')

    print('\n\tReading in file and parsing data ...')
    df, annots = add_metadata(ldata, mdata, sdata, gdata)

    print('\n\tComputing correlations and building plots ... \n\n')
    _ = correlate_and_plot(df, outpre, expVar, resVar, annots)

    print('\n\nComplete success space cadet! Finished without errors.\n\n')


if __name__ == "__main__":
    main()
