#!/usr/bin/env python

''' Statistical Testing for functional count and abundance matrix

This script runs on the two output matrices from 08a step.

Uses generalized linear model (GLM) from statsmodels to perform
ANOVA / permanova and posthoc tests.

Inputs:

    count or abundance matrix from 08a. (Dependent variable)
    Sample names with metadata. (Independent variables)

Converts data matrix to long form data for input into models.

Outputs:

    hmmmm.

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


def matrix_to_long(matrix, metadata):
    # read in files and create a longform dataform

    # Read the data matrix, stack it and reset the index
    m1 = pd.read_csv(matrix, sep='\t', index_col=0).stack().reset_index()
    # Set the column names
    m1.columns = ['Gene', 'Sample', 'Dependent_Value']

    # read the metadata and match to m1 using sample names
    m2 = pd.read_csv(metadata, sep='\t')
    for col in m2.columns:
        d1 = pd.Series(m2[col].values, index=m2.Sample).to_dict()
        d2 = []

        for S in m1.Sample.values:
            d2.append(d1[S])

        m1[col] = d2

    return m1


def exploratory_plots(df, outpre, dtype, expVar):

    for v in expVar:
        # Get unique values of expVar and sort them
        order = df[v].unique().tolist()
        order.sort()

        print('\n\t\t\t\tCreating Boxen plot ...')
        ax = sns.boxenplot(
                    x=v, y='Dependent_Value', order=order, color="b", data=df
                    )
        ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
        plt.tight_layout()
        plt.savefig(f'{outpre}_{dtype}_Boxen.pdf')
        plt.close()

        print('\n\t\t\t\tCreating Histogram-KDE plot ...')
        for D in order:
            dft = df[df[v] == D]
            sns.displot(
                data=dft, x='Dependent_Value', bins=30, kde=True, rug=True
                )
            plt.tight_layout()
            plt.savefig(f'{outpre}_HistKDE_{dtype}_{v}-{D}.pdf')
            plt.close()

    return True
    

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-d1', '--data_matrix',
        help='A data matrix with binary, count, or continuous values!',
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
        '-d3', '--data_matrix_type',
        help='Once of: binary, count, or continuous',
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
        '-e', '--Exploratory_Data_Analysis',
        help='(Optional) Build initial plots to inform analysis decisions!',
        action='store_true',
        required=False
        )
    parser.add_argument(
        '-v', '--Explanatory_Variable_Name',
        help='Column name for the explanatory variable to test',
        metavar='',
        nargs='+',
        type=str,
        required=True
        )

    args=vars(parser.parse_args())

    # Define input variable
    data = args['data_matrix']
    mdata = args['metadata_mapping_file']
    dtype = args['data_matrix_type']
    outpre = args['prefix_for_output_files']
    expVar = args['Explanatory_Variable_Name']

    # Do what you came here to do:
    print('\n\nRunning Script ...\n')

    print('\n\t\tReading in file and parsing data ...')
    df = matrix_to_long(data, mdata)

    if args['Exploratory_Data_Analysis']:
        print('\n\t\tBuilding exploratory plots ...')
        # If -E option build some distribution plots
        _ = exploratory_plots(df, outpre, dtype, expVar)

    print('\n\nComplete success space cadet! Finished without errors.\n\n')


if __name__ == "__main__":
    main()
