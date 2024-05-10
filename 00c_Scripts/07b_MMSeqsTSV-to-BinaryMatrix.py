#!/usr/bin/env python

''' Convert MMSeqs2 Cluster TSV file to binary matrix

The binary matrix stores gene presence/absence data for each sample.
Gene clusters are the rows and samples are the columns.
Zero indicates the gene cluster is absent from a sample.

input: MMSeqs2 Cluster TSV file.
output: Binary matrix as TSV file.

Sample names are inferred from the gene sequence names.
Modify line 93 sample variable to select the right thing.

The cluster key output file just links the cluster number assigned
to rows of the binary matrix output to the representative sequence
for that cluster.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: December, 2021
License :: GNU GPLv3
Copyright 2021 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from collections import defaultdict
import pandas as pd


def generate_binary_matrix(byCluster, bySample):

    '''
    This function takes two dicionarys with cluster data stored
    by cluster representative name and by sample name.
    It iterates throught the set of cluster representative names
    for each sample entry and checks to see if the sample is present
    or absent in that gene cluster. It encodes a binary datatable with
    genes clusters as rows and samples as columns. Zero indicates a 
    gene cluster is absent from a sample.
    '''

    print('\n\nGenerating binary matrix ...')
    #initialize variables
    # get list of all the cluster ids
    cid = list(byCluster.keys())
    # create new list to use for index. cleaner names: Cluster_#
    index = [f'Cluster_{i}' for i, c in enumerate(cid)]
    # store binary encoding in dict of {Sample: [binary]}
    matrix = defaultdict(list)
    # read through data and encode binary info
    for sample, clusters in bySample.items():
        for c in cid:
            if c in clusters: matrix[sample].append(1)
            else: matrix[sample].append(0)
    # create dataframe from binary matrix encoding
    binary_matrix = pd.DataFrame(matrix, index=index)
    # create key dictionary of {index: cid}
    key = dict(zip(index, cid))

    return binary_matrix, key


def parse_mmseqs_cluster_tsv(infile):

    '''
    The mmseqs cluster tsv file is two columns.
    column 1 is the cluster representative sequence.
    This sequence is repeated for every sequence in the cluster.
    Column 2 are sequences in the cluster.
    The general idea with this function is to create 2 dicionaries.
    One dictionary stores data by the cluster represenative name as in
    {cluster: sample}. Which samples are in each cluster.
    The other dicitonary stores data by the sample name.
    {sample: cluster}. Which clusters are in each sample.
    '''

    print('\n\nParsing MMSeqs2 Cluster TSV file ...')
    # initialize variables
    # Store data in two dictionaries.
    byCluster = defaultdict(lambda: defaultdict(int))
    bySample = defaultdict(lambda: defaultdict(int))

    # read through the file and populate the dictionaries
    with open(infile, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            cluster = X[0]
            sample = 'EN' + X[1].split('_')[1]
            # store data in dict of dict
            byCluster[cluster][sample] += 1
            bySample[sample][cluster] += 1

    return byCluster, bySample


def main():
    bySample = defaultdict(lambda: defaultdict(int))

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--mmseqs_cluster_tsv_input_file',
        help='Please specify the mmseqs cluster tsv input file!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--binary_matrix_tsv_output_file',
        help='Please specify the name to use for the output file!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-k', '--GeneCluster_RepGeneName_key',
        help='Please specify output name for the cluster key!',
        metavar=':',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # define the input parameters
    infile = args['mmseqs_cluster_tsv_input_file']
    outfile = args['binary_matrix_tsv_output_file']
    keyout = args['GeneCluster_RepGeneName_key']

    # parse the input file
    byCluster, bySample = parse_mmseqs_cluster_tsv(infile)
    # generate the binary matrix
    binary_matrix, key = generate_binary_matrix(byCluster, bySample)
    # write the binary matrix to a tsv file
    binary_matrix.to_csv(outfile, sep='\t', index=True, header=True)
    # write the key to a csv file
    with open(keyout, 'w') as file:
        for cluster, name in key.items():
            file.write(f'{cluster}, {name}\n')

    print('\n\nComplete success space cadet! Finished without errors.\n\n')

if __name__ == "__main__":
    main()

