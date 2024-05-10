#!/usr/bin/env python

'''Generate Whole Genome Summary Table Across Metagenome Samples.

This script normalizes the TAD80 by the metagenome size in
Giga-Base-Pairs so TAD80 is comparable across metagenomes. This means the
TAD80 reported by this script is calculated as:

        TAD80 / length of metagenome (Gbps) making the units 1/Gbps

This script also optionally normalizes by genome equivalents from
microbe census outputs if they are provided.

ANIrs and Relative Abundance should be comparable without normalization.

This tool takes the following input parameters:

    * input dir - containing *_genome.tsv files
    * out - output file name

This script returns the following files:

    * Summary Data Table

This script requires the following packages:

    * argparse
    * os.listdir
    * os.path.isfile
    * os.path.join
    * collections.defaultdict

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: February 5th, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from os import listdir
from os.path import isfile, join
from collections import defaultdict


def gather_geqs(geq_dir):
    ''' get the geq results for each sample '''

    file_list = [f for f in listdir(geq_dir) if isfile(join(geq_dir, f))]

    geqs = {}

    for file in file_list:
        with open(f'{geq_dir}/{file}', 'r') as f:
            smpl = file.split('.')[0]
            for line in f:
                X = line.rstrip().split('\t')
                if X[0] == 'genome_equivalents:':
                    geq = float(X[1])
            geqs[smpl] = geq

        print(smpl, geq)

    return geqs


def gather_data(gtd, outfile, geqs):
    ''' Read in the files and organize the data '''

    file_list = [f for f in listdir(gtd) if isfile(join(gtd, f))]

    tads = defaultdict(lambda: defaultdict(list))
    gqed = defaultdict(lambda: defaultdict(list))
    breadths = defaultdict(lambda: defaultdict(list))
    anirs = defaultdict(lambda: defaultdict(list))
    rels = defaultdict(lambda: defaultdict(list))
    sample_dict = {}

    for file in file_list:
        with open(f'{gtd}/{file}', 'r') as f:
            header = f.readline()
            X = f.readline().rstrip().split('\t')

            #############################################################
            ####### ADJUST THIS BLOCK FOR YOUR FILE NAMING SCHEME #######
            #############################################################
            file_basename = X[0].split('-')
            smpl = file_basename[0]
            sample_dict[smpl] = ''
            genome = file_basename[1]
            #############################################################
            #############################################################

            tad = float(X[1])
            breadth = X[2]
            anir = X[3][:-1]
            relabnd = X[4][:-1]
            metag_size = X[6]

            tads[genome][smpl] = str(tad / int(metag_size) * 1000000000)
            if geqs:
                gqed[genome][smpl] = str(tad / geqs[smpl])
            breadths[genome][smpl] = breadth
            anirs[genome][smpl] = anir
            rels[genome][smpl] = relabnd


    sample_list = list(sample_dict.keys())
    header = 'Genome\t' + '\t'.join(sample_list) + '\n'

    with open(f'{outfile}_tads.tsv', 'w') as o:

        o.write(header)
        for genome, samples in tads.items():
            lineout = [genome]

            for s in sample_list:
                value = samples[s]
                if isinstance(value, list): value = '0.0'
                lineout.append(value)

            o.write('\t'.join(lineout) + '\n')

    if geqs:
        with open(f'{outfile}_gqed.tsv', 'w') as o:

            o.write(header)
            for genome, samples in gqed.items():
                lineout = [genome]

                for s in sample_list:
                    value = samples[s]
                    if isinstance(value, list): value = '0.0'
                    lineout.append(value)

                o.write('\t'.join(lineout) + '\n')

    with open(f'{outfile}_breadths.tsv', 'w') as o:

        o.write(header)
        for genome, samples in breadths.items():
            lineout = [genome]

            for s in sample_list:
                value = samples[s]
                if isinstance(value, list): value = '0.0'
                lineout.append(value)

            o.write('\t'.join(lineout) + '\n')

    with open(f'{outfile}_anirs.tsv', 'w') as o:

        o.write(header)
        for genome, samples in anirs.items():
            lineout = [genome]

            for s in sample_list:
                value = samples[s]
                if isinstance(value, list): value = '0.0'
                lineout.append(value)

            o.write('\t'.join(lineout) + '\n')

    with open(f'{outfile}_rels.tsv', 'w') as o:

        o.write(header)
        for genome, samples in rels.items():
            lineout = [genome]

            for s in sample_list:
                value = samples[s]
                if isinstance(value, list): value = '0.0'
                lineout.append(value)
                
            o.write('\t'.join(lineout) + '\n')

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-gtd', '--genome_tsv_dir',
        help='Please specify the directory with *_genome.tsv files!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-out', '--output_file_prefix',
        help='Please specify the name to us for the output file prefix.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-gq', '--geq_dir',
        help='(Optional): Directory of Microbe Census results.',
        metavar='',
        type=str,
        required=False
        )
    args=vars(parser.parse_args())

    # Run this scripts main function
    print('\n\nParsing data into summary table ...\n\n')

    fdir = args['genome_tsv_dir']
    outfile = args['output_file_prefix']
    geq_dir = args['geq_dir']

    if fdir[-1] == '/': fdir = fdir[:-1]
    if geq_dir[-1] == '/': geq_dir = geq_dir[:-1]

    if geq_dir:
        geqs = gather_geqs(geq_dir)
    else:
        geqs = None

    _ = gather_data(fdir, outfile, geqs)


if __name__ == "__main__":
    main()