#!/usr/bin/env python

'''Plot histograms from tabular Blast.

Plots histograms from tabular blast file for percent identity, alignment
length, and percent match length (alignment lenght / query sequence length)

This tool takes the following input parameters:

    * Sequences.blast - tabular Blast File

This script returns the following files:

    * 3 Histogram Plots .pdf format

This script requires the following packages:

    * argparse
    * matplotlib

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: January 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''


import argparse
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt


def parse_blast(file):
    """ parse blast file for pidents returns list of floats """

    data = {
            'pid': [],
            'qalen': [],
            'ralen': [],
            'pml': [],
            'qrl': [],
            'nalign': []
            }

    with open(file, 'r') as f:
        for l in f:
            if l.startswith('#'): continue # skip header lines
            X = l.rstrip().split('\t')
            pident = float(X[2])

            qastart = int(X[6]) # query seq alignment start
            qastop = int(X[7]) # query seq alignment stop
            rastart = int(X[8]) # ref seq alignment start
            rastop = int(X[9]) #ref seq alignment stop

            qalen = max([qastart, qastop]) - min([qastart, qastop])
            ralen = max([rastart, rastop]) - min([rastart, rastop])
            pml = qalen / ralen
            qrl = float(X[15]) # Query/read length
            nalignments = int(X[17]) # number of alignments reported for query


            data['pid'].append(pident)
            data['qalen'].append(qalen)
            data['ralen'].append(ralen)
            data['pml'].append(pml)
            data['qrl'].append(qrl)
            data['nalign'].append(nalignments)

    return data


def plot_hist(data, outpre, key):

    # Define the titles
    plot_titles = {
                    'pid': 'Sequence Identity Histogram (%)',
                    'qalen': 'Query Alignment Length Histogram (bp)',
                    'ralen': 'Reference ALignment Length Histogram (bp)',
                    'pml': 'Query/Reference Alignment Ratio (%)',
                    'qrl': 'Query to Query Alignment Ratio (%)',
                    'nalign': 'Alignments Per Query Read'
                    }

    print(f'\t\tPlotting {plot_titles[key]}...')

    # Set the colors
    bar_color = '#2171b5'
    gridM = '#bdbdbd'
    gridm = '#d9d9d9'
    alpha = 0.6

    # Build the plot
    fig, ax = plt.subplots(figsize=(10, 7))

    # Plot titles
    ax.set_title(
        plot_titles[key],
        fontsize=20, y=1.02
        )

    # Plot labels
    ax.set_ylabel('Count', fontsize=14)
    ax.set_xlabel('Value', fontsize=14)

    # Set plot/grid style
    ax.minorticks_on()
    ax.tick_params(
        which='minor', axis='both', left=False, bottom=False
        )
    ax.tick_params(
                which='major', axis='both',
                left=True, bottom=True,
                size=6, width=2, tickdir='inout',
                labelsize=12, zorder=10
                )
    ax.yaxis.grid(
        which="minor", color=gridm, linestyle='--',
        linewidth=1, alpha=0.6, zorder=1
        )
    ax.yaxis.grid(
        which="major", color=gridM, linestyle='--',
        linewidth=1.5, alpha=0.4, zorder=1
        )
    ax.set_axisbelow(True)
    for spine in ax.spines.values(): spine.set_linewidth(2)

    # Plot the data
    ax.hist(
        data[key],
        bins=30,
        #orientation='horizontal',
        rwidth=0.9,
        color=bar_color,
        alpha=alpha,
        )

    # Set plot axis ranges
    #ax.set_xlim(left=0, right=int((max(d['xs'])+min(d['xs']))))

    # adjust layout, save, and close
    #plt.gca().invert_yaxis()
    fig.set_tight_layout(True)
    plt.savefig(f'{outpre}_{key}_histogram.pdf')
    plt.close()


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the input file name!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_prefix',
        help='Please specify the prefix to use for output files!',
        metavar=':',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\nRunning Script ...\n')
    
    print('\tParsing tabular blast ...\n')
    data = parse_blast(args['input_file'])

    print('\tFinished parsing. Building plots ...\n')
    for key in data.keys():
        _ = plot_hist(data, args['output_prefix'], key)

    print('\n\nComplete success space cadet! Hold on to your boots.\n\n')


if __name__ == "__main__":
    main()
