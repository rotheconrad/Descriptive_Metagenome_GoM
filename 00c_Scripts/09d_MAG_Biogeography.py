#!/usr/bin/env python

'''Build heatmap of MAG Biogeography from Coverage Magic Combined Stats

Takes a tsv file of Genome Summary files from 09c_CoverageMagic.py
compiled by 09c_CoverageMagic_CombineGenomeStats.py.

Can be coverage (TAD80) normailzed by sequencing effort of genome
equivalents, breadths, ANIr, or Relative Abundance tsv file.

Writes out heatmap of results in PDF format. Publication ready.

I don't like heatmaps that are too busy. This script sums the rows
(that is genome values across samples) and plots the top 50 rows.
This can be changed with the -n flag.

### Future improvement. Take metadata file with samples by category
and a second file with the colors to use for each category.
Map colors to each category and write out the legends.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Jan 04 2021
License :: GNU GPLv3
Copyright 2021 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, re
from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd
import matplotlib
from matplotlib.patches import Patch
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Set the default sans-serif font to Helvetica
#matplotlib.rcParams['font.sans-serif'] = "Helvetica"
# Set default font to sans-serif
#matplotlib.rcParams['font.family'] = "sans-serif"
# If you don't have Helvetica and need it see this:
# https://fowlerlab.org/2019/01/03/changing-the-sans-serif-font-to-helvetica/
# If you don't need Helvetica font just delete the above code.


def build_heatmaps(data, metadata, metacolors, outpre, sorder, classifications):

    color = 'Greys'

    # if metadata and colors read in and prepare the data.
    if metadata and metacolors:
        metad = pd.read_csv(metadata, sep='\t', header=0, dtype = str)
        metad.sort_values(by='Sample', inplace=True)
        # The indexes need to match metad to df.
        metad.set_index('Sample', inplace=True)
        if sorder: metad = metad.rename(index=sorder)
        # make a dict from Labels: Colors and map to metad values.
        metac = pd.read_csv(metacolors, sep='\t', header=0, dtype = str)
        metac_dict = dict(zip(metac['Labels'], metac['Colors']))
        mappedmeta = metad.stack().map(metac_dict).unstack()

    cbarpos = (0.05, .9, .2, .05) #(left, bottom, width, height),
    plotsize = (25,12)


    dflist = [i for i in data.to_numpy().flatten() if i != 0]
    low = min(dflist)
    high = max(dflist)
    print(f'\t\t\t\tlow: {low}\n\t\t\t\thigh: {high}')

    # plot clustermap
    if metadata and metacolors:
        g = sns.clustermap(
                    data, figsize=plotsize, cmap=color, vmin=low, vmax=high,
                    cbar_pos=cbarpos, row_colors=mappedmeta,
                    cbar_kws={"orientation": "horizontal"},
                    yticklabels=True, xticklabels=True
                    )
    else:
        g = sns.clustermap(
                    data, figsize=plotsize, cmap=color,
                    cbar_pos=cbarpos, vmin=low, vmax=high,
                    cbar_kws={"orientation": "horizontal"},
                    yticklabels=True, xticklabels=True
                    )

    plt.savefig(f'{outpre}_clustermap.pdf', dpi=300)
    plt.close()

    # plot heatmap
    if metadata and metacolors:
        #sns.set(font_scale=0.5)
        g = sns.clustermap(
                    data, figsize=plotsize, cmap=color, vmin=low, vmax=high,
                    cbar_pos=cbarpos, row_colors=mappedmeta,
                    cbar_kws={"orientation": "horizontal"},
                    row_cluster=False, xticklabels=True
                    )
    else:
        g = sns.clustermap(
                    data, figsize=plotsize, cmap=color,
                    cbar_pos=cbarpos, vmin=low, vmax=high,
                    cbar_kws={"orientation": "horizontal"},
                    row_cluster=False, xticklabels=True
                    )

    # Print Mag numbers
    plt.savefig(f'{outpre}_heatmap.pdf', dpi=300)
    plt.close()

    if classifications:
        xlabo = [i.get_text() for i in g.ax_heatmap.get_xticklabels()]
        switch = [classifications[i] for i in xlabo]
        with open(f'{outpre}_Classifications.txt', 'w') as clssfy_out:
                clssfy_out.writelines(f'{clss}\n' for clss in switch)

    ## plot meta data legends #########################################
    if metadata and metacolors:
        for meta in metad.columns:
            labels = natural_sort(metad[meta].unique())

            fig, ax = plt.subplots(figsize=(10,10))
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)

            for label in labels:
                ax.bar(
                    0, 0, color=metac_dict[label], label=label, linewidth=0
                    )

            ax.legend(
                title=meta, title_fontsize='xx-large', loc="center",
                frameon=False, markerscale=5, fontsize='xx-large'
                )

            plt.savefig(f'{outpre}_Legend_{meta}.pdf', dpi=300)
            plt.close()
    ###################################################################

def natural_sort(l): 
    # copied from:
    # https://stackoverflow.com/questions/11150239/python-natural-sorting
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)
    

def get_order(samp_order):
    # read in sample order list
    sorder = {}

    with open(samp_order, 'r') as file:
        header = file.readline()
        for line in file:
            X = line.rstrip().split('\t')
            samp = X[0]
            order = int(X[-1])
            sorder[samp] = order

    return sorder

def get_classifications(infile):
    # parse MAG Classification file
    classifications = {}

    with open(infile, 'r') as file:
        header = file.readline()
        for line in file:
            X = line.rstrip().split('\t')
            MAG = X[0]
            classify = ' '.join(X[1].split('_')[:2])
            AAI = float(X[2])
            classifications[MAG] = f'{MAG}\t{classify}\t({AAI:.2f})'

    return classifications

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the input tsv file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-md', '--metadata',
        help='OPTIONAL: Metadata to color samples by! requires -mc',
        metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-mc', '--metacolors',
        help='OPTIONAL: Metadata to color samples by! requires -md',
        metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-o', '--output_prefix',
        help='Please specify the filename prefix for the plots!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-v', '--value_threshold',
        help='Optional: Set lower bound cutoff for tsv values (default: 0.0)',
        metavar='',
        type=float,
        required=False,
        default=10.0
        )
    parser.add_argument(
        '-s', '--sample_threshold',
        help='Optional: Set sample threshold detection limit (default>=2)',
        metavar='',
        type=int,
        required=False,
        default=2
        )
    parser.add_argument(
        '-r', '--sample_reorder',
        help='Optional: Specify the order of samples in the heatmap',
        metavar='',
        type=str,
        required=False,
        )
    parser.add_argument(
        '-c', '--MAG_classification',
        help='Optional: Specify the Classification for MAGs',
        metavar='',
        type=str,
        required=False,
        )
    parser.add_argument(
        '-n', '--rows_to_plot',
        help='Optional: Set number of rows to plot (default: 50)',
        metavar='',
        type=int,
        required=False,
        default=50
        )
    args=vars(parser.parse_args())
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script ...\n')

    # Check for sample order file and process if true
    sorder = None
    classifications = None
    if args['sample_reorder']: sorder = get_order(args['sample_reorder'])

    if args['MAG_classification']:
        classifications = get_classifications(args['MAG_classification'])

    # Collect stats data
    print('\n\t\tReading files ... ')
    data = pd.read_csv(args['input_file'], sep='\t')
    data = data.set_index('Genome')
    data['Sum'] = data.sum(axis=1)
    data = data.sort_values('Sum', ascending=False).head(args['rows_to_plot'])
    data = data.drop('Sum', axis=1).T
    print(data)
    

    print('\n\t\tPlotting ... ')
    # Plot stats data
    build_heatmaps(
                data,
                args['metadata'],
                args['metacolors'],
                args['output_prefix'],
                sorder,
                classifications
                )

    print('\n\nComplete success space cadet!! We finished without errors.\n\n')


if __name__ == "__main__":
    main()
