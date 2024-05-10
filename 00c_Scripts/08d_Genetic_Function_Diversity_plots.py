#!/usr/bin/env python

'''
Builds plots from the *_Summary.tsv output of 07c_binary_matrix_analysis.py

Used to compare clustering results from multiple differeing clustering
parameters such as nucleotides at 95% and 90% sequence similarity and
amino acids at 70% and 40% sequence similarity.

input: Directory of *_Summary.tsv files
output: Plots in pdf format.

* Name name of the *_Summary.tsv is used in the plots
* Name example:
    * 95_nucleotide_Summary.tsv
    * 90_nucleotide_Summary.tsv
    * 70_aminoacid_Summary.tsv
    * 40_aminoacid_Summary.tsv

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: April 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from os import listdir
from os.path import isfile, join
from collections import defaultdict
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Set the default sans-serif font to Helvetica
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
# Set default font to sans-serif
matplotlib.rcParams['font.family'] = "sans-serif"

class TickRedrawer(matplotlib.artist.Artist):
    # this is just to get the stupid ticks to draw right on the plot
    #https://stackoverflow.com/questions/19677963/
    #matplotlib-keep-grid-lines-behind-the-graph-but-the-y-and-x-axis-above
    """Artist to redraw ticks."""

    __name__ = "ticks"

    zorder = 10

    @matplotlib.artist.allow_rasterization
    def draw(self, renderer: matplotlib.backend_bases.RendererBase) -> None:
        """Draw the ticks."""
        if not self.get_visible():
            self.stale = False
            return

        renderer.open_group(self.__name__, gid=self.get_gid())

        for axis in (self.axes.xaxis, self.axes.yaxis):
            loc_min, loc_max = axis.get_view_interval()

            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                if tick.get_visible() and loc_min <= tick.get_loc() <= loc_max:
                    for artist in (tick.tick1line, tick.tick2line):
                        artist.draw(renderer)

        renderer.close_group(self.__name__)
        self.stale = False


def read_in_files(indir):

    # define dict for longform dataframe of main part of analysis
    data1 = {
            "cluster_type": [], # mmseq param ie 95_nucleotide
            "combination": [], # set of samples analyzed
            "category": [],
            "value": [],
            }

    # fractional data
    data2 = {
            "cluster_type": [], # mmseq param ie 95_nucleotide
            "combination": [], # set of samples analyzed
            "category": [],
            "value": [],
            }

    # define dict for longform dataframe of breakdown part of analysis
    data3 = {
            "cluster_type": [],
            "combination": [],
            "sample": [], # sample for flexible/specific gene breakdowns
            "category": [], # flexible or specific
            "count": []
            }

    # read through input files and populate dicts
    if indir[-1] == '/': indir = indir[:-1]
    file_list = [f for f in listdir(indir) if isfile(join(indir, f))]
    if '.DS_Store' in file_list: file_list.remove('.DS_Store')

    for file in file_list:
        cluster_type = file.split('_Summary.tsv')[0]
        with open(f'{indir}/{file}', 'r') as f:
            header = f.readline()
            for line in f:
                # split each line by tabs
                X = line.rstrip().split('\t')
                # define params
                combination = X[0]
                tgenes = int(X[1])
                sgenes = int(X[2])
                sfrac = float(X[3])
                fgenes = int(X[4])
                ffrac = float(X[5])
                pgenes = int(X[6])
                pfrac = float(X[7])
                samples = X[8]
                flex_break = X[9].split(', ')
                spec_break = X[10].split(', ')

                # Remove two categories showing reduntant info
                if combination in ['Station5A', 'Allminus']: continue

                # append count data to main dict
                count_data = {
                            "total_genes": tgenes,
                            "shared_genes": sgenes,
                            "flexible_genes": fgenes,
                            "specific_genes": pgenes
                            }
                for k, v in count_data.items():
                    data1["cluster_type"].append(cluster_type)
                    data1["combination"].append(combination)
                    data1["category"].append(k)
                    data1["value"].append(v)

                # append fractional data to 2nd dict
                frac_data = {
                            "shared_genes": sfrac,
                            "flexible_genes": ffrac,
                            "specific_genes": pfrac
                            }
                for k, v in frac_data.items():
                    data2["cluster_type"].append(cluster_type)
                    data2["combination"].append(combination)
                    data2["category"].append(k)
                    data2["value"].append(v)

                # append data to breakdon dict
                d = {'flex': flex_break, 'spec': spec_break}
                for name, brk in d.items():
                    for i in brk:
                        x = i.split(':')
                        sample = x[0]
                        count = x[1]

                        data3["cluster_type"].append(cluster_type)
                        data3["combination"].append(combination)
                        data3["sample"].append(sample)
                        data3["category"].append(name)
                        data3["count"].append(count)

    df1 = pd.DataFrame(data1)
    df2 = pd.DataFrame(data2)
    df3 = pd.DataFrame(data3)

    print('Main dataframe:\n', df1)
    print('\n\nFractional dataframe:\n', df2)
    print('\n\nBreakdown dataframe:\n', df3)

    return df1, df2, df3


def genetic_functional_diversity_plot(df, outpre, name):

    horder = ['95_nucleotide', '90_nucleotide', '70_aminoacid', '40_aminoacid']
    #colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3']
    colors = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c']

    g = sns.catplot(
                y='combination', x='value', data=df, hue='cluster_type',
                row='category', kind='bar', ci=None, sharex=False,
                hue_order=horder, palette=colors, orient='h',
                )

    for ax in g.axes.flat:
        ax.ticklabel_format(style='plain', axis='x')

    g.set_xticklabels(rotation=45)

    # adjust layout, save, and close
    g.tight_layout()
    g.savefig(f'{outpre}_{name}.pdf')
    plt.close()

    return True


def breakdown_plot(df, outpre):
    return True


def main():
    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_directory',
        help='Directory with *_Summary.tsv files',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_path_prefix',
        help='What would you like to use for the Output path/prefix?',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n\n')

    # define input params
    indir = args['input_directory']
    outpre = args['output_path_prefix']

    df1, df2, df3 = read_in_files(indir)

    _ = genetic_functional_diversity_plot(df1, outpre, "count_data")
    _ = genetic_functional_diversity_plot(df2, outpre, "fractional_data")
    #_ = breakdown_plot(df3, outpre)


    print('\n\nComplete success space cadet! The script without errors.\n\n')

if __name__ == "__main__":
    main()
