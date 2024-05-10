#!/usr/bin/env python

'''Build bar plots for sequencing and assembly data

Takes TSV file with following columns:
    * the order and columns are important.

    1) Sample Name # these are x-axis labels
    2) Layer # this is sample description / environment etc.
    3) Sequencing Effort (Gbp) # Giga base pairs
    4) Metagenome Reads GC (%)
    5) Average Genome Size (Mbp) # Mega base pairs
    6) NonPareil Diversity
    7) Contigs Assembled
    8) Assembled Length (Mbp) # Mega base pairs
    9) Assembled N50
    10) Assembled GC (%)
    11) Predicted CDS
    12) Reads Mapping to Assembly (%)
    13) Reads Mapping to rMAGs (%)

Outputs bar plots in PDF ready for publication.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: March 2022
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


def Sequencing_Assembly_Plots(infile, outpre):

    # read the tsv file to pandas dataframe
    df = pd.read_csv(infile, sep='\t', index_col=0)
    print('\n\n', df, '\n\n')

    order = [
        'EN21', 'EN51', 'EN81', 'EN22', 'EN52', 'EN82',
        'EN23', 'EN53', 'EN83', 'EN24', 'EN54', 'EN84',
        'EN25', 'EN55', 'EN85', 'EN56', 'EN57', 'EN58', 'EN59'
        ]

    df = df.loc[order]

    columns = df.columns.values.tolist()[1:]

    for i, col in enumerate(columns):
        print(col)
        fig, ax = plt.subplots(figsize=(3,6), dpi=300)

        df.plot.barh(
                y=col, ax=ax, alpha=0.75,
                color="#525252", width=.98, legend=False
                )
        ax.set_ylabel("Sample Name", fontsize=12)
        ax.set_xlabel(col, fontsize=12)
        #ax.minorticks_on()
        ax.tick_params(axis='y', which='minor', left=False)
        ax.tick_params(labelsize=10, direction='inout', width=2, length=6)
        # set grid style
        ax.xaxis.grid(
            which="minor", color='#737373', linestyle='--', linewidth=1
            )
        ax.xaxis.grid(
            which="major", color='#737373', linestyle='--', linewidth=1
            )
        ax.invert_yaxis()
        ax.add_artist(TickRedrawer())
        ax.set_axisbelow(True)

        # adjust layout, save, and close
        fig.set_tight_layout(True)
        plt.savefig(f'{outpre}_{i:.2f}.pdf')
        plt.close()


    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the input file name!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_prefix',
        help='Please specify the prefix to use for the output file name!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...')

    # define parameters
    infile = args["input_file"]
    outpre = args["output_file_prefix"]

    # Run function
    _ = Sequencing_Assembly_Plots(infile, outpre)


if __name__ == "__main__":
    main()
