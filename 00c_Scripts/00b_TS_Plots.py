#!/usr/bin/env python

'''Builds Temp-Salinity plots with density contour.

Input is Metagenome_Samples_MetaData.tsv from GoM project.

Need to install python seawater or gsw package to calculate density:
https://pypi.org/project/seawater/
conda install -c conda-forge seawater or pip install seawater

https://pypi.org/project/gsw/
conda install -c conda-forge gsw or pip install gsw

Colors points by Ocean and labels points by depth.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: May 13 2021
License :: GNU GPLv3
Copyright 2021 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seawater
import random

def rnd():
    x = random.uniform(-.2, .2)
    return x

def TS_plot(infile, outfile, metacolors):
    "reads file builds plots"

    # Load Metadata
    df = pd.read_csv(infile, sep='\t')

    # Load the metadata colors
    metac = pd.read_csv(metacolors, sep='\t', header=0, dtype = str)
    colors = dict(zip(metac['Labels'], metac['Colors']))

    # Adjustments for text labels.
    posx = {
                'Gulf Of Mexico': -0.0,
                'South Pacific': -0.0,
                'North Pacific': -0.0,
                'South Atlantic': 0.0,
                'North Atlantic': 0.0
                }
    posy = {
                'Gulf Of Mexico': -0.0,
                'South Pacific': 0.0,
                'North Pacific': 0.0,
                'South Atlantic': -0.0,
                'North Atlantic': -0.0
                }

    #### Setup Density contours ######################
    temp = df.Temp.to_numpy()
    salt = df.Salinity.to_numpy()

    # Figure out boudaries (mins and maxs)
    smin = salt.min() - (0.01 * salt.min())
    smax = salt.max() + (0.01 * salt.max())
    tmin = temp.min() - (0.1 * temp.max())
    tmax = temp.max() + (0.1 * temp.max())

    # Calculate how many gridcells we need in the x and y dimensions
    xdim = int(round((smax-smin)/0.1+1,0))
    ydim = int(round((tmax-tmin)+1,0))

    # Create empty grid of zeros
    dens = np.zeros((ydim,xdim))
     
    # Create temp and salt vectors of appropiate dimensions
    ti = np.linspace(1,ydim-1,ydim)+tmin
    si = np.linspace(1,xdim-1,xdim)*0.1+smin
     
    # Loop to fill in grid with densities
    for j in range(0,int(ydim)):
        for i in range(0, int(xdim)):
            dens[j,i] = seawater.eos80.dens0(si[i],ti[j])
     
    # Substract 1000 to convert to sigma-t
    dens = dens - 1000
    ###################################################

    # Plot data ***********************************************
    fig, ax = plt.subplots()
    CS = plt.contour(si,ti,dens, linestyles='dashed', colors='k')
    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.0f') # Label every second level
    
    for X in df.Ocean.unique():
        dfx = df[df.Ocean == X]
        x = dfx.Salinity.to_numpy()
        y = dfx.Temp.to_numpy()
        ax.plot(
            x, y, linestyle='None', marker="o", markersize=15,
            color=colors[X], label=X, alpha=0.5
            )

        '''
        for i, txt in enumerate(dfx.Depth):
            ax.annotate(
                f'{int(txt)}m', (x[i]+rnd(), y[i]+rnd()),
                fontsize=20
                )
        '''
     
    ax.set_xlabel('Salinity (ppt)', fontsize=30, fontweight='bold', y=-0.02)
    ax.set_ylabel('Temperature (C)', fontsize=30, fontweight='bold', x=-0.02)

    h, l = ax.get_legend_handles_labels()
    ordered_handles = [h[3], h[2], h[1], h[0], h[4]]
    ordered_labels = [l[3], l[2], l[1], l[0], l[4]]
    ax.legend(
        ordered_handles, ordered_labels,
        bbox_to_anchor=(0.5,1.02), loc="lower center", ncol=5, frameon=False,
        markerscale=2, fontsize='xx-large'
        )
    fig.set_figwidth(16)
    fig.set_figheight(6)
    #plt.subplots_adjust(left=0.15, bottom=0.1, right=0.995, top=0.995)
    plt.tight_layout()
    plt.savefig(f'{outfile}_TS_Plot.pdf', dpi=300)
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
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file',
        help='Please specify the output file name!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-c', '--metacolors',
        help='OPTIONAL: Metadata to color samples by!',
        metavar='',
        type=str,
        required=False
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...')
    TS_plot(args['input_file'], args['output_file'], args['metacolors'])


if __name__ == "__main__":
    main()
