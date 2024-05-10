#!/usr/bin/env python

'''Builds Temp-Salinity plots with density contour.

This version of the script takes two input files.

1) The running CTD data for the stations.
    "Station_Depth_Temp_Salinity.tsv"
2) The data just for the metagenome sampling points. These are plotted
    on top of the running CTD data in (1)
    "Sample_Temp_Salinity_Density.tsv"

Need to install python seawater or gsw package to calculate density:
https://pypi.org/project/seawater/
conda install -c conda-forge seawater or pip install seawater

########################################################################
# Seawater alternative that I found but didn't use.
https://pypi.org/project/gsw/
conda install -c conda-forge gsw or pip install gsw
########################################################################

Colors points by station and highlights the metagenome sample points.

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
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seawater
import random


def TS_plot(Stations, Samples, Outfile):
    "reads file builds plots"

    # Load Data
    stations_df = pd.read_csv(Stations, sep='\t')
    samples_df = pd.read_csv(Samples, sep='\t')

    ####################################################################
    ## Sample 2 has several less salty points in the CTD cast and the
    ## metagenome from sample 2 surface point is low as well.
    ## I made two plots and overlaid them in illustrator because it was
    ## faster than fidgeting with plot overlay params in matplotlib.

    ## Turn this bit on to remove low salinity points for the main plot
    stations_df = stations_df[stations_df["Salinity"] >= 34]
    samples_df = samples_df[samples_df["Salinity"] >= 34]

    ## Turn this bit on to plot the low salinity points
    #stations_df = stations_df[stations_df["Salinity"] < 34]
    #samples_df = samples_df[samples_df["Salinity"] < 34]
    ####################################################################
    # Load the metadata markers and colors
    station_color = "#bdbdbd"
    markers = {8: "s", 5: "o" , 2: "D"} # for stations
    colors = { # for metagenome samples
                "SURFACE": "#F5C714",
                "ML": "#5DC7D3", 
                "DCM": "#70BF52", 
                "AOM": "#E0546C", 
                "OM": "#B24198", 
                "DEEP": "#3D94D1"
            }

    #### Setup Density contours ######################
    temp = stations_df["Temperature"].to_numpy()
    salt = stations_df["Salinity"].to_numpy()

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
    
    station_labels = stations_df["Station"].unique()
    layer_labels = samples_df["Layer"].unique()

    for X in station_labels:
        dfx = stations_df[stations_df["Station"] == X]
        x = dfx["Salinity"].to_numpy()
        y = dfx["Temperature"].to_numpy()
        ax.plot(
            x, y, linestyle='None', marker=markers[X], markersize=5,
            color=station_color, label=X, alpha=0.25
            )

    for X in layer_labels:
        dfx = samples_df[samples_df["Layer"] == X]

        for Y in station_labels:
            dfy = dfx[dfx["Station"] == Y]

            x = dfy["Salinity"].to_numpy()
            y = dfy["Temp"].to_numpy()
            ax.plot(
                x, y, linestyle='None', marker=markers[Y], markersize=10,
                color=colors[X], label=X, alpha=0.5)

    ax.set_xlabel('Salinity (PSU)', fontsize=12, y=-0.02)
    ax.set_ylabel('Temperature (Celsius)', fontsize=12, x=-0.02)

    #h, l = ax.get_legend_handles_labels()
    #ordered_handles = [h[3], h[2], h[1], h[0], h[4]]
    #ordered_labels = [l[3], l[2], l[1], l[0], l[4]]

    ax.legend(
        #ordered_handles, ordered_labels,
        bbox_to_anchor=(1.2, 0.5), loc="center", ncol=1, frameon=False,
        markerscale=1, fontsize=12
        )
    fig.set_figwidth(8)
    fig.set_figheight(6)
    #plt.subplots_adjust(left=0.15, bottom=0.1, right=0.995, top=0.995)
    plt.tight_layout()
    plt.savefig(Outfile)
    plt.close() 

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-f1', '--Station_SalTemp_file',
        help='Please specify the input file name!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-f2', '--Sample_SalTemp_file',
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
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...')

    # define params

    Stations = args["Station_SalTemp_file"]
    Samples = args["Sample_SalTemp_file"]
    Outfile = args["output_file"]
    _ = TS_plot(Stations, Samples, Outfile)


if __name__ == "__main__":
    main()
