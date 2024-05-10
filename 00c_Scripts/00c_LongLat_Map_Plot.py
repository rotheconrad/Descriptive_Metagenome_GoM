#!/usr/bin/env python

'''Plot Longitude and Latitude points on a world map.

This script uses GeoPandas: https://geopandas.org/

# Installation:
> conda install -c plotly plotly_express
> conda install -c conda-forge python-kaleido

Input:
    - Metadata tsv file with 3 columns: Sample, Long, Lat.
    - Bounding Box boundaries for the map. 
      i.e. min and max x and y axis values and long and lat.

Output:
    - Plot as PDF.

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
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def LonLat_map_plot(infile, outfile, xmin, xmax, ymin, ymax):

    # read the tsv file to pandas dataframe
    df = pd.read_csv(infile, sep='\t')
    print('\n\n', df, '\n\n')

    # Plot the points
    x = 1
    point_size = [x, x, x]
    fig = px.scatter_geo(
                    df, lat=df.Lat, lon=df.Lon,
                    text=df.Station, size=point_size,
                    )

    # define the bounding box if params provided
    if xmin and xmax and ymin and ymax:
        fig.add_traces(
            go.Scattergeo(lon=[xmin, xmax], lat=[ymin, ymax],
            mode='markers', marker=dict(size=2, color='rgba(0,0,0,0)'),
            opacity=0
            ))

    # add some plot modifications
    fig.update_geos(
                visible=True, scope='north america', resolution=50,
                fitbounds='locations',
                showcountries=True, countrycolor="Black",
                showsubunits=True, subunitcolor="grey",
                showcoastlines=True, showlakes=True,
                showrivers=True,
                )

    # save the plot to file
    fig.write_image(outfile)
    
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
        '-o', '--output_file',
        help='Please specify the output file name!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-xmin', '--longitude_minimum',
        help='OPTIONAL: x-axis minimum!',
        metavar='',
        type=float,
        required=False
        )
    parser.add_argument(
        '-xmax', '--longitude_maximum',
        help='OPTIONAL: x-axis maximum!',
        metavar='',
        type=float,
        required=False
        )
    parser.add_argument(
        '-ymin', '--latitude_minimum',
        help='OPTIONAL: x-axis minimum!',
        metavar='',
        type=float,
        required=False
        )
    parser.add_argument(
        '-ymax', '--latitude_maximum',
        help='OPTIONAL: x-axis maximum!',
        metavar='',
        type=float,
        required=False
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...')

    # define parameters
    infile = args["input_file"]
    outfile = args["output_file"]
    xmin = args["longitude_minimum"]
    xmax = args["longitude_maximum"]
    ymin = args["latitude_minimum"]
    ymax = args["latitude_maximum"]

    # Run function
    _ = LonLat_map_plot(infile, outfile, xmin, xmax, ymin, ymax)


if __name__ == "__main__":
    main()
