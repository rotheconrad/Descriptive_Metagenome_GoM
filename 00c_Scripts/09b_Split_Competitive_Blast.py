#!/usr/bin/env python

'''Split up competitive blast results

Uses second column of blast results to select a unique MAG identifier.
Uses underscore "_" as a delimeter and selects first two positions.
    Ex:
       if the sample name is EN_21_trim_45_scaffold_2506,
       then EN_21_trim_45 is the unique MAG identifier.

This can be changed on line

Input:
    A blast file from competitive mapping

Output:
    A new blast file for each unique MAG identitier

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
from collections import defaultdict

def split_competitive_blast(infile, outpre):

    # initial dict of lists to store results per unique MAG ID
    data = defaultdict(list)

    # read through file, retrieve unique MAG ID, and populate data dict
    with open(infile, 'r') as file:
        for line in file:
            # skip header lines
            if line.startswith('#'): continue
            # retrieve unique MAG ID
            MAG = '_'.join(line.split('\t')[1].split('_')[:5])
            # add line to data dict
            data[MAG].append(line)

    # write output file for each unique MAG ID
    for MAG, DATA in data.items():
        with open(f'{outpre}-{MAG}.blst', 'w') as out:
            for line in DATA:
                out.write(line)

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the tabular magic blast file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_prefix',
        help='Please specify the prefix for the output files!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')
    
    _ = split_competitive_blast(args['input_file'], args['output_prefix'])

    print('\n\nComplete success space cadet! Hold on to your boots.\n\n')


if __name__ == "__main__":
    main()
