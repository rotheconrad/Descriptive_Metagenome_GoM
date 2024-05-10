#!/usr/bin/env python

'''Calculates GC for each contig and plots GC distribution.

Optional: Provide contig name to highlight in the distribution plot.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: September 08, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def get_GC(string):
    seq = string.upper()
    C = seq.count('C')
    G = seq.count('G')
    GC = (G + C) / len(seq) * 100
    return GC


def get_average_GC(infile):

    gc_dist = []

    with open(infile, 'r') as f:
        for name, seq in read_fasta(f):
            GC = get_GC(seq)
            gc_dist.append(GC)

    average_GC = sum(gc_dist) / len(gc_dist)

    print(infile, average_GC)

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify an input fasta file!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    infile = args['input_file']

    _ = get_average_GC(infile)


if __name__ == "__main__":
    main()


