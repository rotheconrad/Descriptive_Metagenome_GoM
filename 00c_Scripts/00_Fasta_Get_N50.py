#!/usr/bin/env python

'''Returns N50 of a fasta file. For Metagenome assembly or MAG.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: January 06, 2020
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


def get_N50(length_list):
    half = sum(length_list) / 2
    descend_list = sorted(length_list, reverse=True)
    n = 0
    for i in descend_list:
        n += i
        if n >= half:
            return i


def return_N50(infile):

    length_list = []

    with open(infile, 'r') as f:
        for name, seq in read_fasta(f):
            length_list.append(len(seq))

    N50 = get_N50(length_list)

    print(infile, N50)

    return True

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--in_file',
        help='Please specify an input file!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    infile = args["in_file"]

    _ = return_N50(infile)


if __name__ == "__main__":
    main()


