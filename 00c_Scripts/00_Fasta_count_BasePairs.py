#!/usr/bin/env python

''' Count total base pairs in a fasta file.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: June 21st, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
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


def Fasta_count_bps(infile):

    bp_count = 0

    with open(infile, 'r+') as f:
        for name, seq in read_fasta(f):
            bp_count += len(seq)

    print(f'{infile}\t{bp_count}')

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the fasta file to rename deflines!',
        metavar=':',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Count base pairs in file
    bp_count = Fasta_count_bps(args['input_file'])


if __name__ == "__main__":
    main()

