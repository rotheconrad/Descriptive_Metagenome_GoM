#!/usr/bin/env python

''' Filters and renames contigs in fasta file

Removes contigs < user defined length in base pairs (bp) and renames
contigs in sequential order while appending the user provided prefix.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: February 7th, 2020
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


def Filter_contig_fasta(infile, outfile, prefix, minlen):

    X = infile.split('.')

    passed, failed = 0, 0

    with open(infile, 'r') as f, open(outfile, 'w') as o:

        for name, seq in read_fasta(f):
            
            contig_length = len(seq)

            if contig_length < minlen:
                failed += 1

            elif contig_length >= minlen:
                passed += 1
                line_out = f'>{prefix}_contig_{passed}\n{seq}\n'
                o.write(line_out)

            else:
                print('\n\nERROR: Something is wrong')

    print(
        f'\n\n\t\tContigs Passed: {passed}\n'
        f'\t\tContigs Removed: {failed}\n\n'
        )

    return True


def main():
    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--fasta_input_file',
        help='Please specify the fasta input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file',
        help='Please specify the fasta input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-p', '--contig_prefix',
        help='Please specify the prefix for contigs!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-m', '--minimum_contig_length',
        help='Please specify the minimum contig length in bp (default=1000)!',
        metavar=':',
        type=int,
        default=1000,
        required=False
        )
    args=vars(parser.parse_args())

    infile = args['fasta_input_file']
    outfile = args['output_file']
    prefix = args['contig_prefix']
    minlen = args['minimum_contig_length']

    _ = Filter_contig_fasta(infile, outfile, prefix, minlen)
    
if __name__ == "__main__":
    main()

