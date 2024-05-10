#!/usr/bin/env python

''' Append sample name to fasta sequences

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: July 9th, 2021
License :: GNU GPLv3
Copyright 2021 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, subprocess

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

def Fasta_rename_sequences(infile, sample_name):

    outfile = sample_name + '.rename'

    with open(infile, 'r+') as f, open(outfile, 'w') as o:

        for name, seq in read_fasta(f):
            newName = f'>{sample_name}_{name[1:]}\n{seq}\n'
            o.write(newName)

    _ = subprocess.run(['mv', outfile, infile])

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
    parser.add_argument(
        '-s', '--sample_name',
        help='Please specify the sample name to append to sequences!',
        metavar=':',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Run the renaming function:
    # Print renaming file
    Fasta_rename_sequences(
                    args['input_file'],
                    args['sample_name']
                    )
    
if __name__ == "__main__":
    main()

