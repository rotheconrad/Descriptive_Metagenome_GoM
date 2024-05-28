#!/usr/bin/env python

''' Subsample a paired read interleaved fasta file

Reads the fasta file into memory and writes output files in interleaved
fasta format for each percent to subsample value requested.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: March 31st, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, random
from collections import defaultdict


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


def interleaved_toDict(fasta):

    data = defaultdict(list)
    readCount = 1
    pairCount = 1
    readPair = f'readPair_{pairCount}'

    # Read the fasta file into a dictionary of readPair = [read1, read2]
    with open(fasta, 'r') as file:
        for name, seq in read_fasta(file):

            if readCount % 2 == 0:

                entry = f'>{readPair}_read2\n{seq}\n'
                data[readPair].append(entry)

                pairCount += 1
                readPair = f'readPair_{pairCount}'

                if pairCount % 1000000 == 0:
                    print(f'\t... Read pairs processed: {pairCount}')

            else:

                entry = f'>{readPair}_read1\n{seq}\n'
                data[readPair].append(entry)

            readCount += 1

    totalReadPairs = len(data)

    print(f'\n\tNumber of read pairs processed: {totalReadPairs}')

    return data, totalReadPairs


def subsample_fasta(data, totalReadPairs, subsample, prefix):

    # for percent submitted, select random subsample and write to file
    for percent in subsample:
        print(f'\n\tSubsampling at {percent}% ...')
        subsample_percent = int(totalReadPairs * percent / 100)
        random.seed()
        subsampleDict = dict(random.sample(data.items(), subsample_percent))

        outFile = f'{prefix}_sbsmpl{percent}.fa'

        print(
            f'\tNumber of read pairs selected for {percent}% subsample: '
            f'{len(subsampleDict)}\n\tWritting reads to file: {outFile}\n'
            )

        with open(outFile, 'w') as fout:
            for readPair, entries in subsampleDict.items():
                for entry in entries:
                    fout.write(entry)

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-s', '--percents_to_subsample',
        help=
            'List of whole integer values as percents to subsample '
            '(ie: -s 20 40 60 80',
        metavar='',
        type=int,
        nargs='+',
        required=True
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the paired read interleaved fasta file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_prefix',
        help='Please specify the prefix to use for the output files!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    print('\n\nRunning Script ...')
    print('\n\nReading paired read interleaved fasta file ...\n')
    data, totalReadPairs = interleaved_toDict(args['input_file'])
    print('\n\nSubsampling data and writing output files ...')
    _ = subsample_fasta(
                    data,
                    totalReadPairs,
                    args['percents_to_subsample'],
                    args['output_prefix']
                    )


if __name__ == "__main__":
    main()