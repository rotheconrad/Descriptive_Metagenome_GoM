#!/usr/bin/env python

''' Metagenome Functional Abundance Analysis

This script combines the following four outputs:
    
    1) MMSeqs gene clustering, and representative genes
    2) Output files from MMSeqs cluster analysis custom python scripts
    3) Microbe Annotator MMSeqs cluster representative gene annotations
    4) Coverage Magic Gene TADs

The goal is to analyze functional differences between metagenomes.

This script will look at functional differences in three ways:

    1) Count the number of genes annotated as a particular function or
        functional category. Counts number of genes in each mmseqs cluster
        and combines clusters with same function/functional category.
        (optional) accepts list of user defined genes/functions.
    ?? 2) KEGG Metabolic completeness comparison of Microbe Annotatr output ??
    3) Coverage/Abundance of particular functions or functional categories.
        Normalizes coverage by GEQs from microbe census.

This script will output several tsv files and plots:

    1) list of output files

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
import Functional_Abundance_Functions as faf

# Set the default sans-serif font to Helvetica
# matplotlib.rcParams['font.sans-serif'] = "Helvetica"
# Set default font to sans-serif
# matplotlib.rcParams['font.family'] = "sans-serif"
# If you don't have Helvetica and need it see this:
# https://fowlerlab.org/2019/01/03/changing-the-sans-serif-font-to-helvetica/
# If you don't need Helvetica font just delete the above code.

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-m1', '--mmseqs_cluster_tsv',
        help='MMSeqs2 cluster file in tsv format!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-m2', '--mmseq_analysis_core_csv',
        help='MMSeqs2 cluster file in tsv format!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-a', '--MicrobeAnnotator_Annot',
        help='Microbe Annotator .annot file for rep genes!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-c', '--CoverageMagic_GeneTAD_Dir',
        help='Directory of CoverageMAgic Results with GeneTAD for each sample!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-mc', '--MicrobeCensus_Files_Dir',
        help='Directory of MircrobeCensus results for each sample!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-gc', '--Gene_Categories',
        help='(Optional) File with list of one gene category keyword per line!',
        metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-o', '--prefix_for_output_files',
        help='Prefix to use for the output files!',
        metavar='',
        type=str,
        required=True
        )

    args=vars(parser.parse_args())

    # Define input variable
    mmseq_cluster = args['mmseqs_cluster_tsv']
    mmseq_core = args['mmseq_analysis_core_csv']
    microbe_annotations = args['MicrobeAnnotator_Annot']
    coverage_magic = args['CoverageMagic_GeneTAD_Dir']
    microbe_census_dir = args['MicrobeCensus_Files_Dir']
    gene_categories_list = args['Gene_Categories']
    outpre = args['prefix_for_output_files']

    # Do what you came here to do:
    print('\n\nRunning Script ...\n')

    ###########################################################################
    ###### STEP 0 - Collect Data ##############################################
    ###########################################################################

    print('\t### STEP 0 ##########\n')

    # Parse MMSeqs Cluster file into two dictionaries. byCluster and bySample.
    print('\t\tParsing MMSeqs2 cluster tsv file ...\n')
    byCluster, bySample = faf.parse_mmseqs_cluster_tsv(mmseq_cluster)

    # Parse MMSeqs cluster analysis core csv file into dictionary
    print('\t\tParsing MMSeqs2 cluster analysis core csv file ...\n')
    coreClusters = faf.parse_mmseq_core_csv(mmseq_core)

    # Parse Microbe Annotator .annot file into dictionary
    print('\t\tParsing Microbe Annotator .annot file ...\n')
    repGeneAnnot = faf.parse_MicrobeAnnotator_annots(microbe_annotations)

    # Use microbe census to normalize gene coverage by genome equivalents
    microbecensus = microbe_census_dir
    print('\t\tParsing Microbe Census files ...\n')
    geqs = faf.parse_microbe_census(microbecensus)

    # Parse coverage magic geneTAD files into dictionary
    print('\t\tParsing CoverageMagic geneTAD files ...\n')
    CovMag = faf.parse_CoverageMagic_geneTAD(coverage_magic, geqs)

    geneCatList = args['Gene_Categories']
    geneCats = None
    if geneCatList:
        print('\t\tParsing Microbe Census files ...\n')
        geneCats = faf.parse_gene_categories(geneCatList)

    ###########################################################################
    ###### STEP 1 - Functional Counts #########################################
    ###########################################################################

    print('\t### STEP 1 ##########\n')


    # Count the gene functions in each sample by gene name, ID, or category
    print('\t\tCounting Functions in each sample ...\n')
    (
    geneNameCounts,
    clstrNameCounts,
    geneIDCounts,
    geneCatCounts,
    geneAbundance
                    ) = faf.get_funcational_abundance(
                                                    repGeneAnnot,
                                                    bySample,
                                                    byCluster,
                                                    geneCats,
                                                    outpre,
                                                    CovMag
                                                    )

    ###########################################################################
    ###### STEP 2 - Analyze and Plot ##########################################
    ###########################################################################

    print('\t### STEP 2 ##########\n')

    # do some kind of analysis and write output
    print('\t\tAnalyzing functional counts ...\n')
    #df = faf.analyze_functional_counts(geneNameCounts, geneIDCounts, geneCatCounts)

    # write some kind of plot
    print('\t\tPlotting functional analysis summary ...\n')
    #_ = faf.plot_functional_summary(df)
    # Normalize by shared gene average coverage
    # Find which shared gene clusters are most stable across samples
    # Look at relative abundance. Number of genes in cluster.
    # Look at difference between relative abundance within each cluster.

    print('\n\nComplete success space cadet! Finished without errors.\n\n')


if __name__ == "__main__":
    main()
