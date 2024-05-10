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
        Normalizes by the average coverage of genes shared across all
        metagenome samples in the set.

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

from os import listdir
from os.path import isfile, join
from collections import defaultdict
import pandas as pd

###############################################################################
########## STEP 0 - Collect Data ##############################################
###############################################################################

def parse_mmseqs_cluster_tsv(infile):

    '''
    The mmseqs cluster tsv file is two columns.
    column 1 is the cluster representative sequence.
    This sequence is repeated for every sequence in the cluster.
    Column 2 are sequences in the cluster.
    The general idea with this function is to create 2 dicionaries.
    One dictionary stores data by the cluster represenative name as in
    {cluster: sample}. Which samples are in each cluster.
    The other dicitonary stores data by the sample name.
    {sample: cluster}. Which clusters are in each sample.
    '''

    # initialize variables
    # Store data in two dictionaries.
    byCluster = defaultdict(dict)
    bySample = defaultdict(dict)

    # read through the file and populate the dictionaries
    with open(infile, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            cluster = X[0]
            gene = X[1]
            sample = 'EN' + gene.split('_')[1]
            # store data in dicts
            byCluster[cluster][gene] = ''
            bySample[sample][cluster] = ''

    return byCluster, bySample

###############################################################################

def parse_mmseq_core_csv(infile):

    '''
    This file comes from the custom mmseqs analysis script.
    Core gene csv files are output for different sample combinations
    The input file here is the core gene csv for whichever combination
    '''

    # initialize dict to store {core_gene: ''}
    # we'll gene abundance info downstream
    coreClusters = {}

    # read through input file and population dict
    with open(infile, 'r') as file:
        for line in file:
            gene = line.rstrip().split(', ')[1]
            coreClusters[gene] = ''

    return coreClusters

###############################################################################

def parse_MicrobeAnnotator_annots(infile):

    '''
    This file is a tsv with several columns of relavant information
    This function will place info in a dict with gene names as the key
    The gene names are the representative genes from mmseq clustering.
    '''

    #initialize dict to store gene annot info
    repGeneAnnot = {}
    databases = defaultdict(int)
    # read through input file and populate dict
    with open(infile, 'r') as file:
        header = file.readline()
        for line in file:
            X = line.rstrip().split('\t')

            repGene = X[0]
            protID = X[1]
            KEGG = X[3]

            if KEGG != 'NA': geneID, geneName = X[3], X[4]
            elif protID != 'NA': geneID, geneName = X[1], X[2]
            else: print(X)

            repGeneAnnot[repGene] = [geneID, geneName]

    return repGeneAnnot

###############################################################################

def parse_microbe_census(indir):

    '''
    Get genome equivalents for each sample
    '''

    # initialize dict to store coverage results
    geqs = {}

    # read through input files and populate dict
    if indir[-1] == '/': indir = indir[:-1]
    file_list = [f for f in listdir(indir) if isfile(join(indir, f))]
    if '.DS_Store' in file_list: file_list.remove('.DS_Store')

    for file in file_list:
        with open(f'{indir}/{file}', 'r') as f:
            smpl = 'EN' + file.split('.')[0].split('_')[1]
            for line in f:
                X = line.rstrip().split('\t')
                if X[0] == 'genome_equivalents:':
                    geq = float(X[1])
            geqs[smpl] = geq

    return geqs

###############################################################################

def parse_CoverageMagic_geneTAD(indir, geqs):

    '''
    These are 3 column tsv files for each sample.
    columns are gene name, tad, genelength
    '''

    # initialize dict to store coverage results
    CovMag = defaultdict(lambda: defaultdict(list))

    # read through input files and populate dict
    if indir[-1] == '/': indir = indir[:-1]
    file_list = [f for f in listdir(indir) if isfile(join(indir, f))]
    if '.DS_Store' in file_list: file_list.remove('.DS_Store')

    for file in file_list:
        sample = 'EN' + file.split('_')[1]
        print(f'\t\t\t\tSample {sample} ...')
        geq = geqs[sample]
        with open(f'{indir}/{file}', 'r') as f:
            header = f.readline()
            for line in f:
                X = line.rstrip().split('\t')
                gene = X[0]
                tad = float(X[1]) / geq
                length = int(X[2])
                CovMag[sample][gene] = [tad, length]

    return CovMag

###############################################################################

def parse_gene_categories(geneCatList):

    '''
    Read in gene category list to python list
    '''

    geneCats = []

    with open(geneCatList, 'r') as file:
        for line in file:
            geneCats.append(line.rstrip())

    return geneCats

###############################################################################
########## STEP 1 - Functional Counts #########################################
###############################################################################

def compile_dfs(gene_list, data_dict, samples, outpre, label):

    data = {'Genes': gene_list}

    for samp in samples:
        samp_list = []
        dd = data_dict[samp]
        for gene in gene_list:
            samp_list.append(dd[gene])
        data[samp] = samp_list

    df = pd.DataFrame(data)
    df = df.set_index('Genes')

    df.to_csv(f'{outpre}_{label}.tsv', sep='\t')

    return df

###############################################################################

def get_funcational_abundance(
                    repGeneAnnot, bySample, byCluster, geneCats, outpre, CovMag
                    ):

    '''
    Gets functional abundance in a few ways:

    1) Counts number of genes with same gene name in each sample
    2) Increases count in (1) by number of genes within each cluster
    3) Sum cluster TAD80 for each gene name (function) in each sample
    4) Counts number of genes with same gene ID in each sample
    5) Longform tsv file with abundance, sample, and annotation of each gene

    Writes tsv file of counts with genes-rows, samples-columns
    '''

    # Counts for each cluster based on gene annotation name
    geneNameCounts = defaultdict(lambda: defaultdict(int))
    # increase geneNameCounts by all genes within each cluster for each sample
    clstrNameCounts = defaultdict(lambda: defaultdict(int))
    # abundance for each cluster based on TAD80 / Genome Equivalents
    geneAbundance = defaultdict(lambda: defaultdict(int))
    # counts by geneID (KEGG KO, SwissProt, or RefSeq)
    geneIDCounts = defaultdict(lambda: defaultdict(int))
    # Counts by user defined categories
    geneCatCounts = defaultdict(lambda: defaultdict(int))
    # Long form datatable of with tad80/geq.
    longform = {'Gene': [], 'Sample': [], 'Annotation': [], 'Abundance': []}

    names_d = {}
    ids_d = {}
    samples = []

    # if the sample has a gene in the repgene cluster count the annotation
    # else set the count to zero for that annotation in that sample.
    for sample, genes in bySample.items():
        print(f'\t\t\t\tSample {sample} ...')
        samples.append(sample)
        for gene, name in repGeneAnnot.items():
            geneID = name[0]
            geneName = name[1].split(' [')[0]
            names_d[geneName] = ''
            ids_d[geneID] = ''

            if gene in genes:
                geneNameCounts[sample][geneName] += 1
                geneIDCounts[sample][geneID] += 1

                # increase geneNameCounts by all genes within the cluster
                # Sum gene abundance for each cluster
                s_count = 0
                clstr_abudnance = 0

                for g in byCluster[gene]:
                    s = 'EN' + g.split('_')[1]
                    if s == sample:
                        s_count += 1
                        tad80 = CovMag[sample][g][0]
                        clstr_abudnance += tad80

                        #update longform
                        longform['Gene'].append(g) 
                        longform['Sample'].append(sample)
                        longform['Annotation'].append(geneName)
                        longform['Abundance'].append(tad80)

                clstrNameCounts[sample][geneName] += s_count
                geneAbundance[sample][geneName] += clstr_abudnance

            else:
                if geneName not in geneNameCounts[sample]:
                    geneNameCounts[sample][geneName] = 0
                    clstrNameCounts[sample][geneName] = 0
                if geneName not in geneAbundance[sample]:
                    geneAbundance[sample][geneName] = 0
                if geneID not in geneIDCounts[sample]:
                    geneIDCounts[sample][geneID] = 0

    print('\t\tCompiling gene name counts ...\n')
    names = list(names_d.keys())
    name_df = compile_dfs(names, geneNameCounts, samples, outpre, 'GeneNames')
    print('\t\tCompiling gene name counts plus clustered gene counts ...\n')
    clstr_df = compile_dfs(names, clstrNameCounts, samples, outpre, 'ClstrNames')
    print('\t\tCompiling genee abundance ...\n')
    abndnc_df = compile_dfs(names, geneAbundance, samples, outpre, 'geneAbunance')
    print('\t\tCompiling gene ID counts ...\n')
    ids = list(ids_d.keys())
    id_df = compile_dfs(ids, geneIDCounts, samples, outpre, 'GeneIDS')
    if geneCats:
        print('\t\tCompiling gene category counts ...\n')
        cat_df = compile_dfs(geneCats, geneCatCounts, samples, outpre, 'GeneCats')
    else:
        cat_df = None

    # write out longform
    long_df = pd.DataFrame(longform)
    long_df.to_csv(f'{outpre}_longform.tsv', sep='\t', index=False)

    return name_df, clstr_df, id_df, cat_df, abndnc_df

###############################################################################
########## STEP 2 - Exploratory Plots #########################################
###############################################################################

def plot_functional_summary(df):

    '''
    Counts number of genes with same gene ID in each sample
    Writes tsv file of counts with genes-rows, samples-columns
    Look for significant differences and plot them
    '''

    return True

###############################################################################

def plot_functional_summary(df):

    '''
    Counts number of genes with same gene ID in each sample
    Writes tsv file of counts with genes-rows, samples-columns
    Look for significant differences and plot them
    '''

    return True

###############################################################################


