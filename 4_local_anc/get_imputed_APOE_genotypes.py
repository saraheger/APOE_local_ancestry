#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# E280A snp array imputed
# get APOE genotypes

import os
import pandas as pd
import numpy as np

REDO = False

def run_plink(cmd, output):
    if REDO or not os.path.isfile(output+'.log'):
        print (cmd)
        os.system(cmd)
    else:
        print (output+'.log already exists')
    return output

if __name__ ==  '__main__':
    DIR1 = '/scratch/eger/projects/E280A_SNP_arrays'
    DIR2 = DIR1+'/APOE_genotypes'
    if not os.path.isdir(DIR2):
        os.makedirs(DIR2)

    cohort = 'Esc-Imputed-HiQual'
    geno = os.path.join(DIR1, cohort+'.hg38')

    # make extract file
    extract = ['rs429358', 'rs7412']
    df0 = pd.DataFrame()
    df0['SNP'] = extract
    fname_ext = os.path.join(DIR2, 'extract_APOE_genotype_SNPs.tsv')
    df0.to_csv(fname_ext, sep='\t', header=False, index=False)
    
    output = os.path.join(DIR2, 'APOE_genotype_SNPs')
    cmd = ' '.join(['plink',
                    '--bfile %s'% geno, # original plink files
                    '--extract %s' % fname_ext, # .tsv with SNP names
                    '--recode A', # this changes Cs and Ts to 0, 1, 2
                    '--out %s' % output,]) # APOE output file
    run_plink(cmd, output)

    df1 = pd.read_csv(output+'.raw', sep='\s+')
    df1['rs429358'] = df1['rs429358_C']
    df1['rs7412'] = df1['rs7412_T']
    
    # make subject ID the index
    df1.index = df1.IID
    # drop any subjects that don't have both SNPs
    df1 = df1[['rs429358', 'rs7412']].dropna()
    
    # Rules from https://www.snpedia.com/index.php/APOE
    # this finds the row index of subjects with each genotype
    # creates a dictionary where the key is the genotype 
    # and the value is the indexes of subjects with that genotype
    ids = {}
    ids['E1/E1'] = df1.loc[(df1.rs429358 == 2) & (df1.rs7412 == 2)].index
    ids['E1/E2'] = df1.loc[(df1.rs429358 == 1) & (df1.rs7412 == 2)].index
    # Note we assume E2/E4 is more likely than E1/E3
    ids['E2/E4'] = df1.loc[(df1.rs429358 == 1) & (df1.rs7412 == 1)].index
    ids['E1/E4'] = df1.loc[(df1.rs429358 == 2) & (df1.rs7412 == 1)].index
    ids['E2/E2'] = df1.loc[(df1.rs429358 == 0) & (df1.rs7412 == 2)].index
    ids['E2/E3'] = df1.loc[(df1.rs429358 == 0) & (df1.rs7412 == 1)].index
    ids['E3/E3'] = df1.loc[(df1.rs429358 == 0) & (df1.rs7412 == 0)].index
    ids['E3/E4'] = df1.loc[(df1.rs429358 == 1) & (df1.rs7412 == 0)].index
    ids['E4/E4'] = df1.loc[(df1.rs429358 == 2) & (df1.rs7412 == 0)].index
    
    # create a new dataframe
    df_APOE = pd.DataFrame()
    # list of subjects
    subjs = []
    # list of subjects' APOE genotypes
    apoes = []
    
    # loop through each genotype (key) in the dictionary
    for key in ids.keys():
        # prints how many subjects have the given genotype
        print ('APOE', key, len(ids[key]))
        # adds the subjects to the subject list
        subjs += list(ids[key])
        # adds their APOEs to the APOE list
        apoes += [key]*len(ids[key])
    # makes two columns in the new dataframe from the two lists
    df_APOE['SUBJID'] = subjs
    df_APOE['APOE'] = apoes
    # makes the index of the dataframe the SUBJIDs
    df_APOE.index = df_APOE.SUBJID
    
    # make a .tsv of the dataframe
    output = os.path.join(DIR2, 'APOE_genotypes_imputed.tsv')
    df_APOE.to_csv(output, header=False, index=False, sep='\t')
    print(output)
