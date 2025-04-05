#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
merge all unrelated genomes from ancestral populations with E280As

Did not use --keep-allele-order, which would have preserved more matching SNPs
This is because the VCF allele orders are more similar than the minor 
allele assignments with unbalanced numbers of genomes per ancestry. Flag will
likely add another ~25k variants.
"""

import os
import numpy as np
import pandas as pd

REDO = True

def run_plink(cmd, output):
    if REDO or not os.path.isfile(output+'.log'):
        print (cmd)
        os.system(cmd)
    else:
        print (output+'.log already exists')
    return output

if __name__ ==  '__main__':
    DIR1 = '/scratch/eger/projects/hgdp_1kg'
    DIR2 = '/scratch/eger/projects/E280A_SNP_arrays/E280A_340'
    DIR3 = '/scratch/eger/projects/E280A_SNP_arrays/E280A_340/GENESIS_v1'
    OUTDIR = DIR2+'/merge_hgdp_1kg'
    if not os.path.isdir(OUTDIR):
        os.makedirs(OUTDIR)

    # two version of reference panel samples
    keep1 = os.path.join(DIR1, 'GENESIS_v2', 'hgdp_1kg.AFR_AMR_EUR_unrelated_sampleIDs.txt')
    keep2 = os.path.join(DIR1, 'GENESIS_v2', 'hgdp_1kg.AFR_AMR_EUR_unrelated_downsampled_sampleIDs.txt')
    
    sample_map = os.path.join(DIR1, 'GENESIS_v2', 'hgdp_1kg.AFR_AMR_EUR_unrelated_samplemap.txt')
    df = pd.read_csv(sample_map, sep='\t', header=None)
    
    # very related and distantly related
    rels = os.path.join(DIR3, 'E280A_340.hg38_pcair_rels.txt')
    unrels = os.path.join(DIR3, 'E280A_340.hg38_pcair_unrels.txt')
    
    geno1 = os.path.join(DIR1, 'hgdp_1kg.hg38')
    geno2 = os.path.join(DIR2, 'E280A_340.hg38') # unimputed

    # ref panel subset
    output = os.path.join(DIR1, 'hgdp_1kg.AFR_AMR_EUR_unrelated.auto.hg38')
    cmd = ' '.join(['/home/eger/software/PLINK1.9/plink',
                    '--bfile %s' % geno1,
                    '--keep %s' % keep1,
                    '--make-bed',
                    '--out %s' % output])
    geno1 = run_plink(cmd, output)
    
    ## extract the positions that are in both
    df2 = pd.read_csv(geno2+'.bim', sep='\s+', header=None,
                        names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
    df2['END'] = df2['BP']+1
    range = os.path.join(OUTDIR, 'E280A_340_positions.txt')
    df2[['CHR', 'BP', 'END', 'SNP']].to_csv(range, sep='\t', header=False, index=False)

    output = os.path.join(OUTDIR, 'hgdp_1kg.AFR_AMR_EUR_unrelated.auto.hg38.pos')
    cmd = ' '.join(['/home/eger/software/PLINK1.9/plink',
                    '--bfile %s' % geno1,
                    '--extract range %s' % range,
                    '--make-bed',
                    '--out %s' % output])
    geno1 = run_plink(cmd, output)

    # add the rsIDs
    df1 = pd.read_csv(geno1+'.bim', sep='\s+', header=None,
                        names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
    df1['ID'] = df1.apply(lambda row: 
                    str(row['CHR'])+':'+str(row['BP'])+':'+row['A2']+':'+row['A1'], axis=1)
    df2['ID'] = df2.apply(lambda row: 
                    str(row['CHR'])+':'+str(row['BP'])+':'+row['A2']+':'+row['A1'], axis=1)
    df1['SNP'] = df1['ID'].map(dict(zip(df2['ID'], df2['SNP']))).fillna('.')                      
    
    df1[['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2']].to_csv(
        geno1+'.bim' ,sep='\t', header=False, index=False)

    # 474725 variants
    extract = os.path.join(OUTDIR, 'hgdp_1kg_E280A_340_variants.txt')
    df1[df1['SNP'] != '.'][['SNP']].to_csv(extract, sep='\t', header=False, index=False) 

    # extract the SNPs from geno1
    output = os.path.join(OUTDIR, 'hgdp_1kg.AFR_AMR_EUR_unrelated.auto.hg38.ext')
    cmd = ' '.join(['/home/eger/software/PLINK1.9/plink',
                    '--bfile %s' % geno1,
                    '--extract %s' % extract,
                    '--make-bed',
                    '--out %s' % output])
    geno1a = run_plink(cmd, output)

    ## change the FID to anc in the .fam files
    df1 = pd.read_csv(geno1a+'.fam', sep='\s+', header=None,
                        names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'])
    df1['FID'] = df1['IID'].map(dict(zip(df[0], df[1])))
    df1.to_csv(geno1a+'.fam', sep=' ', index=False, header=False)
    
    # extract the SNPs from geno2
    output = os.path.join(OUTDIR, 'E280A_340.hg38.ext')
    cmd = ' '.join(['/home/eger/software/PLINK1.9/plink',
                    '--bfile %s' % geno2,
                    '--extract %s' % extract,
                    '--make-bed',
                    '--out %s' % output])
    geno2 = run_plink(cmd, output)

    ## merge the two datasets
    output = os.path.join(OUTDIR, 'hgdp_1kg.AFR_AMR_EUR_unrelated_E280A_340.auto.hg38')
    cmd = ' '.join(['/home/eger/software/PLINK1.9/plink',
                    '--bfile %s' % geno1a,
                    '--bmerge %s' % geno2,
                    '--make-bed',
                    '--out %s' % output])
    run_plink(cmd, output)
    # 474725 variants and 1138 people