#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ADMIXTURE prep
# match the SNPs
# LD prune

import os
import numpy as np
import pandas as pd

REDO = False

def run_plink(cmd, output):
    if REDO or not os.path.isfile(output+'.log'):
        print (cmd)
        os.system(cmd)
    else:
        print (output+'.log already exists')
    return output

if __name__ ==  '__main__':
    DIR0 = '/home/eger/projects/E280A'
    DIR1 = '/scratch/eger/projects/hgdp_1kg'
    DIR2 = '/scratch/eger/projects/E280A_SNP_arrays'
    DIR3 = DIR2+'/E280A_340/GENESIS_v1'
    DIR4 = DIR2+'/E280A_340/merge_hgdp_1kg'
    OUTDIR = DIR2+'/E280A_340/ADMIXTURE_v1'

    if not os.path.isdir(OUTDIR):
        os.makedirs(OUTDIR)
    
    # fam file from Hudson Alpha
    fam_340_file = os.path.join(DIR0, '2021.11.30-E280A-340_mac6.fam')
    
    # 309 reference samples for ADMIXTURE
    keep1 = os.path.join(DIR1, 'GENESIS_v2', 
            'hgdp_1kg.AFR_AMR_EUR_unrelated_downsampled_sampleIDs.txt')
    
    # 340 samples used in the paper
    keep2 = os.path.join(DIR2, 'E280A_340', 'E280A_340_VCF_IDs.txt')

    # ancestry assignments for reference samples
    sample_map = os.path.join(DIR1, 'GENESIS_v2', 
            'hgdp_1kg.AFR_AMR_EUR_unrelated_downsampled_samplemap.txt')
    
    # very related and distantly related
    rels = os.path.join(DIR3, 'E280A_340.hg38_pcair_rels.txt')
    unrels = os.path.join(DIR3, 'E280A_340.hg38_pcair_unrels.txt')

    # plink files from the VCFs
    geno1 = os.path.join(DIR1, 'hgdp_1kg.hg38')
    geno2 = os.path.join(DIR2, 'CLM-Esc-unimputed-HiQual.hg38') 

    # all unimputed positions
    range = os.path.join(DIR4, 'E280A_340_positions.txt')

    ## Extract unimputed positions & samples for analysis
    output = os.path.join(OUTDIR, 'hgdp_1kg.AFR_AMR_EUR_unrelated_downsampled.auto.hg38.pos')
    cmd = ' '.join(['/home/eger/software/PLINK1.9/plink',
                    '--bfile %s' % geno1,
                    '--keep %s' % keep1,
                    '--keep-allele-order',
                    '--extract range %s' % range,
                    '--make-bed',
                    '--out %s' % output])
    geno1 = run_plink(cmd, output)

    output = os.path.join(OUTDIR, 'E280A_340.auto.hg38')
    cmd = ' '.join(['/home/eger/software/PLINK1.9/plink',
                    '--bfile %s' % geno2,
                    '--keep %s' % keep2,
                    '--keep-allele-order',
                    '--make-bed',
                    '--out %s' % output])
    geno2 = run_plink(cmd, output)
    
    # add the rsIDs
    df1 = pd.read_csv(geno1+'.bim', sep='\s+', header=None,
                        names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
    df2 = pd.read_csv(geno2+'.bim', sep='\s+', header=None,
                        names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
    
    df1['ID'] = df1.apply(lambda row: 
                    str(row['CHR'])+':'+str(row['BP'])+':'+row['A2']+':'+row['A1'], axis=1)
    df2['ID'] = df2.apply(lambda row: 
                    str(row['CHR'])+':'+str(row['BP'])+':'+row['A2']+':'+row['A1'], axis=1)

    print(len(set(df1['ID']) & set(df2['ID']))) # 506694
    
    df1['SNP'] = df1['ID'].map(dict(zip(df2['ID'], df2['SNP']))).fillna('.')                      
    
    df1[['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2']].to_csv(
        geno1+'.bim' ,sep='\t', header=False, index=False)

    extract = os.path.join(OUTDIR, 'shared_variants.txt')
    df1[df1['SNP'] != '.'][['SNP']].to_csv(extract, sep='\t', header=False, index=False) 

    ## Extract the shared SNPs from the reference
    output = os.path.join(OUTDIR, 'hgdp_1kg.AFR_AMR_EUR_unrelated_downsampled.auto.hg38.ext')
    cmd = ' '.join(['/home/eger/software/PLINK1.9/plink',
                    '--bfile %s' % geno1,
                    '--extract %s' % extract,
                    '--keep-allele-order',
                    '--make-bed',
                    '--out %s' % output])
    geno1 = run_plink(cmd, output)

    ## LD pruning of reference samples
    prune = os.path.join(OUTDIR, 'hgdp_1kg.AFR_AMR_EUR_unrelated_downsampled.auto.hg38.ext.indep_SNPs')
    cmd = ' '.join(['/home/eger/software/PLINK1.9/plink',
                    '--bfile %s' % geno1,
                    '--indep-pairwise 50 10 0.1',
                    '--out %s' % prune])
    run_plink(cmd, prune)

    output = os.path.join(OUTDIR, 'hgdp_1kg.AFR_AMR_EUR_unrelated_downsampled.auto.hg38.ext.pruned')
    cmd = ' '.join(['/home/eger/software/PLINK1.9/plink',
                    '--bfile %s' % geno1,
                    '--extract %s' % prune+'.prune.in',
                    '--keep-allele-order',
                    '--make-bed',
                    '--out %s' % output])
    geno1 = run_plink(cmd, output) # 108115 variants
    
    ## Extract the same SNPs from the study samples
    output = os.path.join(OUTDIR, 'E280A_340.auto.hg38.ext.pruned')
    cmd = ' '.join(['/home/eger/software/PLINK1.9/plink',
                    '--bfile %s' % geno2,
                    '--extract %s' % prune+'.prune.in',
                    '--keep-allele-order',
                    '--make-bed',
                    '--out %s' % output])
    geno2 = run_plink(cmd, output)

    ## Update the .fam files
    df = pd.read_csv(sample_map, sep='\t', header=None)
    df1 = pd.read_csv(geno1+'.fam', sep='\s+', header=None,
                        names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'])
    df1['FID'] = df1['IID'].map(dict(zip(df[0], df[1])))
    df1.to_csv(geno1+'.fam', sep=' ', index=False, header=False)

    df = pd.read_csv(fam_340_file, sep='\s+', header=None,
                      names=['FID', 'IID', 'FA', 'MO', 'SEX', 'PHENO'])
    df2 = pd.read_csv(geno2+'.fam', sep='\s+', header=None,
                      names=['FID', 'IID', 'FA', 'MO', 'SEX', 'PHENO'])
    for col in ['FID', 'FA', 'MO', 'SEX', 'PHENO']:
        df2[col] = df2['IID'].map(dict(zip(df['IID'], df[col])))
    df2.to_csv(geno2+'.fam', sep=' ', index=False, header=False)

    ## Optional step: split study samples into 2 groups
    output = os.path.join(OUTDIR, 'E280A_rels.auto.hg38.ext.pruned')
    cmd = ' '.join(['/home/eger/software/PLINK1.9/plink',
                    '--bfile %s' % geno2,
                    '--keep %s' % rels,
                    '--keep-allele-order',
                    '--make-bed',
                    '--out %s' % output])
    run_plink(cmd, output)

    output = os.path.join(OUTDIR, 'E280A_unrels.auto.hg38.ext.pruned')
    cmd = ' '.join(['/home/eger/software/PLINK1.9/plink',
                    '--bfile %s' % geno2,
                    '--keep %s' % unrels,
                    '--keep-allele-order',
                    '--make-bed',
                    '--out %s' % output])
    run_plink(cmd, output)
    
