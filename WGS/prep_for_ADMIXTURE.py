#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# panel created on server
# /home/eger/references/gnomad/resources/hgdp_1kg/ADMIXTURE_99AFR_99AMR_99EUR

# srun -p interactive --pty --time=24:00:00 --mem=50G --cpus-per-task 10 --ntasks 2 bash
# ml load cluster/plink/1.90

# ADMIXTURE prep
# match the SNPs
# LD prune

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
    DIR1 = '/cluster/home/seger/references/gnomad/resources/hgdp_1kg/ADMIXTURE_99AFR_99AMR_99EUR'
    DIR2 = '/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25/global_ancestry/ADMIXTURE'

    # plink files from the VCFs
    geno1 = os.path.join(DIR1, 'hgdp_1kg.99AFR_99AMR_99EUR.auto.hg38.maf010')
    geno2 = os.path.join(DIR2, 'redlat_tangl.auto.qceed.include_anc') 

    # read in the bims
    df1 = pd.read_csv(geno1+'.bim', sep='\s+', header=None,
                        names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])


    extract = os.path.join(DIR2, 'ref_variants.txt')
    df1[['SNP']].to_csv(extract, header=False, index=False)

    ## Extract the same variants from the study samples
    output = os.path.join(DIR2, 'redlat_tangl.auto.qceed.include_anc.ref_variants')
    cmd = ' '.join(['plink',
                    '--bfile %s' % geno2,
                    '--extract %s' % extract,
                    '--keep-allele-order',
                    '--make-bed',
                    '--out %s' % output])
    geno2 = run_plink(cmd, output)

    df2 = pd.read_csv(geno2+'.bim', sep='\s+', header=None,
                        names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])

    # shared variants
    extract = os.path.join(DIR2, 'shared_variants.txt')
    df2[['SNP']].to_csv(extract,  header=False, index=False)

    ## Extract from ref
    output = os.path.join(DIR2, 'hgdp_1kg.99AFR_99AMR_99EUR.auto.hg38.maf010.shared')
    cmd = ' '.join(['plink',
                    '--bfile %s' % geno1,
                    '--keep-allele-order',
                    '--extract %s' % extract,
                    '--make-bed',
                    '--out %s' % output])
    geno1 = run_plink(cmd, output)

    ## LD pruning of reference samples
    prune = os.path.join(DIR2, 'hgdp_1kg.99AFR_99AMR_99EUR.auto.hg38.maf010.shared.indep_SNPs')
    cmd = ' '.join(['plink',
                    '--bfile %s' % geno1,
                    '--indep-pairwise 50 10 0.1',
                    '--out %s' % prune])
    run_plink(cmd, prune)

    output = os.path.join(DIR2, 
        'hgdp_1kg.99AFR_99AMR_99EUR.auto.hg38.maf010.shared.pruned')
    cmd = ' '.join(['plink',
                    '--bfile %s' % geno1,
                    '--extract %s' % prune+'.prune.in',
                    '--keep-allele-order',
                    '--make-bed',
                    '--out %s' % output])
    geno1 = run_plink(cmd, output) # 125039 variants

    ## Extract the same variants from the study samples
    output = os.path.join(DIR2, 'redlat_tangl.auto.qceed.include_anc.pruned')
    cmd = ' '.join(['plink',
                    '--bfile %s' % geno2,
                    '--extract %s' % prune+'.prune.in',
                    '--keep-allele-order',
                    '--make-bed',
                    '--out %s' % output])
    geno2 = run_plink(cmd, output)

