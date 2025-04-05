#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# conda activate rfmix
# 442164

# genetic map came from this post
# https://www.biostars.org/p/9517179/
# https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/

# https://console.cloud.google.com/storage/browser/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2

# using the downsampled reference panel

import os
import numpy as np
import pandas as pd
from multiprocessing import Pool

def run(cmd):
    print(cmd)
    os.system(cmd)

if __name__ ==  '__main__':
    tasks = 6
    
    RESOURCE_DIR = '/home/eger/references/gnomad/resources/1kg/1000G_genetic_map'
    DIR1 = '/scratch/eger/projects/hgdp_1kg'
    DIR2 = '/scratch/eger/projects/E280A_SNP_arrays'
    OUTDIR1 = DIR2+'/RFmix_LAI_v2'
    OUTDIR2 = DIR2+'/RFmix_LAI_v2/Tractor'
    for dir in [OUTDIR1, OUTDIR2]:
        if not os.path.isdir(dir):
            os.makedirs(dir)
    os.chdir(OUTDIR1)

    # files
    genetic_map = os.path.join(RESOURCE_DIR, 'genetic_map_hg38_withX.txt.gz')
    sample_map = os.path.join(DIR1, 'GENESIS_v2',
                'hgdp_1kg.AFR_AMR_EUR_unrelated_downsampled_samplemap.txt')
    ref_vcf = os.path.join(DIR1, 'hgdp1kgp_autosome.filtered.SNV_INDEL.phased.shapeit5.vcf.gz')
    query_bcf = os.path.join(DIR2, 'Esc-Imputed-HiQual_R2-0.3.bcf')
    query_vcf = os.path.join(DIR2, 'Esc-Imputed-HiQual_R2-0.3.vcf.gz')

    # edit the genetic map file
    rfmix_genetic_map = os.path.join(RESOURCE_DIR, 'genetic_map_hg38_withX_map_only.txt')
    if not os.path.exists(rfmix_genetic_map):
        df2 = pd.read_csv(genetic_map, sep='\s+')
        df2['chr'] = df2['chr'].apply(lambda x: 'chrX' if x == 23 else 'chr'+str(x))
        df2[['chr', 'position', 'Genetic_Map(cM)']].to_csv(
            rfmix_genetic_map, header=False, index=False, sep='\t')

    # run RFmix
    cmds = []
    for ch in ['chr'+str(i) for i in np.arange(1, 23)]:
        ref_bcf = os.path.join('/home/eger/references/gnomad/resources/hgdp_1kg',
                               'hgdp1kgp_'+ch+'.filtered.SNV_INDEL.phased.shapeit5.bcf')
        out = 'E280A_LAI.'+ch
        if not os.path.exists(OUTDIR1+'/'+out+'.msp.tsv'):
            cmd = ' '.join(['rfmix',
                            '-f %s' % query_bcf,
                            '-r %s' % ref_bcf,
                            '-m %s' % sample_map,
                            '-g %s' % rfmix_genetic_map,
                            '-o %s' % out,
                            '--chromosome=%s' % ch]) 
            cmds.append(cmd)

    pool = Pool(processes = tasks)
    pool.map(run, cmds)
    pool.close()
    pool.join()

    os.chdir(OUTDIR2)
    
    # soft link phased VCF by chr
    for ch in ['chr'+str(i) for i in np.arange(1, 23)]:
        chr_vcf = os.path.join(DIR2, 'RFmix_LAI_v1/Tractor', 'Esc-Imputed-HiQual_R2-0.3.'+ch+'.vcf.gz')
        if os.path.exists(chr_vcf):
            os.symlink(chr_vcf, os.path.join(OUTDIR2, 'Esc-Imputed-HiQual_R2-0.3.'+ch+'.vcf.gz'))

    # run Tractor
    cmds = []
    for ch in ['chr'+str(i) for i in np.arange(1, 23)]:
        msp = os.path.join(OUTDIR1, 'E280A_LAI.'+ch+'.msp.tsv')
        chr_vcf = os.path.join(OUTDIR2, 'Esc-Imputed-HiQual_R2-0.3.'+ch+'.vcf.gz')
        if not os.path.exists(os.path.join(OUTDIR2, 'Esc-Imputed-HiQual_R2-0.3.'+ch+'.anc0.dosage.txt')):
            cmd = ' '.join(['python',
                            '/home/eger/software/Tractor/scripts/extract_tracts.py',
                            '--msp %s' % msp,
                            '--vcf %s' % chr_vcf,
                            '--num-ancs 3',
                            '--output-dir %s' % OUTDIR2]) 
            cmds.append(cmd)

    pool = Pool(processes = tasks)
    pool.map(run, cmds)
    pool.close()
    pool.join()
