#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# create an reference panel for RFmixv2

# list of AMR populations to use for reference panel
# https://github.com/broadinstitute/gnomad_local_ancestry/issues/29
# 1kg: PEL
# HGDP: Colombian, Karitiana, Maya, Pima, Surui
# Final reference is 241 samples: gs://gnomad-mwilson/lai/sample_selection/ref_to_run.tsv

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
    cohort = 'hgdp_1kg'
    DIR1 = '/home/eger/references/gnomad/resources/hgdp_1kg'
    OUTDIR1 = '/scratch/eger/projects/'+cohort
    if not os.path.isdir(OUTDIR1):
        os.makedirs(OUTDIR1)

    anc_dict = {'Central South Asia (HGDP)': 'SAS',
               'Oceania (SGDP),Oceania (HGDP)': 'OCE',
               'Middle East (HGDP)': 'MDE',
               'Oceania (HGDP)': 'OCE',
               'Europe (HGDP)': 'EUR',
               'Africa (HGDP)': 'AFR',
               'East Asia (HGDP)': 'EAS',
               'America (HGDP)': 'AMR'}
    
    #### ALL SAMPLES ###
    # list of BCFs for a single VCF
    file_list = os.path.join(OUTDIR1, 'hgdp_1kg_chr_file.list')
    fnames = [DIR1+'/hgdp1kgp_chr'+str(i)+'.filtered.SNV_INDEL.phased.shapeit5.bcf' for i in np.arange(1, 23)]
    pd.DataFrame(fnames).to_csv(file_list, header=False, index=False)

    # concatenate BCFs
    vcf = os.path.join(OUTDIR1, 'hgdp1kgp_autosome.filtered.SNV_INDEL.phased.shapeit5.vcf.gz')
    cmd = 'bcftools concat -f '+file_list+' -O z -o '+vcf
    print(cmd)
    os.system(cmd)

    # index the VCF
    cmd = 'bcftools index '+vcf
    print(cmd)
    os.system(cmd)
    
    # convert to plink
    output = os.path.join(OUTDIR1, cohort+'.hg38')
    cmd = ' '.join(['/home/eger/software/PLINK1.9/plink',
                    '--vcf %s'% vcf,
                    '--const-fid', 
                    '--keep-allele-order',
                    '--allow-no-sex',
                    '--make-bed',
                    '--out %s' % output])
    geno = run_plink(cmd, output)
    # 76,413,521 variants

    # hwe, maf, geno
    output = os.path.join(OUTDIR1, cohort+'.auto.hg38.geno001.maf005.hwe10-5')
    cmd = ' '.join(['/home/eger/software/PLINK1.9/plink',
                    '--bfile %s'% geno,
                    '--geno 0.01',
                    '--maf 0.05',
                    '--hwe 0.00001',
                    '--allow-no-sex',
                    '--make-bed',
                    '--out %s' % output])
    geno = run_plink(cmd, output)
    # 3,531,419 variants
