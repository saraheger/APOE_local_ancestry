#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# srun -p interactive --pty --time=24:00:00 --ntasks=2 --mem=20G bash

import os
import pandas as pd
import numpy as np

if __name__ == '__main__':
    MAINDIR = '/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25/LAI'
    DIR1 = MAINDIR+'/APOE'

    # AFR=0, AMR=1, EUR=2
    # dos0 = os.path.join(DIR1, 'redlat_tangl.APOE.anc0.dosage.txt')
    # dos1 = os.path.join(DIR1, 'redlat_tangl.APOE.anc1.dosage.txt')
    # dos2 = os.path.join(DIR1, 'redlat_tangl.APOE.anc2.dosage.txt')

    # hap0 = os.path.join(DIR1,'redlat_tangl.APOE.anc0.hapcount.txt')
    # hap1 = os.path.join(DIR1, 'redlat_tangl.APOE.anc1.hapcount.txt')
    # hap2 = os.path.join(DIR1, 'redlat_tangl.APOE.anc2.hapcount.txt')
    
    # dos0 = os.path.join(DIR1, 'redlat_tangl.AFR_AMR_EUR.APOE.anc0.dosage.txt')
    # dos1 = os.path.join(DIR1, 'redlat_tangl.AFR_AMR_EUR.APOE.anc1.dosage.txt')
    # dos2 = os.path.join(DIR1, 'redlat_tangl.AFR_AMR_EUR.APOE.anc2.dosage.txt')

    # hap0 = os.path.join(DIR1,'redlat_tangl.AFR_AMR_EUR.APOE.anc0.hapcount.txt')
    # hap1 = os.path.join(DIR1, 'redlat_tangl.AFR_AMR_EUR.APOE.anc1.hapcount.txt')
    # hap2 = os.path.join(DIR1, 'redlat_tangl.AFR_AMR_EUR.APOE.anc2.hapcount.txt')

    dos0 = os.path.join(DIR1, 'redlat_tangl.include_anc.APOE.anc0.dosage.txt')
    dos1 = os.path.join(DIR1, 'redlat_tangl.include_anc.APOE.anc1.dosage.txt')
    dos2 = os.path.join(DIR1, 'redlat_tangl.include_anc.APOE.anc2.dosage.txt')

    hap0 = os.path.join(DIR1,'redlat_tangl.include_anc.APOE.anc0.hapcount.txt')
    hap1 = os.path.join(DIR1, 'redlat_tangl.include_anc.APOE.anc1.hapcount.txt')
    hap2 = os.path.join(DIR1, 'redlat_tangl.include_anc.APOE.anc2.hapcount.txt')

    # read files
    df0d = pd.read_csv(dos0, sep='\t')
    df1d = pd.read_csv(dos1, sep='\t')
    df2d = pd.read_csv(dos2, sep='\t')

    df0h = pd.read_csv(hap0, sep='\t')
    df1h = pd.read_csv(hap1, sep='\t')
    df2h = pd.read_csv(hap2, sep='\t')

    # new dataframes
    df0 = df0d[['SAMPLE']].copy()
    df0['APOE'] = ''

    df1 = df0d[['SAMPLE']].copy()

    snps = ['chr19:44908684:T:C', 'chr19:44908822:C:T']

    for snp in snps:
        df0[snp+'_AFR_dos'] = df0d[snp]
        df0[snp+'_AMR_dos'] = df1d[snp]
        df0[snp+'_EUR_dos'] = df2d[snp]

        df0[snp+'_AFR_hap'] = df0h[snp]
        df0[snp+'_AMR_hap'] = df1h[snp]
        df0[snp+'_EUR_hap'] = df2h[snp]
    
    # sum dosages to get APOE genotypes
    df1['rs429358'] = df0[['chr19:44908684:T:C_AFR_dos', 'chr19:44908684:T:C_AMR_dos',
                        'chr19:44908684:T:C_EUR_dos']].sum(axis=1)
    df1['rs7412'] = df0[['chr19:44908822:C:T_AFR_dos', 'chr19:44908822:C:T_AMR_dos',
                        'chr19:44908822:C:T_EUR_dos']].sum(axis=1)

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

    # loop through each genotype (key) in the dictionary
    for key,vals in ids.items():
        # prints how many subjects have the given genotype
        print ('APOE', key, len(ids[key]))

        # add the genotypes to the dataframe
        df0.loc[vals, 'APOE'] = key

    # separate alleles
    df0['APOE_H1'] = df0['APOE'].apply(lambda x: x.split('/')[0])
    df0['APOE_H2'] = df0['APOE'].apply(lambda x: x.split('/')[1])

    # check that 2 SNPs have matching anc for each individual
    for anc in ['AFR', 'AMR', 'EUR']:
        df_hap = df0[df0['chr19:44908684:T:C_'+anc+'_hap'] != df0['chr19:44908822:C:T_'+anc+'_hap']]
        if not df_hap.empty:
            print(df_hap)

    # recode as anc labels
    df0['LAI'] = ''
    df0.loc[(df0['chr19:44908684:T:C_AFR_hap'] == 2), 'LAI'] = 'AFR/AFR'
    df0.loc[(df0['chr19:44908684:T:C_AFR_hap'] == 1) & (df0['chr19:44908684:T:C_AMR_hap'] == 1), 'LAI'] = 'AFR/AMR'
    df0.loc[(df0['chr19:44908684:T:C_AFR_hap'] == 1) & (df0['chr19:44908684:T:C_EUR_hap'] == 1), 'LAI'] = 'AFR/EUR'
    df0.loc[(df0['chr19:44908684:T:C_AMR_hap'] == 2), 'LAI'] = 'AMR/AMR'
    df0.loc[(df0['chr19:44908684:T:C_AMR_hap'] == 1) & (df0['chr19:44908684:T:C_EUR_hap'] == 1), 'LAI'] = 'AMR/EUR'
    df0.loc[(df0['chr19:44908684:T:C_EUR_hap'] == 2), 'LAI'] = 'EUR/EUR'

    # both haps are the same ancestry
    df0['LA_H1'] = np.where(df0['LAI'].isin(['EUR/EUR', 'AMR/AMR', 'AFR/AFR']), 
                            df0['LAI'].apply(lambda x: x.split('/')[0]), np.nan)
    df0['LA_H2'] = np.where(df0['LAI'].isin(['EUR/EUR', 'AMR/AMR', 'AFR/AFR']), 
                            df0['LAI'].apply(lambda x: x.split('/')[-1]), np.nan)

    # both alleles are the same
    df0['LA_H1'] = np.where(df0['APOE'].isin(['E3/E3', 'E4/E4']), 
                            df0['LAI'].apply(lambda x: x.split('/')[0]), df0['LA_H1'])
    df0['LA_H2'] = np.where(df0['APOE'].isin(['E3/E3', 'E4/E4']), 
                            df0['LAI'].apply(lambda x: x.split('/')[-1]), df0['LA_H2'])

    # different alleles and different ancestries
    # ['E2/E3', 'E3/E4', 'E2/E4']
    for i, row in df0[df0['LA_H1'].isna()].iterrows():
        ancestries = row['LAI'].split('/')
        
        if row['APOE'] == 'E2/E3':
            # E2 anc
            if (row['chr19:44908822:C:T_AFR_dos'] == 1):
                df0.loc[i, 'LA_H1'] = 'AFR'
            elif (row['chr19:44908822:C:T_AMR_dos'] == 1):
                df0.loc[i, 'LA_H1'] = 'AMR'
            else:
                df0.loc[i, 'LA_H1'] = 'EUR'
            # E3 anc
            if df0.loc[i, 'LA_H1'] == ancestries[0]:
                df0.loc[i, 'LA_H2'] = ancestries[1]
            else:
                df0.loc[i, 'LA_H2'] = ancestries[0]
        
        else: #E2/E4 or E3/E4
            # E4 anc
            if (row['chr19:44908684:T:C_AFR_dos'] == 1):
                df0.loc[i, 'LA_H2'] = 'AFR'
            elif (row['chr19:44908684:T:C_AMR_dos'] == 1):
                df0.loc[i, 'LA_H2'] = 'AMR'
            else:
                df0.loc[i, 'LA_H2'] = 'EUR'
            # E2 or E3 anc
            if df0.loc[i, 'LA_H2'] == ancestries[0]:
                df0.loc[i, 'LA_H1'] = ancestries[1]
            else:
                df0.loc[i, 'LA_H1'] = ancestries[0]

    df0 = df0[['SAMPLE', 'APOE', 'APOE_H1', 'LA_H1', 'APOE_H2', 'LA_H2',
        'chr19:44908684:T:C_AFR_dos',
       'chr19:44908684:T:C_AMR_dos', 'chr19:44908684:T:C_EUR_dos',
       'chr19:44908684:T:C_AFR_hap', 'chr19:44908684:T:C_AMR_hap',
       'chr19:44908684:T:C_EUR_hap', 'chr19:44908822:C:T_AFR_dos',
       'chr19:44908822:C:T_AMR_dos', 'chr19:44908822:C:T_EUR_dos',
       'chr19:44908822:C:T_AFR_hap', 'chr19:44908822:C:T_AMR_hap',
       'chr19:44908822:C:T_EUR_hap']]

    # out1 = os.path.join(DIR1, 'redlat_tangl.APOE.LAI.genotypes.tsv')
    # out1 = os.path.join(DIR1, 'redlat_tangl.AFR_AMR_EUR.APOE.LAI.genotypes.tsv')
    out1 = os.path.join(DIR1, 'redlat_tangl.include_anc.APOE.LAI.genotypes.tsv')
    df0.to_csv(out1, sep='\t', index=False)
    print(out1)

    pd.pivot_table(df0[['APOE_H1', 'LA_H1']], 
                    index='APOE_H1', columns='LA_H1', 
                    aggfunc=len, fill_value=0, margins=True)

    pd.pivot_table(df0[['APOE_H2', 'LA_H2']], 
                    index='APOE_H2', columns='LA_H2', 
                    aggfunc=len, fill_value=0, margins=True)