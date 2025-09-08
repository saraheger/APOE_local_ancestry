#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd

if __name__ ==  '__main__':
    include_anc = ['AMR', 'CAH', 'EUR', 'AAC', 'AFR']
    exclude_anc = ['AJ', 'MDE', 'CAS', 'EAS']

    # missingness > 5%
    miss_ids = ['HALB3079914',
                'HALB3079915',
                'HALB3079917',
                'HALB3079916',
                'HALB3079918',
                'HALB3076440',
                'HALB3015251']

    DIR1 = '/cluster/projects/ADFTD/batch_calls/redlat_plus_tangl_10-25'
    DIR2 = '/cluster/home/seger/projects/psomagen3'
    DIR3 = '/cluster/home/seger/projects/redlat_tangl/redlat_plus_tangl_5-9-25/ADMIXTURE'

    fname0 = os.path.join(DIR1, 'tangl_ids.txt')
    fname1 = os.path.join(DIR1, 'redlat_ids.txt')
    fname2 = os.path.join(DIR2, 'AN00020190.ids')
    fname3 = os.path.join(DIR3, 'redlat_genotools_combined_ages.csv')
    
    df0 = pd.read_csv(fname0, header=None, names=['IID'])
    df1 = pd.read_csv(fname1, header=None, names=['IID'])
    df2 = pd.read_csv(fname2, sep='\t')
    df2.columns = ['PSOMAGENID', 'IID']
    df3 = pd.read_csv(fname3)

    # redlat missing anc
    df4 = pd.DataFrame({'IID': sorted(set(df1['IID']) - set(df3['IID']) - set(df2['IID']))})
    
    out1 = os.path.join(DIR3, 'redlat_ids_missing_genotools_anc.txt')
    df4.to_csv(out1, header=False, index=False)
    print(out1)

    # tangl + psomagen3 + redlat include_anc individuals
    redlat_include_anc = set(df3[df3['GenoTools_Label'].isin(include_anc)]['IID']) # 2626
    df5 = pd.DataFrame({'IID': sorted(set(df0['IID']) | set(df2['IID']) | redlat_include_anc)}) # 3566
    df5['FID'] = 0
    df5 = df5[~df5['IID'].isin(miss_ids)] # 3565

    out2 = os.path.join(DIR3, 'include_anc_ids.fam')
    df5[['FID', 'IID']].to_csv(out2, sep='\t', header=False, index=False)
    print(out2)