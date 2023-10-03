# convert INTAA data into table with aa pair, score and count columns

#--------------------------------------------------------------------------------------

import pandas as pd
from collections import defaultdict
from Bio.PDB import *
import numpy as np

#--------------------------------------------------------------------------------------

df = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/rbi/intaa/rbi_intaa_combined.csv')
df.drop_duplicates(inplace = True)

df['r1'] = df['PDB_ID_x'].astype(str) + " " + df['trbi'].astype(str) + " " + df['prbi'].astype(str) + " " + df['Energy'].astype(str)
aa1 = list(df['r1'])
aa = pd.DataFrame(list(set(aa1)))
aa.columns = ['aa']

# Split a single column into multiple using str.split()
aa[['pdb', 'trbi', 'prbi', 'energy']] = aa.aa.str.split(expand = True)


def heatmap_df(aa):
    l, m, n = [], [], []
    c = len(aa)
    
    for j in range(0, c):
        for i in range(j, c):
            if aa[i] in Polypeptide.standard_aa_names and aa[j] in Polypeptide.standard_aa_names:
                m.append(aa[i])
                l.append(aa[j])
                n.append(str(aa[j]) + " " + str(aa[i]))
    
    data = [l, m, n]
    df = pd.DataFrame(data)
    df = df.transpose()
    df.columns = ['AA1', 'AA2', 'pair']
            
    return df

def aa_pair(df, aa_order):
    aa1 = df['trbi']
    aa2 = df['prbi']
    n = len(df)
            
    df['joined'] = aa1 + " " + aa2
    df['reverse'] = aa2 + " " + aa1
    
    aadf = heatmap_df(aa_order)
    pairlist = list(aadf['pair'])
    
    p_dict = defaultdict(list)
    # #initializing all key-value pairs
    # for i in range(0, len(pairlist)):
    #     p_dict[pairlist[i]] = []
        
    for i in range(len(df)):
        pair1 = df.iloc[i]['joined']
        try_pair = df.iloc[i]['reverse']
        if pair1 in pairlist:
            p_dict[pair1].append(df.iloc[i]['energy'])
        else:
            p_dict[pair1] = 'NA'
        # else:
        #     try_pair = df.iloc[i, -1]
        #     if try_pair in pairlist:
        #         if try_pair in p_dict:
        #             p_dict[try_pair]+=1
        #         else:
        #             p_dict[try_pair] = 1    
       
    return p_dict

aa_order = ['GLY', 'ALA', 'VAL', 'LEU', 'MET', 'ILE', 'PHE', 'TYR', 'TRP', 'SER', 'THR', 'CYS', 'PRO', 'ASN', 'GLN', 'LYS', 'ARG', 'HIS', 'ASP', 'GLU', 'MOD']
result = aa_pair(aa, aa_order)

c, n, d = 0, 0, 0
new_dict = defaultdict()
for i in list(result.keys()):
    c, n = 0, 0
    a = result[i]
    d += 1
    if a != 'NA':
        for idx, item in enumerate(a):
            if item != 'NA':
                c += float(item)
                n += 1
        if n != 0:
            new_dict[i] = c/n
        else:
            new_dict[i] = 'NA'
    else:
        new_dict[i] = 'NA'
# print(len(list(result.keys())))
pd.DataFrame.from_dict(new_dict, orient = 'index').to_csv('/home/user/revathy/bifur/programs/new_code/test_files/intaa_rbi.csv')
