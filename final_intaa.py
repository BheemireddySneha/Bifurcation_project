# convert Consurf data into table with aa, score and count columns
#--------------------------------------------------------------------------------------
import pandas as pd
from collections import defaultdict
from Bio.PDB import *
import numpy as np
#--------------------------------------------------------------------------------------

# for RBIs
df = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/Analysis/INTAA_final_rbi_csv.csv')
# df = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/Analysis/INTAA_inter_final.csv')
df.drop(['Unnamed: 0'], axis = 1, inplace = True)

df['r1'] = df['PDB_ID'].astype(str) + " " + df['res_name1'].astype(str) + " " + df['res_name2'].astype(str) + " " + df['Energy'].astype(str)
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

    pair = [i + " " + j for i, j in zip(list(df['trbi']), list(df['prbi']))]
    try_pair = [j + " " + i for i, j in zip(list(df['trbi']), list(df['prbi']))]

    aadf = heatmap_df(aa_order)
    pairlist = list(aadf['pair'])
    
    p_dict = {}
    # initializing all key-value pairs
    for i in range(0, len(pairlist)):
        p_dict[pairlist[i]] = []
        
    for i in range(len(pair)):
        pair1 = pair[i]
        pair2 = try_pair[i]
        if pair1 in pairlist:
            if pair1 in p_dict:
                p_dict[pair1].append(df.iloc[i]['energy'])
            elif pair2 in p_dict:
                p_dict[pair2].append(df.iloc[i]['energy'])
            else:
                p_dict[pair1] = list(df.iloc[i]['energy'])
                p_dict[pair2] = list(df.iloc[i]['energy'])
    # aa1 = df['trbi']
    # aa2 = df['prbi']
            
    # df['joined'] = aa1 + " " + aa2
    # df['reverse'] = aa2 + " " + aa1
    
    # aadf = heatmap_df(aa_order)
    # pairlist = list(aadf['pair'])

    # p_dict = defaultdict(list)
    # #initializing all key-value pairs
    # for i in range(0, len(pairlist)):
    #     p_dict[pairlist[i]] = []
        
    # for i in range(len(df)):
    #     pair1 = df.iloc[i]['joined']
    #     try_pair = df.iloc[i]['reverse']
    #     if pair1 in pairlist:
    #         if pair1 in p_dict:
    #             p_dict[pair1].append(df.iloc[i]['energy'])
    #             # p_dict[try_pair].append(df.iloc[i]['energy'])
    #         else:
    #             p_dict[pair1] = list(df.iloc[i]['energy'])
    #             # p_dict[try_pair] = list(df.iloc[i]['energy'])

    #     else:
    #         if try_pair in pairlist:
    #             if try_pair in p_dict:
    #                 p_dict[try_pair].append(df.iloc[i]['energy'])
    #                 # p_dict[pair1].append(df.iloc[i]['energy'])
    #             else:
    #                 p_dict[try_pair] = list(df.iloc[i]['energy'])
    #                 # p_dict[pair1] = list(df.iloc[i]['energy'])
       
    return p_dict

aa_order = ['GLY', 'ALA', 'VAL', 'LEU', 'MET', 'ILE', 'PHE', 'TYR', 'TRP', 'SER', 'THR', 'CYS', 'PRO', 'ASN', 'GLN', 'LYS', 'ARG', 'HIS', 'ASP', 'GLU']
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
                try:
                    c += float(item)
                    n += 1
                except:
                    continue
            else:
                c = 0
                n = 0
        if n != 0:
            new_dict[i] = c/n
        else:
            new_dict[i] = 0
    else:
        new_dict[i] = 0

pd.DataFrame.from_dict(new_dict, orient = 'index').to_csv('/home/user/revathy/bifur/programs/new_code/test_files/final/intaa_rbi3.csv')
# pd.DataFrame.from_dict(new_dict, orient = 'index').to_csv('/home/user/revathy/bifur/programs/new_code/test_files/final/intaa_inter.csv')