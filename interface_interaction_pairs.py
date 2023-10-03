from Bio.PDB import *
import pandas as pd
import os, time
from collections import defaultdict
from alive_progress import alive_bar

#------------------------------------------------------------------------------------------

aa_order = ['GLY', 'ALA', 'VAL', 'LEU', 'MET', 'ILE', 'PHE', 'TYR', 'TRP', 'SER', 'THR', 'CYS', 'PRO', 'ASN', 'GLN', 'LYS', 'ARG', 'HIS', 'ASP', 'GLU']

hs = defaultdict()

def heatmap_df(aa):
    
    l, m, n = [], [], []
    c = len(aa)
    
    for j in range(0, c):
        for i in range(0, c):
            if aa[i] in Polypeptide.standard_aa_names and aa[j] in Polypeptide.standard_aa_names:
                m.append(aa[i])
                l.append(aa[j])
                n.append(str(aa[j]) + " " + str(aa[i]))
    
    data = [l, m, n]
    df = pd.DataFrame(data)
    df = df.transpose()
    df.columns = ['AA1', 'AA2', 'pair']
            
    return df

def aa_pair(list1, list2, aa_order):
    
    aa1 = [i[2:5] for i in list1]
    aa2 = [i[2:5] for i in list2]
            
    joined = [i + " " + j for i, j in zip(aa1, aa2)]
    reverse = [j + " " + i for i, j in zip(aa1, aa2)]

    aadf = heatmap_df(aa_order)
    pairlist = list(aadf['pair'])
    
    p_dict = defaultdict()
    # initializing all key-value pairs
    for i in range(0, len(pairlist)):
        p_dict[pairlist[i]] = 0
        
    for i in range(len(joined)):
        pair1 = joined[i]
        try_pair = reverse[i]
        if pair1 in pairlist:
            if pair1 in p_dict:
                p_dict[pair1] += 1
                p_dict[try_pair] += 1
            else:
                p_dict[pair1] = 1
                p_dict[try_pair] = 1
        else:
            try_pair = reverse[i]
            if try_pair in pairlist:
                if try_pair in p_dict:
                    p_dict[try_pair] += 1
                    p_dict[pair1] += 1
                else:
                    p_dict[try_pair] = 1    
                    p_dict[pair1] = 1
    
    return p_dict

    
os.chdir("/home/user/revathy/bifur/programs/new_code/test_files/interface/Interface_interactions/")
with alive_bar(682, bar='solid', theme = 'smooth') as bar:
    for file in os.listdir("/home/user/revathy/bifur/programs/new_code/test_files/interface/Interface_interactions/"):
        filename = os.fsdecode(file)                        # precautionary : decoding filenames if encoded (like in Windows)
        if filename.endswith(".csv"):
            df = pd.read_csv(filename)

            res1 = list(set(list(df[df.columns[1:3]].apply(lambda x: ' '.join(x.astype(str)), axis=1))))
            res2 = list(set(list(df[df.columns[4:6]].apply(lambda x: ' '.join(x.astype(str)), axis=1))))

            pair_true_aa = aa_pair(res1, res2, aa_order)
            
            for i in list(pair_true_aa.keys()):
                if i in list(hs.keys()):
                    hs[i] = int(hs[i]) + int(pair_true_aa[i])
                    
                else:
                    hs[i] = int(pair_true_aa[i])
        bar()
            
pd.DataFrame.from_dict(hs, orient = 'index').to_csv("/home/user/revathy/bifur/programs/new_code/test_files/final/interface_interaction_pairs.csv")