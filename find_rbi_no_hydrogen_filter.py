#!/usr/bin/python3

#-------------------------------------------------------------------------------
import pandas as pd
from collections import defaultdict
from optparse import OptionParser
import time

start_time = time.time()
#-------------------------------------------------------------------------------
# Refine listed interactions with more stringent criteria: remove glycines, check for beta-carbons, etc.
def refine_interactions(df):

    # ensuring that no glycine residues are counted as tRBIs
    df = df[df['trbi'] != 'GLY']

    # ensuring that no CB atom is involved in any interaction other than a hydrophobic interaction
    df.reset_index(drop = True, inplace = True)

    for index, row in df.iterrows():
        try:
            if 'hydrophobic interaction' not in row['inttype']:
                if row['aa1'] == 'CB' or row['aa2'] == 'CB':
                    df.drop(df.index[index], inplace = True)
        except:
            continue     

    # further refining interaction types depending on the interacting atoms                     
    df.reset_index(drop = True, inplace = True)
    for i, row in df.iterrows():
        if 'disulphide bond' in df['inttype'][i]:
            if df['trbi'][i] == 'MET' or df['prbi'][i] == 'MET':
                df.drop([i], inplace = True)
                i = i+1
                continue

    return df

#-------------------------------------------------------------------------------
# Removing all main_chain_atom-main_chain_atom interactions
def remove_mainchain_interactions(input_df):

    atomA = input_df['Atom_type1'].tolist()
    atomB = input_df['Atom_type2'].tolist()
    main_chain_atom = ['O', 'C', 'CA', 'N']
    c = len(input_df)

    n_dict = defaultdict(list)
    n_dict_keys = input_df.columns.tolist()

    for row in range(0, c):
        if atomA[row] in main_chain_atom and atomB[row] in main_chain_atom:
            continue
        else:
            for col in n_dict_keys:
                r = input_df.iloc[row][col]
                n_dict[col].append(r)

    df2 = pd.DataFrame(n_dict)
    df2.drop_duplicates(subset=["Chain1", "AA1Res_num1", "Atom_type1", "Chain2", "AA2Res_num2", "Atom_type2"], keep="first", inplace = True)

    return df2

#-------------------------------------------------------------------------------

usage = "USAGE: python3 find_rbi_no_hydrogen_filter.py --f1 file"
parser = OptionParser(usage=usage)
parser.add_option("--f1", help="file", dest="f1")
(options, args) = parser.parse_args()

df = pd.read_csv(options.f1)

# remove mc-mc interactions
df2 = remove_mainchain_interactions(df)

#-------------------------------------------------------------------------------
'''
defining true RBIs:
    1. participate in > 1 main-chain-side-chain or side-chain-side-chain interaction
    2. have a side-chain participating component
'''
#-------------------------------------------------------------------------------
def true_rbi(df):

    c = len(list(df.index.values)) - 1

    res = defaultdict()
    df_dict = defaultdict()

    # include all interface residues in a dictionary {residue: [interacting residues]}
    while c > -1:
        
        c1 = df.iloc[c]["Chain1"] + " " + df.iloc[c]["AA1Res_num1"]
        c2 = df.iloc[c]["Chain2"] + " " + df.iloc[c]["AA2Res_num2"]

        if c1  not in res:
            res[c1] = [c2]
        else:
            res[c1].append(c2)
        
        if c2 not in res:
            res[c2] = [c1]
        else:
            res[c2].append(c1) 

        c -= 1

    for i in res.keys():
        s = list(set(res[i]))
        res[i] = s

    # ordering the dictionary alphabetically
    res = dict(sorted(res.items()))

    # remove dictionary items that have only one value
    res = {m:n for m,n in res.items() if len(n) != 1}

    for i in res.keys():
        
        chain = i[:1]
        resnum = i[2:]
        
        # here the first residue is a tRBI
        df1 = df[(df['Chain1'] == chain) & (df['AA1Res_num1'] == resnum)]

        # here the second residue is a tRBI (need to rearrange columns)
        df2 = df[(df['Chain2'] == chain) & (df['AA2Res_num2'] == resnum)]
        cols = df2.columns.tolist()
        cols = cols[:1] + cols[4:7] + cols[1:4] + cols[7:]
        df2 = df2[cols]
        
        dff = pd.concat([df1, df2])
        df_dict[i] = dff
        
    if len(df_dict.keys()) > 1:
        data = pd.concat(list(df_dict.values()))
    elif len(df_dict.keys()) == 1:
        data = pd.DataFrame(list(df_dict.items())[0][1])
    else:
        quit()
    
    data.columns = ['PDB_ID', 'chain1', 'trbi', 'aa1', 'chain2', 'prbi', 'aa2', 'dist', 'inttype']
    
    # tRBIs must have a side-chain participating
    data = data[data['aa1'] != 'O']
    data = data[data['aa1'] != 'C']
    data = data[data['aa1'] != 'CA']
    data = data[data['aa1'] != 'N']

    # tRBIs must interact with 2 different chains
    trbi_list = list(set(list(data['chain1'].astype(str) + " " + data['trbi'])))
    for i in range(len(trbi_list)):
        l = trbi_list[i].split(' ')
        res = l[1] + " " + l[2]
        d = data[data['trbi'] == res]
        chains_list = set(list(d['chain2']))
        if len(chains_list) < 2:
            data = data[data['trbi'] != res]

    return data

#-------------------------------------------------------------------------------

tr_rbi = true_rbi(df2)

if len(tr_rbi) > 0:

    outfilename1 = "/home/user/revathy/bifur/Programs/new_code/test_files/rbi/" + options.f1[-12:-8] + "_refined_rbi.csv"
    outfilename2 = "/home/user/revathy/bifur/Programs/new_code/test_files/rbi/" + options.f1[-12:-8] + "_unique_rbi_list.csv"

    tr_rbi[['trbi', 'res_num1']] = tr_rbi['trbi'].str.split(' ',expand=True)
    tr_rbi[['prbi', 'res_num2']] = tr_rbi['prbi'].str.split(' ',expand=True)
    cols = ['PDB_ID', 'chain1', 'trbi', 'res_num1', 'aa1', 'chain2', 'prbi', 'res_num2', 'aa2', 'dist', 'inttype']
    tr_rbi = tr_rbi[cols]

    trrbi = refine_interactions(tr_rbi)

    trrbi.sort_values(by = ['chain1', 'trbi', 'res_num1'], ascending = [True, True, True])

    # tr_rbi.to_csv('refined_rbi.csv', index = False)
    trrbi.drop_duplicates(subset = ['chain1', 'trbi', 'res_num1', 'aa1', 'chain2', 'prbi', 'res_num2', 'aa2'], keep = "first", inplace = True)
    trrbi.to_csv(outfilename1, index = False)

    # drop duplicates in rbi list: each amino acid pair will appear once, irrespective of interacting atoms
    trrbi.drop_duplicates(subset = ['chain1', 'trbi', 'res_num1', 'chain2', 'prbi', 'res_num2'], keep = "first", inplace = True)
    trrbi.to_csv(outfilename2, index = False)


print("--- %s seconds ---" % (time.time() - start_time))
