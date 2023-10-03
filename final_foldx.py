# to get hs prop
#------------------------------------------------------------------------------------------
import os, time
import pandas as pd
import numpy as np
from collections import defaultdict
from alive_progress import alive_bar
#------------------------------------------------------------------------------------------

# defining amino acid propensity with the following example:
#    a = no. of ALA residues that are RBIs
#    b = no. of ALA residues in dataset
#    c = no. of RBIs in dataset
#    d = no. of residues in dataset
#    aa propensity for ALA = (a/b)/(c/d)

def aa_propensity(inputArray, og_df, list_order):

    result = defaultdict(int)
    cd = (len(inputArray)/len(og_df))
    
    for element in list_order:
        a = inputArray.count(element)
        b = og_df.count(element)
        print(a, b, cd)
        if cd > 0:
            if b > 0:
                ab = a/b
                aaprop = ab/cd
                result[element] = aaprop
    
    return result

interface_aa = defaultdict()
interface_residues, trbi, prbi = [], [], []
aa_order = ['GLY', 'ALA', 'VAL', 'LEU', 'MET', 'ILE', 'PHE', 'TYR', 'TRP', 'SER', 'THR', 'CYS', 'PRO', 'ASN', 'GLN', 'LYS', 'ARG', 'HIS', 'ASP', 'GLU']

# df = pd.read_csv("/home/user/revathy/bifur/programs/new_code/test_files/Analysis/combined_hs_trbi.csv")
# df.drop_duplicates(inplace = True)
# pdbs = list(set(list(df['PDB_ID'])))

l = os.listdir("/home/user/revathy/bifur/programs/new_code/test_files/interface/rbi")

os.chdir("/home/user/revathy/bifur/programs/new_code/test_files/interface/rbi")
# with alive_bar(529, bar='solid', theme = 'smooth') as bar:
#     for i in range(len(pdbs)):
#         pdb = pdbs[i]
#         hs = df[df['PDB_ID'] == pdb]
#         res = list(hs.aa.str.split(" ", expand = True)[0])
#         trbi.extend(res)

#         # get residues from corresponding interface file
#         interface_file = [x for x in l if pdb in x]
#         if len(interface_file) > 0:
#             f_name = "/home/user/revathy/bifur/programs/new_code/test_files/interface/rbi/" + interface_file[0]
#             interface_df = pd.read_csv(f_name)
            
#             res1 = np.unique(list(interface_df[interface_df.columns[1:3]].apply(lambda x: ' '.join(x.astype(str)), axis=1)))
#             res2 = np.unique(list(interface_df[interface_df.columns[4:6]].apply(lambda x: ' '.join(x.astype(str)), axis=1)))
#             res = np.unique(np.hstack([res1, res2]))        

#             for i in res:
#                 r = i.split(" ")[1]
#                 interface_residues.append(r)
        
#         bar()

# trbi_aa_prop = aa_propensity(trbi, interface_residues, aa_order)
# pd.DataFrame.from_dict(trbi_aa_prop, orient='index').to_csv('/home/user/revathy/bifur/programs/new_code/test_files/final/hs_trbi_aa_prop.csv')

# for pRBIs
df = pd.read_csv("/home/user/revathy/bifur/programs/new_code/test_files/Analysis/combined_hs_prbi.csv")
df.drop_duplicates(inplace = True)
pdbs = list(set(list(df['PDB_ID'])))
interface_residues = []

with alive_bar(614, bar='solid', theme = 'smooth') as bar:
    for i in range(len(pdbs)):
        pdb = pdbs[i]

        hs = df[df['PDB_ID'] == pdb]
        res = list(hs.aa.str.split(" ", expand = True)[0])
        prbi.extend(res)

        # get residues from corresponding interface file
        interface_file = [x for x in l if pdb in x]
        if len(interface_file) > 0:
            f_name = "/home/user/revathy/bifur/programs/new_code/test_files/interface/rbi/" + interface_file[0]
            interface_df = pd.read_csv(f_name)
            res1 = np.unique(list(interface_df[interface_df.columns[1:3]].apply(lambda x: ' '.join(x.astype(str)), axis=1)))
            res2 = np.unique(list(interface_df[interface_df.columns[4:6]].apply(lambda x: ' '.join(x.astype(str)), axis=1)))
            res = np.unique(np.hstack([res1, res2]))        

            for i in res:
                r = i.split(" ")[1]
                interface_residues.append(r)
        bar()

prbi_aa_prop = aa_propensity(prbi, interface_residues, aa_order)
pd.DataFrame.from_dict(prbi_aa_prop, orient='index').to_csv('/home/user/revathy/bifur/programs/new_code/test_files/final/hs_prbi_aa_prop.csv')

# for interface
# df = pd.read_csv("/home/user/revathy/bifur/programs/new_code/test_files/Analysis/combined_hs_inter.csv")
# df.drop_duplicates(inplace = True)
# pdbs = list(set(list(df['PDB_ID'])))

# interface_residues = []
# inter = []

# with alive_bar(626, bar='solid', theme = 'smooth') as bar:
#     for i in range(len(pdbs)):
#         pdb = pdbs[i]

#         hs = df[df['PDB_ID'] == pdb]
#         res = list(hs.aa.str.split(" ", expand = True)[0])
#         inter.extend(res)

#         # get residues from corresponding interface file
#         interface_file = [x for x in l if pdb in x]
#         if len(interface_file) > 0:
#             f_name = "/home/user/revathy/bifur/programs/new_code/test_files/interface/rbi/" + interface_file[0]
#             interface_df = pd.read_csv(f_name)
#             res1 = np.unique(list(interface_df[interface_df.columns[1:3]].apply(lambda x: ' '.join(x.astype(str)), axis=1)))
#             res2 = np.unique(list(interface_df[interface_df.columns[4:6]].apply(lambda x: ' '.join(x.astype(str)), axis=1)))
#             res = np.unique(np.hstack([res1, res2]))        

#             for i in res:
#                 r = i.split(" ")[1]
#                 interface_residues.append(r)
#         bar()

# inter_aa_prop = aa_propensity(inter, interface_residues, aa_order)
# pd.DataFrame.from_dict(inter_aa_prop, orient='index').to_csv('/home/user/revathy/bifur/programs/new_code/test_files/final/hs_inter_aa_prop.csv')