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
        if b > 0:
            ab = a/b
            aaprop = ab/cd
            result[element] = aaprop
    
    return result

interface_aa = defaultdict()
interface_residues, trbi, prbi = [], [], []

l = os.listdir("/home/user/revathy/bifur/programs/new_code/test_files/interface/Interface_interactions/")

os.chdir("/home/user/revathy/bifur/programs/new_code/test_files/rbi")
with alive_bar(635, bar='solid', theme = 'smooth') as bar:
    for file in os.listdir("/home/user/revathy/bifur/programs/new_code/test_files/rbi/"):
        filename = os.fsdecode(file)                                # precautionary : decoding filenames if encoded (like in Windows)
        if filename.endswith(".csv"):
            df = pd.read_csv(file)

            df_trbi = list(set(list(df[df.columns[1:4]].apply(lambda x: ' '.join(x.astype(str)), axis=1))))
            for i in df_trbi:
                r = i.split(" ")[1]
                trbi.append(r)

            df_prbi = list(set(list(df[df.columns[5:8]].apply(lambda x: ' '.join(x.astype(str)), axis=1))))
            for i in df_prbi:
                r = i.split(" ")[1]
                prbi.append(r)

            pdb = filename[0:4]
            # get residues from corresponding interface file
            interface_file = [x for x in l if pdb in x]
            f_name = "/home/user/revathy/bifur/programs/new_code/test_files/interface/rbi/" + interface_file[0]
            interface_df = pd.read_csv(f_name)
            res1 = np.unique(list(interface_df[interface_df.columns[1:3]].apply(lambda x: ' '.join(x.astype(str)), axis=1)))
            res2 = np.unique(list(interface_df[interface_df.columns[4:6]].apply(lambda x: ' '.join(x.astype(str)), axis=1)))
            res = np.unique(np.hstack([res1, res2]))

            for i in res:
                r = i.split(" ")[1]
                interface_residues.append(r)
        bar()

aa_order = ['GLY', 'ALA', 'VAL', 'LEU', 'MET', 'ILE', 'PHE', 'TYR', 'TRP', 'SER', 'THR', 'CYS', 'PRO', 'ASN', 'GLN', 'LYS', 'ARG', 'HIS', 'ASP', 'GLU']

trbi_aa_prop = aa_propensity(trbi, interface_residues, aa_order)
prbi_aa_prop = aa_propensity(prbi, interface_residues, aa_order)

pd.DataFrame.from_dict(trbi_aa_prop, orient='index').to_csv('/home/user/revathy/bifur/programs/new_code/test_files/final/nofilter_trbi_aa_prop.csv')
pd.DataFrame.from_dict(prbi_aa_prop, orient='index').to_csv('/home/user/revathy/bifur/programs/new_code/test_files/final/nofilter_prbi_aa_prop.csv')