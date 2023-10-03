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

# For interface residues, we can only calculate relative frequency.

interface_aa = defaultdict()
res = []

os.chdir("/home/user/revathy/bifur/programs/new_code/test_files/interface/Interface_interactions/")
with alive_bar(682, bar='solid', theme = 'smooth') as bar:
    for file in os.listdir("/home/user/revathy/bifur/programs/new_code/test_files/interface/Interface_interactions/"):
        filename = os.fsdecode(file)                                # precautionary : decoding filenames if encoded (like in Windows)
        if filename.endswith(".csv"):
            df = pd.read_csv(file)

            res1 = list(set(list(df[df.columns[1:3]].apply(lambda x: ' '.join(x.astype(str)), axis=1))))
            res2 = list(set(list(df[df.columns[4:6]].apply(lambda x: ' '.join(x.astype(str)), axis=1))))

            res1.extend(res2)

            for i in res1:
                r = i.split(" ")[1]
                if r not in interface_aa:
                    interface_aa[r] = 1
                else:
                    interface_aa[r] += 1
        bar()

total = sum(list(interface_aa.values()))
print(total)
# quit()

aa_order = ['GLY', 'ALA', 'VAL', 'LEU', 'MET', 'ILE', 'PHE', 'TYR', 'TRP', 'SER', 'THR', 'CYS', 'PRO', 'ASN', 'GLN', 'LYS', 'ARG', 'HIS', 'ASP', 'GLU', 'MOD']

aa_freq = defaultdict()

for i in aa_order:
    if i in interface_aa:
        aa_freq[i] = int(interface_aa[i])/total
    else:
        aa_freq[i] = 0

df = pd.DataFrame.from_dict(aa_freq, orient='index')
df.to_csv('/home/user/revathy/bifur/programs/new_code/test_files/final/interface_aa_count.csv')