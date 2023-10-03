# sort nox data
#-----------------------------------------------------
import pandas as pd
import numpy as np
from collections import defaultdict
#-----------------------------------------------------

df = pd.read_csv("/home/user/revathy/bifur/programs/new_code/test_files/Analysis/nox_final_rbi_resname.csv")
df.drop(['Unnamed: 0'], axis = 1, inplace = True)

pdbs = np.unique(list(df['pdb']))
interface = defaultdict()
partial, obligate, non_obligate = 0, 0, 0
for i in range(len(pdbs)):
    dff = df[df['pdb'] == pdbs[i]]
    if 'False' in list(dff['Obligatory'].astype(str)):
        if 'True' in list(dff['Obligatory'].astype(str)):
            interface[pdbs[i]] = 'partial'
            partial += 1
        else:
            interface[pdbs[i]] = 'non-obligate'
            non_obligate += 1
    else:
        interface[pdbs[i]] = 'obligate'
        obligate += 1

# data = pd.DataFrame.from_dict(interface, orient = 'index').to_csv('/home/user/revathy/bifur/programs/new_code/test_files/final/nox.csv')
print('partial', partial)
print('non-obligate', non_obligate)
print('obligate', obligate)