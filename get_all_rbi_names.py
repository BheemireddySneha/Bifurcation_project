# to get a list of all PDBs with RBIs
#----------------------------
import os
import pandas as pd
#----------------------------

pdb = []
l = os.listdir('/home/user/revathy/bifur/programs/new_code/test_files/rbi/naccess/')

os.chdir('/home/user/revathy/bifur/programs/new_code/test_files/rbi/naccess/')
for i in range(len(l)):
    df = pd.read_csv(l[i])
    if len(df) > 0:
        pdb.append(l[i][0:4])

# print(pdb)
pd.DataFrame(pdb).to_csv('/home/user/revathy/bifur/programs/new_code/test_files/all_rbi_pdbs.csv')
