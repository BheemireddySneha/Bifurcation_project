# convert Consurf data into table with aa, score and count columns

#--------------------------------------------------------------------------------------

import pandas as pd

#--------------------------------------------------------------------------------------

df = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/interface/consurf/interface_consurf_combined.csv')
df.drop(['PDB_ID_x.1', 'Consurf_1.1', 'Consurf_2.1'], axis = 1, inplace = True)
df.drop_duplicates(inplace = True)

# drop rows where consurf scores are not available
df = df[df['Consurf_1'].notna()]
df = df[df['Consurf_2'].notna()]

df['r1'] = df['PDB_ID_x'].astype(str) + " " + df['res1'].astype(str) + " " + df['Consurf_1'].astype(str)
df['r2'] = df['PDB_ID_x'].astype(str) + " " + df['res2'].astype(str) + " " + df['Consurf_2'].astype(str)

aa1 = list(df['r1'])
aa2 = list(df['r2'])
aa3 = aa1 + aa2

aa = pd.DataFrame(list(set(aa3)))
aa.columns = ['aa']

# Split a single column into multiple using str.split()
aa[['pdb', 'chain', 'res_name', 'res_num', 'consurf_score']] = aa.aa.str.split(expand = True)
aa.drop(['aa', 'pdb', 'chain', 'res_num'], axis = 1, inplace = True)

consurf = aa.groupby(['res_name', 'consurf_score']).size()
pd.DataFrame(consurf).to_csv('test_files/inter_consurf_no_na.csv')
