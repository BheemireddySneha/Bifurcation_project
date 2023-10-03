# convert Consurf data into table with aa, score and count columns for trbis

#--------------------------------------------------------------------------------------

import pandas as pd

#--------------------------------------------------------------------------------------

df_t = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/rbi/consurf/trbi_consurf_combined.csv')
df_p = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/rbi/consurf/prbi_consurf_combined.csv')

# drop rows where consurf scores are not available
df_t = df_t[df_t['Consurf_1'].notna()]
df_p = df_p[df_p['Consurf_2'].notna()]

df_t['trbi'] = df_t['PDB_ID_x'].astype(str) + " " + df_t['t'].astype(str) + " " + df_t['Consurf_1'].astype(str)
df_p['prbi'] = df_p['PDB_ID_x'].astype(str) + " " + df_p['p'].astype(str) + " " + df_p['Consurf_2'].astype(str)

# for RBIs
aa1 = list(df_t['trbi'])
aa = pd.DataFrame(list(set(aa1)))
aa.columns = ['aa']

# Split a single column into multiple using str.split()
aa[['pdb', 'chain', 'res_name', 'res_num', 'consurf_score']] = aa.aa.str.split(expand = True)
aa.drop(['aa', 'pdb', 'chain', 'res_num'], axis = 1, inplace = True)
consurf = aa.groupby(['res_name', 'consurf_score']).size()
pd.DataFrame(consurf).to_csv('test_files/trbi_consurf_no_na.csv')

# for pRBIs
aa2 = list(df_p['prbi'])
aa = pd.DataFrame(list(set(aa2)))
aa.columns = ['aa']

# Split a single column into multiple using str.split()
aa[['pdb', 'chain', 'res_name', 'res_num', 'consurf_score']] = aa.aa.str.split(expand = True)
aa.drop(['aa', 'pdb', 'chain', 'res_num'], axis = 1, inplace = True)

consurf = aa.groupby(['res_name', 'consurf_score']).size()
pd.DataFrame(consurf).to_csv('test_files/prbi_consurf_no_na.csv')
