# convert Consurf data into table with aa, score and count columns

#--------------------------------------------------------------------------------------

import pandas as pd

#--------------------------------------------------------------------------------------

# for interface
df = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/Analysis/consurf_final_inter.csv')
df.drop(['Unnamed: 0.1', 'Unnamed: 0'], axis = 1, inplace = True)
df.drop_duplicates(inplace = True)

# drop rows where consurf scores are not available
df = df[df['Consurf_1'].notna()]
df = df[df['Consurf_2'].notna()]

# remove all scores greater than 9
df = df[df['Consurf_1'] <= 9.0]
df = df[df['Consurf_2'] <= 9.0]

df['r1'] = df['PDB_ID'].astype(str) + " " + df['Chain1'].astype(str) + " " + df['AA1Resname'].astype(str) + " " + df['AA1Res_num'].astype(str) + " " + df['Consurf_1'].astype(str)
df['r2'] = df['PDB_ID'].astype(str) + " " + df['Chain2'].astype(str) + " " + df['AA2Resname'].astype(str) + " " + df['AA2Res_num'].astype(str) + " " + df['Consurf_2'].astype(str)

aa1 = list(df['r1'])
aa2 = list(df['r2'])
aa3 = aa1 + aa2

aa = pd.DataFrame(list(set(aa3)))
aa.columns = ['aa']

# Split a single column into multiple using str.split()
aa[['pdb', 'chain', 'res_name', 'res_num', 'consurf_score']] = aa.aa.str.split(expand = True)
aa.drop(['aa', 'pdb', 'chain', 'res_num'], axis = 1, inplace = True)

consurf = aa.groupby(['res_name', 'consurf_score']).size()
pd.DataFrame(consurf).to_csv('test_files/final/inter_consurf_no_na.csv')

#-----------------------------------------
# for RBIs and pRBIs
df_t = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/Analysis/consurf_final_rbi.csv')

# drop rows where consurf scores are not available
df_t = df_t[df_t['Consurf_1'].notna()]
df_t = df_t[df_t['Consurf_2'].notna()]

df_t['t'] = df_t['PDB_ID'].astype(str) + " " + df_t['chain1'].astype(str) + " " + df_t['trbi'].astype(str) + " " + df_t['res_num1'].astype(str) + " " + df_t['Consurf_1'].astype(str)

# for RBIs
aa1 = list(df_t['t'])
aa = pd.DataFrame(list(set(aa1)))
aa.columns = ['aa']

# Split a single column into multiple using str.split()
aa[['pdb', 'chain', 'res_name', 'res_num', 'consurf_score']] = aa.aa.str.split(expand = True)
aa.drop(['aa', 'pdb', 'chain', 'res_num'], axis = 1, inplace = True)
consurf = aa.groupby(['res_name', 'consurf_score']).size()
pd.DataFrame(consurf).to_csv('test_files/final/trbi_consurf_no_na.csv')

# for pRBIs
df_t['p'] = df_t['PDB_ID'].astype(str) + " " + df_t['chain2'].astype(str) + " " + df_t['prbi'].astype(str) + " " + df_t['res_num2'].astype(str) + " " + df_t['Consurf_2'].astype(str)
aa2 = list(df_t['p'])
aa = pd.DataFrame(list(set(aa2)))
aa.columns = ['aa']

# Split a single column into multiple using str.split()
aa[['pdb', 'chain', 'res_name', 'res_num', 'consurf_score']] = aa.aa.str.split(expand = True)
aa.drop(['aa', 'pdb', 'chain', 'res_num'], axis = 1, inplace = True)

consurf = aa.groupby(['res_name', 'consurf_score']).size()
pd.DataFrame(consurf).to_csv('test_files/final/prbi_consurf_no_na.csv')

