import pandas as pd
#----------------------------------------------

df_ppcheck = pd.read_csv("/home/user/revathy/bifur/programs/new_code/test_files/ppcheck_hotspot/results/extracted_results/combined_ppcheck_hotspots.csv")
df_foldx_trbi = pd.read_csv("/home/user/revathy/bifur/programs/new_code/test_files/Analysis/combined_hs_trbi.csv")

df_ppcheck['res'] = df_ppcheck['PDB_ID'].astype(str) + " " + df_ppcheck['aa'].astype(str)
df_foldx_trbi['res'] = df_foldx_trbi['PDB_ID'].astype(str) + " " + df_foldx_trbi['aa'].astype(str)

s1 = pd.merge(df_ppcheck, df_foldx_trbi, how='inner', on=['res'])
s1.drop(['Unnamed: 0', 'res', 'PDB_ID_y', 'aa_y'], axis = 1, inplace = True)
s1.columns = ['PDB_ID', 'aa', 'chain', 'energy_change']
# print(s1.head())
# print("ppcheck: ", len(df_ppcheck))
print("FOR RBI")
print("foldx: ", len(df_foldx_trbi))
print("intersect: ", len(s1))


df_foldx_prbi = pd.read_csv("/home/user/revathy/bifur/programs/new_code/test_files/Analysis/combined_hs_prbi.csv")
df_foldx_prbi['res'] = df_foldx_prbi['PDB_ID'].astype(str) + " " + df_foldx_prbi['aa'].astype(str)

s1 = pd.merge(df_ppcheck, df_foldx_prbi, how='inner', on=['res'])
s1.drop(['Unnamed: 0', 'res', 'PDB_ID_y', 'aa_y'], axis = 1, inplace = True)
s1.columns = ['PDB_ID', 'aa', 'chain', 'energy_change']
# print(s1.head())
# print("ppcheck: ", len(df_ppcheck))
print("FOR INTERFACE RESIDUES")
print("foldx: ", len(df_foldx_prbi))
print("intersect: ", len(s1))


df_foldx_inter = pd.read_csv("/home/user/revathy/bifur/programs/new_code/test_files/Analysis/combined_hs_inter.csv")
df_foldx_inter['res'] = df_foldx_inter['PDB_ID'].astype(str) + " " + df_foldx_inter['aa'].astype(str)

s1 = pd.merge(df_ppcheck, df_foldx_inter, how='inner', on=['res'])
s1.drop(['Unnamed: 0', 'res', 'PDB_ID_y', 'aa_y'], axis = 1, inplace = True)
s1.columns = ['PDB_ID', 'aa', 'chain', 'energy_change']
# print(s1.head())
# print("ppcheck: ", len(df_ppcheck))
print("FOR PRBI")
print("foldx: ", len(df_foldx_inter))
print("intersect: ", len(s1))