import pandas as pd

df_rbi = pd.read_csv("/home/user/revathy/bifur/programs/new_code/test_files/all_interface_pdbs.csv")
df_rbi.drop(['Unnamed: 0'], axis = 1, inplace = True)
df_rbi.columns = ['pdb']
l = df_rbi['pdb'].values.tolist()
print(l)