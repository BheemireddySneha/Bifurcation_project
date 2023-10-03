# for binning data: https://stackoverflow.com/questions/45273731/binning-a-column-with-pandas

# convert FoldX data into table with aa, score and count columns for trbis

#--------------------------------------------------------------------------------------

import pandas as pd

#--------------------------------------------------------------------------------------

df = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/Analysis/combined_hs_inter.csv')
df = df[df['PDB_ID'] != 'PDB_ID']
print(len(df))
quit()

# Split a single column into multiple using str.split()
df[['res_name', 'res_num']] = df.aa.str.split(expand = True)
df.drop(['aa', 'res_num'], axis = 1, inplace = True)

pd.DataFrame(df).to_csv('test_files/hs_foldx_prbi.csv')
