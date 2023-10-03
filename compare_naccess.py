# cross-check interface, RBI and pRBI residues with NACCESS records and remove buried residues
#------------------------------------------------------------------------------------------

import pandas as pd
import os, time
from alive_progress import alive_bar
#------------------------------------------------------------------------------------------

# interface file with NACCESS results
df = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/Analysis/naccess_final_inter.csv')
df.reset_index()
df['res'] = df['chain_final'].astype(str) + " " + df['Resname'] + " " + df['Resnum'].astype(str)
nac_inter = df[df['Buried'] == True]

# Interface residue files
os.chdir("/home/user/revathy/bifur/programs/new_code/test_files/interface/rbi/")

for file in os.listdir("/home/user/revathy/bifur/programs/new_code/test_files/interface/rbi/"):
    filename = os.fsdecode(file)                        # precautionary : decoding filenames if encoded (like in Windows)

    df_i = pd.read_csv(filename)
    pdb = filename[10:14]
    nac_inter_df = nac_inter[nac_inter['Pdb_id'] == pdb]

    if len(nac_inter_df) == 0:
        continue
    else:
        df_i['res1'] = df_i['Chain1'] + " " + df_i['AA1Res_num1']
        df_i['res2'] = df_i['Chain2'] + " " + df_i['AA2Res_num2']

        dff = pd.merge_ordered(df_i, nac_inter_df, left_on = 'res1', right_on = 'res', how = 'inner')
        dfff = pd.merge_ordered(dff, nac_inter_df, left_on = 'res2', right_on = 'res', how = 'inner')
        dfff.drop(['res1', 'res2', 'Unnamed: 0_x', 'Pdb_id_x', 'chain_final_x', 'Resnum_x', 'Resname_x', 'Buried_x', 'res_x', 'Unnamed: 0_y', 'Pdb_id_y', 'chain_final_y', 'Resnum_y', 'Resname_y', 'Buried_y', 'res_y'], axis = 1, inplace = True)

    dfff.to_csv("/home/user/revathy/bifur/programs/new_code/test_files/interface/naccess/" + filename)

# RBIs and pRBIs

df_t = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/Analysis/naccess_final_RBI.csv')
df_t.reset_index()
df_t['res'] = df_t['chain_final'].astype(str) + " " + df_t['Resname'] + " " + df_t['Resnum'].astype(str)
nac_trbi = df_t[df_t['Buried'] == True]
print("nac-trbi", len(nac_trbi))

df_p = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/Analysis/naccess_final_iRBI.csv')
df_p.reset_index()
df_p['res'] = df_p['chain_final'].astype(str) + " " + df_p['Resname'] + " " + df_p['Resnum'].astype(str)
nac_prbi = df_p[df_p['Buried'] == True]
print("nac-prbi", len(nac_prbi))

os.chdir("/home/user/revathy/bifur/programs/new_code/test_files/rbi/")

for file in os.listdir("/home/user/revathy/bifur/programs/new_code/test_files/rbi/"):
    filename = os.fsdecode(file)                        # precautionary : decoding filenames if encoded (like in Windows)
    if filename.endswith(".csv"):
        df_i = pd.read_csv(filename)
        pdb = filename[0:4]
        nac_t_df = nac_trbi[nac_trbi['Pdb_id'] == pdb]
        nac_p_df = nac_prbi[nac_prbi['Pdb_id'] == pdb]

        if len(nac_t_df) == 0 or len(nac_p_df) == 0:
            continue
        else:
            df_i['res1'] = df_i['chain1'].astype(str) + " " + df_i['trbi'].astype(str) + " " + df_i['res_num1'].astype(str)
            df_i['res2'] = df_i['chain2'].astype(str) + " " + df_i['prbi'].astype(str) + " " + df_i['res_num2'].astype(str)

            dff = pd.merge_ordered(df_i, nac_t_df, left_on = 'res1', right_on = 'res', how = 'inner')
            dfff = pd.merge_ordered(dff, nac_p_df, left_on = 'res2', right_on = 'res', how = 'inner')
            dfff.drop(['res1', 'res2', 'Unnamed: 0_x', 'Pdb_id_x', 'chain_final_x', 'Resnum_x', 'Resname_x', 'Buried_x', 'res_x', 'Unnamed: 0_y', 'Pdb_id_y', 'chain_final_y', 'Resnum_y', 'Resname_y', 'Buried_y', 'res_y'], axis = 1, inplace = True)
            print(dfff)
            quit()
        # dfff.to_csv("/home/user/revathy/bifur/programs/new_code/test_files/rbi/naccess/" + filename)





