# to crosscheck FoldX hs results with other residues
#------------------------------------------------------------------------------------------
import pandas as pd
import os, time
from alive_progress import alive_bar
#------------------------------------------------------------------------------------------

# for interface residues
hs = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/Analysis/combined_hs_inter.csv').reset_index()
hs.drop_duplicates(inplace = True)

os.chdir("/home/user/revathy/bifur/programs/new_code/test_files/interface/naccess/")
with alive_bar(597, bar='solid', theme = 'smooth') as bar:
    for file in os.listdir("/home/user/revathy/bifur/programs/new_code/test_files/interface/naccess/"):
        filename = os.fsdecode(file)                        # precautionary : decoding filenames if encoded (like in Windows)

        df = pd.read_csv(filename).reset_index()
        pdb = filename[10:14]

        if pdb in hs['PDB_ID'].values:
            hs_r = hs[hs['PDB_ID'] == pdb]

            df1 = pd.concat([df['AA1Res_num1'], df['AA2Res_num2']]).rename('aa')
            df1.drop_duplicates(inplace = True)

            dff = pd.merge(hs_r, df1, on = 'aa', how = 'inner')
            dff = dff.loc[(dff['energy_change'].astype(float) > float(2)) | (dff['energy_change'].astype(float) < float(-2))]
            if len(dff) > 0:
                dff['res'] = dff.aa.str.split(' ', expand = True)[0]
                dff.to_csv('/home/user/revathy/bifur/programs/new_code/test_files/interface/hotspot/' + filename)
            else:
                print(filename)
    bar()


# for RBIs and pRBIs
hs_t = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/Analysis/combined_hs_trbi.csv').reset_index()
hs_t.drop_duplicates(inplace = True)

hs_p = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/Analysis/combined_hs_prbi.csv').reset_index()
hs_p.drop_duplicates(inplace = True)

os.chdir("/home/user/revathy/bifur/programs/new_code/test_files/rbi/naccess/")
with alive_bar(588, bar='solid', theme = 'smooth') as bar:
    for file in os.listdir("/home/user/revathy/bifur/programs/new_code/test_files/rbi/naccess/"):
        filename = os.fsdecode(file)                        # precautionary : decoding filenames if encoded (like in Windows)

        df = pd.read_csv(filename).reset_index()

        if len(df) > 0:
            pdb = filename[0:4]

            if pdb in hs_t['PDB_ID'].values:
                hs_r = hs_t[hs_t['PDB_ID'] == pdb]

                trbi = pd.Series(list(set(list(df[df.columns[4:6]].apply(lambda x: ' '.join(x.astype(str)), axis=1))))).rename('aa')
                prbi = pd.Series(list(set(list(df[df.columns[8:10]].apply(lambda x: ' '.join(x.astype(str)), axis=1))))).rename('aa')

                df_t = pd.merge(hs_r, trbi, on = 'aa', how = 'inner')
                df_t = df_t.loc[(df_t['energy_change'].astype(float) > float(2)) | (df_t['energy_change'].astype(float) < float(-2))]
                if len(df_t) > 0:
                    df_t['aa'] = df_t.aa.str.split(' ', expand = True)[0]
                    df_t.to_csv('/home/user/revathy/bifur/programs/new_code/test_files/rbi/trbi_hotspot/' + filename)
                else:
                    print("len(df_t) < 0", filename)
                    
                df_p = pd.merge(hs_r, prbi, on = 'aa', how = 'inner')
                df_p = df_p.loc[(df_p['energy_change'].astype(float) > float(2)) | (df_p['energy_change'].astype(float) < float(-2))]
                if len(df_p) > 0:
                    df_p['aa'] = df_p.aa.str.split(' ', expand = True)[0]
                    df_p.to_csv('/home/user/revathy/bifur/programs/new_code/test_files/rbi/prbi_hotspot/' + filename)
                else:
                    print("len(df_p) < 0", filename)
        else:
            print("len(df) < 0", filename)
        bar()