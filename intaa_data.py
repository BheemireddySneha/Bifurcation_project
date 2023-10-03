# to crosscheck INTAA results with other residues
#------------------------------------------------------------------------------------------
import pandas as pd
import os, time
from alive_progress import alive_bar
#------------------------------------------------------------------------------------------

# for interface residues
intaa = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/Analysis/INTAA_inter_final.csv').reset_index()

os.chdir("/home/user/revathy/bifur/programs/new_code/test_files/interface/naccess/")
with alive_bar(597, bar='solid', theme = 'smooth') as bar:
    for file in os.listdir("/home/user/revathy/bifur/programs/new_code/test_files/interface/naccess/"):
        filename = os.fsdecode(file)                        # precautionary : decoding filenames if encoded (like in Windows)

        df = pd.read_csv(filename).reset_index()
        pdb = filename[10:14]

        if pdb in intaa['PDB_ID'].values:
            intaa_r = intaa[intaa['PDB_ID'] == pdb]

            df['res1'] = df[df.columns[3:5]].apply(lambda x: ' '.join(x.astype(str)), axis=1)
            df['res2'] = df[df.columns[6:8]].apply(lambda x: ' '.join(x.astype(str)), axis=1)
            df['res'] = df['res1'].astype(str) + " " + df['res2'].astype(str)

            intaa_r['res1'] = intaa_r['Chain_1'].astype(str) + " " + intaa_r['res_name1'].astype(str) + " " + intaa_r['res_1'].astype(str)
            intaa_r['res2'] = intaa_r['Chain_2'].astype(str) + " " + intaa_r['res_name2'].astype(str) + " " + intaa_r['res_2'].astype(str)
            intaa_r['res'] = intaa_r['res1'].astype(str) + " " + intaa_r['res2'].astype(str)

            dff = pd.merge(intaa_r, df, on = 'res', how = 'inner')
            if len(dff) > 0:
                dff.drop(['index_x', 'Unnamed: 0_x', 'Chain_1', 'res_1', 'res_name1', 'Chain_2', 'res_2', 'res_name2', 'res1_x', 'res2_x', 'index_y', 'Unnamed: 0_y', 'PDB_ID_y', 'Chain1', 'AA1Res_num1', 'Atom_type1', 'Chain2', 'AA2Res_num2', 'Atom_type2', 'res1_y', 'res2_y', 'Distance', 'IntType'], axis = 1, inplace = True)
                dff['aa1'] = dff.res.str.split(' ', expand = True)[1]
                dff['aa2'] = dff.res.str.split(' ', expand = True)[4]
                dff.drop_duplicates(inplace = True)

                if len(dff) > 0:
                    dff.to_csv('/home/user/revathy/bifur/programs/new_code/test_files/interface/intaa/' + filename)
                else:
                    print(filename)
        bar()

# for RBIs and pRBIs
intaa = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/Analysis/INTAA_final_rbi_csv.csv').reset_index()

os.chdir("/home/user/revathy/bifur/programs/new_code/test_files/rbi/naccess/")
with alive_bar(588, bar='solid', theme = 'smooth') as bar:
    for file in os.listdir("/home/user/revathy/bifur/programs/new_code/test_files/rbi/naccess/"):
        filename = os.fsdecode(file)                        # precautionary : decoding filenames if encoded (like in Windows)

        df = pd.read_csv(filename).reset_index()
        pdb = filename[0:4]

        if len(df) > 0:
            if pdb in intaa['PDB_ID'].values:
                intaa_r = intaa[intaa['PDB_ID'] == pdb]

                df['res1'] = df[df.columns[3:6]].apply(lambda x: ' '.join(x.astype(str)), axis=1)
                df['res2'] = df[df.columns[7:10]].apply(lambda x: ' '.join(x.astype(str)), axis=1)
                df['res'] = df['res1'].astype(str) + " " + df['res2'].astype(str)

                intaa_r['res1'] = intaa_r['Chain_1'].astype(str) + " " + intaa_r['res_name1'].astype(str) + " " + intaa_r['res_1'].astype(str)
                intaa_r['res2'] = intaa_r['Chain_2'].astype(str) + " " + intaa_r['res_name2'].astype(str) + " " + intaa_r['res_2'].astype(str)
                intaa_r['res'] = intaa_r['res1'].astype(str) + " " + intaa_r['res2'].astype(str)

                dff = pd.merge(intaa_r, df, on = 'res', how = 'inner')
                if len(dff) > 0:
                    dff.drop(['index_x', 'Unnamed: 0_x', 'Chain_1', 'res_1', 'res_name1', 'Chain_2', 'res_2', 'res_name2', 'res1_x', 'res2_x', 'index_y', 'Unnamed: 0_y', 'PDB_ID_y', 'chain1', 'res_num1', 'aa1', 'chain2', 'res_num2', 'aa2', 'res1_y', 'res2_y', 'dist', 'inttype'], axis = 1, inplace = True)
                    dff.drop_duplicates(inplace = True)

                    if len(dff) > 0:
                        dff.to_csv('/home/user/revathy/bifur/programs/new_code/test_files/rbi/intaa/' + filename)
                    else:
                        print(filename)
        bar()