# to crosscheck cs results with other residues
#------------------------------------------------------------------------------------------
import pandas as pd
import os, time
from alive_progress import alive_bar
#------------------------------------------------------------------------------------------

# for interface residues
# cs = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/Analysis/consurf_final_inter.csv').reset_index()

# os.chdir("/home/user/revathy/bifur/programs/new_code/test_files/interface/naccess/")
# with alive_bar(597, bar='solid', theme = 'smooth') as bar:
#     for file in os.listdir("/home/user/revathy/bifur/programs/new_code/test_files/interface/naccess/"):
#         filename = os.fsdecode(file)                        # precautionary : decoding filenames if encoded (like in Windows)

#         df = pd.read_csv(filename).reset_index()
#         pdb = filename[10:14]

#         if pdb in cs['PDB_ID'].values:
#             cs_r = cs[cs['PDB_ID'] == pdb]

#             df['res1'] = df[df.columns[3:5]].apply(lambda x: ' '.join(x.astype(str)), axis=1)
#             df['res2'] = df[df.columns[6:8]].apply(lambda x: ' '.join(x.astype(str)), axis=1)

#             cs_r['res1'] = cs_r['Chain1'].astype(str) + " " + cs_r['AA1Resname'].astype(str) + " " + cs_r['AA1Res_num'].astype(str)
#             cs_r['res2'] = cs_r['Chain2'].astype(str) + " " + cs_r['AA2Resname'].astype(str) + " " + cs_r['AA2Res_num'].astype(str)

#             df1 = pd.merge(cs_r, df, on = 'res1', how = 'inner')
#             df2 = pd.merge(cs_r, df, on = 'res2', how = 'inner')
#             dff = pd.concat([df1, df2], axis = 1)

#             if len(dff) > 0:
#                 dff.drop(['index_x', 'Unnamed: 0.1', 'Unnamed: 0_x', 'Chain1_x', 'AA1Resname', 'AA1Res_num', 'Chain2_x', 'AA2Resname', 'AA2Res_num',
#                         'index_y', 'Unnamed: 0_y', 'Chain1_y', 'AA1Res_num1', 'Atom_type1', 'Chain2_y',
#                         'AA2Res_num2', 'Atom_type2', 'Distance', 'IntType', 'res2_y', 'index_x',
#                         'Unnamed: 0.1', 'Unnamed: 0_x', 'Chain1_x', 'AA1Resname', 'PDB_ID_y', 
#                         'res2_x', 'PDB_ID_y', 'res1_y', 'AA1Res_num', 'Chain2_x', 'AA2Resname', 'AA2Res_num', 'res1_x', 'index_y', 'Unnamed: 0_y',
#                         'Chain1_y', 'AA1Res_num1', 'Atom_type1', 'Chain2_y', 'AA2Res_num2',
#                         'Atom_type2', 'Distance', 'IntType'], axis = 1, inplace = True)

#                 dff['aa1'] = dff.res1.str.split(' ', expand = True)[1]
#                 dff['aa2'] = dff.res2.str.split(' ', expand = True)[1]
#                 dff.drop_duplicates(inplace=True)
#                 if len(dff) > 0:
#                     dff.to_csv('/home/user/revathy/bifur/programs/new_code/test_files/interface/consurf/' + filename)
#                 else:
#                     print(filename)
#             else:
#                 print(filename)

#         bar()


# for interface residues
cs = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/Analysis/consurf_final_rbi.csv').reset_index()

os.chdir("/home/user/revathy/bifur/programs/new_code/test_files/rbi/naccess/")
with alive_bar(588, bar='solid', theme = 'smooth') as bar:
    for file in os.listdir("/home/user/revathy/bifur/programs/new_code/test_files/rbi/naccess/"):
        filename = os.fsdecode(file)                        # precautionary : decoding filenames if encoded (like in Windows)

        df = pd.read_csv(filename).reset_index()
        pdb = filename[0:4]

        if len(df) > 0:
            if pdb in cs['PDB_ID'].values:
                cs_r = cs[cs['PDB_ID'] == pdb]

                df['t'] = df[df.columns[3:6]].apply(lambda x: ' '.join(x.astype(str)), axis=1)
                df['p'] = df[df.columns[7:10]].apply(lambda x: ' '.join(x.astype(str)), axis=1)

                cs_r['t'] = cs_r[cs_r.columns[3:6]].apply(lambda x: ' '.join(x.astype(str)), axis = 1)
                cs_r['p'] = cs_r[cs_r.columns[7:10]].apply(lambda x: ' '.join(x.astype(str)), axis = 1)

                df1 = pd.merge(cs_r, df, on = 't', how = 'inner')
                df1.drop(['index_x', 'Unnamed: 0_x', 'chain1_x', 'res_num1_x', 'aa1_x', 'chain2_x', 'prbi_x', 'res_num2_x', 'aa2_x',
                        'dist_x', 'inttype_x', 'Consurf_2', 'p_x', 'index_y', 'Unnamed: 0_y', 'PDB_ID_y', 'chain1_y', 'trbi_y', 'res_num1_y', 'aa1_y',
                        'chain2_y', 'prbi_y', 'res_num2_y', 'aa2_y', 'dist_y', 'inttype_y', 'p_y'], axis = 1, inplace = True)
                df1.drop_duplicates(inplace = True)
                if len(df1) > 0:
                    df1.to_csv('/home/user/revathy/bifur/programs/new_code/test_files/rbi/consurf/trbi_' + filename)
                else:
                    print("len(trbi): ", filename)

                df2 = pd.merge(cs_r, df, on = 'p', how = 'inner')
                df2.drop(['index_x', 'Unnamed: 0_x', 'chain1_x', 'trbi_x', 'res_num1_x', 'aa1_x', 'chain2_x', 'res_num2_x', 'aa2_x',
                        'dist_x', 'inttype_x', 'Consurf_1', 't_x', 'index_y', 'Unnamed: 0_y', 'PDB_ID_y', 'chain1_y', 'trbi_y', 'res_num1_y', 'aa1_y',
                        'chain2_y', 'prbi_y', 'res_num2_y', 'aa2_y', 'dist_y', 'inttype_y', 't_y'], axis = 1, inplace = True)
                df2.drop_duplicates(inplace = True)
                if len(df2) > 0:
                    df2.to_csv('/home/user/revathy/bifur/programs/new_code/test_files/rbi/consurf/prbi_' + filename)
                else:
                    print("len(prbi): ", filename)
        else:
            print(filename)

        bar()