# crosscheck IC data with other residues

#--------------------------------------------------------------------------------------

import pandas as pd
import os

#--------------------------------------------------------------------------------------

# for interface files
ic = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/Analysis/IC_final_inter.csv')
ic['res'] = ic.key.str.split(' ', expand = True)[0]
ic['tmp'] = ic.key.str.split(' ', expand = True)[1]
ic['chain'] = ic.tmp.str.split(':', expand = True)[0]
ic['resnum'] = ic.tmp.str.split(':', expand = True)[1]
ic['aa'] = ic['chain'].astype(str) + " " + ic['res'] + " " + ic['resnum'].astype(str)
ic.drop(['tmp', 'chain', 'resnum', 'Unnamed: 0'], axis = 1, inplace = True)

os.chdir("/home/user/revathy/bifur/programs/new_code/test_files/interface/naccess/")
# with alive_bar(597, bar='solid', theme = 'smooth') as bar:
for file in os.listdir("/home/user/revathy/bifur/programs/new_code/test_files/interface/naccess/"):
    filename = os.fsdecode(file)                        # precautionary : decoding filenames if encoded (like in Windows)

    df = pd.read_csv(filename).reset_index()
    pdb = filename[10:14]

    df['res1'] = df['Chain1'].astype(str) + " " + df['AA1Res_num1']
    df['res2'] = df['Chain2'].astype(str) + " " + df['AA2Res_num2']

    if pdb in ic['PDB_ID'].values:
        ic_r = ic[ic['PDB_ID'] == pdb]

        df1 = pd.concat([df['res1'], df['res2']]).rename('aa')
        df1.drop_duplicates(inplace = True)
        print(df1)
        dff = pd.merge(ic_r, df1, on = 'aa', how = 'inner')
        print(dff)
        quit()
        if len(dff) > 0:
            dff['res'] = dff.aa.str.split(' ', expand = True)[0]
            dff.to_csv('/home/user/revathy/bifur/programs/new_code/test_files/interface/hotspot/' + filename)
        else:
            print(filename)
    # bar()