import pandas as pd
import os

l = os.listdir("/home/user/revathy/bifur/programs/new_code/test_files/ppcheck_hotspot/results/extracted_results/")
os.chdir("/home/user/revathy/bifur/programs/new_code/test_files/ppcheck_hotspot/results/extracted_results/")

pdbs, aa, chain = [], [], []

for i in range(len(l)):

    pdbid = l[i][0:4]
    # df = pd.read_csv(l[i])
    
    try:
        df = pd.read_csv(l[i])
        for j in range(len(df)):
            pdbs.append(pdbid)
            aa.append(df.iloc[j, 2] + " " + df.iloc[j, 1].astype(str))
            chain.append(df.iloc[j, -1])
    except:
        continue

df = pd.DataFrame([pdbs, aa, chain]).T
df.columns = ['PDB_ID', 'aa', 'chain']
df['res'] = df['PDB_ID'].astype(str) + " " + df['aa'].astype(str) + " " + df['chain'].astype(str)
df.drop_duplicates(['res'], inplace = True)
df.drop(['res'], axis = 1, inplace = True)
df.to_csv("/home/user/revathy/bifur/programs/new_code/test_files/ppcheck_hotspot/results/extracted_results/combined_ppcheck_hotspots.csv")

