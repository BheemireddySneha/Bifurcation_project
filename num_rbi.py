import os, time
import pandas as pd
from collections import defaultdict
from alive_progress import alive_bar
import numpy as np

#------------------------------------------------------------------------------------------

l = os.listdir("/home/user/revathy/bifur/programs/new_code/test_files/interface/rbi/")
pdb_dict = defaultdict()

os.chdir("/home/user/revathy/bifur/programs/new_code/test_files/rbi")
with alive_bar(635, bar='solid', theme = 'smooth') as bar:
    for file in os.listdir("/home/user/revathy/bifur/programs/new_code/test_files/rbi/"):
        filename = os.fsdecode(file)                                # precautionary : decoding filenames if encoded (like in Windows)
        if filename.endswith(".csv"):
            df = pd.read_csv(file)
            pdb = filename[0:4]
            n_trbi = np.unique(list(df[df.columns[1:4]].apply(lambda x: ' '.join(x.astype(str)), axis=1)))
            n_prbi = np.unique(list(df[df.columns[5:8]].apply(lambda x: ' '.join(x.astype(str)), axis=1)))

            pdb_dict[pdb] = [len(n_trbi), len(n_prbi), float(len(n_trbi)/len(n_prbi))]

        bar()

pd.DataFrame.from_dict(pdb_dict, orient='index').to_csv('/home/user/revathy/bifur/programs/new_code/test_files/final/num_rbi.csv')