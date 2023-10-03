# make the table for chain list

#--------------------------------------------------------------------------------------

import pandas as pd
from collections import defaultdict

#--------------------------------------------------------------------------------------

df = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/chains_list', sep = ' ')
df.columns = ['pdb_id', 'chains']

chains = list(df['chains'])
chain_dict = defaultdict()
for i in range(len(chains)):
    num = chains[i]
    if num not in chain_dict:
        chain_dict[num] = chains.count(num)

pd.DataFrame(chain_dict.items(), columns = ['chain_num', 'count']).to_csv('test_files/chain_count.csv')

