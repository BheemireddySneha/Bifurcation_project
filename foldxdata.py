#------------------------------------------------------------------------------------------
import pandas as pd
from collections import defaultdict
#------------------------------------------------------------------------------------------

# defining amino acid propensity with the following example:
#    a = no. of ALA residues that are RBIs
#    b = no. of ALA residues in dataset
#    c = no. of RBIs in dataset
#    d = no. of residues in dataset
#    aa propensity for ALA = (a/b)/(c/d)

def aa_propensity(inputArray, og_df, list_order):

    result = defaultdict(int)
    cd = (len(inputArray)/len(og_df))
    
    for element in list_order:
        a = inputArray.count(element)
        b = og_df.count(element)
        print(a, b, cd)
        if b > 0:
            ab = a/b
            aaprop = ab/cd
            result[element] = aaprop
    
    return result

interface_aa = defaultdict()
interface_residues, trbi, prbi = [], [], []

# df = pd.read_csv("/home/user/revathy/bifur/programs/new_code/test_files/all_rbi_combined.csv")
# df.drop_duplicates(inplace = True)
df = pd.read_csv("/home/user/revathy/bifur/programs/new_code/test_files/interface/naccess/all_interface_combined.csv")
df['aa1'] = df['Chain1'].astype(str) + " " + df['AA1Res_num1'].astype(str)
df['aa2'] = df['Chain2'].astype(str) + " " + df['AA2Res_num2'].astype(str)
df.drop(['Chain1', 'AA1Res_num1', 'Atom_type1', 'Chain2', 'AA2Res_num2', 'Atom_type2', 'Distance', 'IntType'], axis = 1, inplace = True)
df.drop_duplicates(inplace = True)

df_t = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/trbi_hotspot_combined.csv')
df_t.drop_duplicates(inplace = True)
df_p = pd.read_csv('/home/user/revathy/bifur/programs/new_code/test_files/prbi_hotspot_combined.csv')
df_p.drop_duplicates(inplace = True)

trbi_hs = list(df_t['aa'])
prbi_hs = list(df_p['aa'])

# trbi = list(df['trbi'])
# prbi = list(df['prbi'])

aa1 = list(df.aa1.str.split(' ', expand = True)[1])
aa2 = list(df.aa2.str.split(' ', expand = True)[1])
aa = aa1 + aa2

aa_order = ['GLY', 'ALA', 'VAL', 'LEU', 'MET', 'ILE', 'PHE', 'TYR', 'TRP', 'SER', 'THR', 'CYS', 'PRO', 'ASN', 'GLN', 'LYS', 'ARG', 'HIS', 'ASP', 'GLU']
# print('trbi', len(trbi))
# print('prbi', len(prbi))
# print(len([element for element in trbi if element in prbi]))
# quit()

trbi_aa_prop = aa_propensity(trbi_hs, aa, aa_order)
prbi_aa_prop = aa_propensity(prbi_hs, aa, aa_order)

pd.DataFrame.from_dict(trbi_aa_prop, orient='index').to_csv('/home/user/revathy/bifur/programs/new_code/test_files/interface_trbi_hs_aa_prop.csv')
pd.DataFrame.from_dict(prbi_aa_prop, orient='index').to_csv('/home/user/revathy/bifur/programs/new_code/test_files/interface_prbi_hs_aa_prop.csv')