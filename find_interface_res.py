# This script identifies all interface residues interacting within a 5 angstrom cut-off radius.
# It also assigns interaction types based on the amino acids and specific atoms involved.

#------------------------------------------------------------------------------------------

from Bio.PDB import *
import numpy as np
import pandas as pd
from optparse import OptionParser
import time

#------------------------------------------------------------------------------------------
start_time = time.time()

df = pd.DataFrame()
usage = "USAGE: python3 find_interface_res.py --f1 file"
parser = OptionParser(usage=usage)
parser.add_option("--f1", help="molecule", dest="f1")
(options, args) = parser.parse_args()

# load the PDB file
file_ext = str(options.f1)[5:]
if file_ext == 'pdb' or file_ext == 'ent':
    p = PDBParser(QUIET = True)
elif file_ext == 'cif':
    p = MMCIFParser(QUIET = True)
structure = p.get_structure('molecule', options.f1)
mol = str(options.f1)[0:4]

# checking if number of chains is less than 3
chain_list = Selection.unfold_entities(structure, 'C')
if len(chain_list) < 3:
    print("only " + str(len(chain_list)) + " chain(s) in " + options.f1 + ", at least 3 chains required")
    quit()

atom_list = Selection.unfold_entities(structure, 'A')
atoms = Selection.unfold_entities(structure, 'A')

# find all atoms within 5 angstrom of each other
NS =  NeighborSearch(atom_list)
Chain_pairs = NS.search_all(5, 'C')

pairs = []
for i in range(0,len(Chain_pairs)):
    y = []
    y.append(Chain.Chain.get_id(Chain_pairs[i][0]))
    y.append(Chain.Chain.get_id(Chain_pairs[i][1]))
    pairs.append(tuple(y))

def contact_pairs(chain_1, chain_2, mol_name):

    pdb_id, dist, aa1, aa2, chainA, chainB = [], [], [], [], [], []
    aares_num1, aares_num2 = [], []

    chain_A = structure[0][str(chain_1)]
    chain_B = structure[0][str(chain_2)]
    atomsA = Selection.unfold_entities(chain_A, 'A')
    atomsB = Selection.unfold_entities(chain_B, 'A')

    for i in range(0, len(atomsA)):
        for j in range(0, len(atomsB)):
            A_name = atomsA[i].get_parent().get_resname()
            B_name = atomsB[j].get_parent().get_resname()

            if (is_aa(A_name)) and (is_aa(B_name)):
                d = atomsA[i] - atomsB[j]
                dist.append(d)
                aa1.append(atomsA[i].get_name())
                aares_num1.append(A_name+" "+str(atomsA[i].get_parent().get_id()[1]))
                chainA.append(chain_1)
                aa2.append(atomsB[j].get_name())
                aares_num2.append(B_name+" "+str(atomsB[j].get_parent().get_id()[1]))
                chainB.append(chain_2)
                pdb_id.append(mol_name)

    df1 = pd.DataFrame([pdb_id, chainA,aares_num1,aa1,chainB,aares_num2,aa2,dist]).transpose()
    interactatom1 = df1[df1[7] < 5]
    interactatom1.columns = ['PDB_ID', 'Chain1', 'AA1Res_num1','Atom_type1','Chain2', 'AA2Res_num2','Atom_type2','Distance']
    return interactatom1

contact_list = {}
for j in range(0,len(pairs)):
    iid = pairs[j][0]+pairs[j][1]
    contact_list[iid] = contact_pairs(pairs[j][0],pairs[j][1],mol)

chains=[]
for i in chain_list:
    chains.append(Chain.Chain.get_id(i))


def comm_inter(list_ch, ch_id):
    high = {}
    a1, a2 = list_ch, ch_id
    aino1, aino2, aino3 = [], [], []
    for inde in range(0, len(list_ch)):
        number = list_ch[inde].index(ch_id)
        if number == 0:
            col1 = contact_list[list_ch[inde]].AA1Res_num1
        else:
            col1 = contact_list[list_ch[inde]].AA2Res_num2

        for inde2 in range(0, len(list_ch)):
            if inde != inde2 and inde < inde2:
                number2 = list_ch[inde2].index(ch_id)
                if number2 == 0:
                    col2 = contact_list[list_ch[inde2]].AA1Res_num1
                else:
                    col2 = contact_list[list_ch[inde2]].AA2Res_num2
                    peru = list_ch[inde]+list_ch[inde2]
                aino1.append(list_ch[inde])
                aino2.append(list_ch[inde2])
                aino3.append(np.intersect1d(col1, col2))
                data2 = [aino1,aino2,aino3]
                high = pd.DataFrame([aino1, aino2, aino3]).transpose()

    return high

dict_new ={}
for ch in chains:
    lis_ch = []
    for k in contact_list.keys():
        if ch in k:
            lis_ch.append(k)
    dict_new[ch] = (comm_inter(lis_ch, ch))
rows_con = []

for chai in chains:
    ii = dict_new[chai]
    if len(ii) == 0:
        continue
    else:
        pass

    for entry in range(0, len(ii[2])):
        count_row = 0
        for contact in ii[2][entry]:
            rows = []
            for go in range(0,2):
                aa_resnum = contact_list[ii[go][count_row]]
                try:
                    res_num = aa_resnum.loc[aa_resnum['AA1Res_num1']==contact]
                    if len(res_num)!=0:
                        rows.append(res_num)
                    else:
                        res_num = aa_resnum.loc[aa_resnum['AA2Res_num2']==contact]
                        rows.append(res_num)
                except:
                    continue
            try:
                join_rows = pd.concat(rows)
            except:
                continue

            rows_con.append(join_rows)
        count_row+=1

if len(rows_con) < 1:
    print(" No bifurcated interactions found!!!")
    quit()
else :
    df = pd.concat(rows_con)

# defining interaction types
def inttype(df):

    aro = ['TRP', 'TYR', 'PHE']
    nonp = ['ALA', 'VAL', 'LEU', 'ILE', 'TRP', 'PRO', 'MET', 'PHE', 'TYR']
    sulph = ['MET', 'CYS']
    amino = ['GLN', 'ASN']
    anionic = ['HIS', 'GLU', 'ASP']
    cationic = ['ARG', 'LYS', 'HIS']

    df['IntType'] = ""

    for i in range(0, len(df)):
        r = df.iloc[i, :]
        dist = r['Distance']
        aa1 = df.iloc[i, 2][0:3]
        atom1 = df.iloc[i, 3][0:1]
        aa2 = df.iloc[i, 5][0:3]
        atom2 = df.iloc[i, 6][0:1]
        c = 0

        if dist < 5 or dist == 5:
            if ((atom1 == 'C' and aa1 in aro and aa2 in cationic) or (atom2 == 'C' and aa2 in aro and aa1 in cationic)):
                df.iloc[i,-1] = df.iloc[i,-1] + "cation-pi interaction, "
                c+=1
            if ((atom1 == 'C' and aa1 in aro and aa2 in anionic) or (atom2 == 'C' and aa2 in aro and aa1 in anionic)):
                df.iloc[i,-1] = df.iloc[i,-1] + "anion-pi interaction, "
                c+=1
            if (aa1 in anionic and aa2 in cationic) or (aa1 in cationic and aa2 in anionic):
                df.iloc[i,-1] = df.iloc[i,-1] + "salt-bridge, "
                c+=1
            if aa1 in nonp and aa2 in nonp:
                df.iloc[i,-1] = df.iloc[i,-1] + "hydrophobic interaction, "
                c+=1
            #if dist<5.3 or dist==5.3:
            if (aa1 in sulph and atom1 == 'S' and aa2 in aro and atom2 == 'C') or (aa1 in aro and atom1 == 'C' and aa2 in sulph and atom2 == 'S'):
                df.iloc[i,-1] = df.iloc[i,-1] + "pi-sulphur bond, "
                c+=1

        if dist>4.5:
            if atom1 == 'C' and atom2 == 'C' and aa1 in aro and aa2 in aro:
                df.iloc[i,-1] = df.iloc[i,-1] + "pi-pi interaction, "
                c+=1

        if dist<4.5 or dist==4.5:
            if ((atom1 == 'C' and aa1 in aro and aa2 in amino) or (atom2 == 'C' and aa2 in aro and aa1 in amino)):
                df.iloc[i,-1] = df.iloc[i,-1] + "amino-pi interaction, "
                c+=1
        if dist<4 or dist==4:
            if atom1 == 'S' or atom2 == 'S':
                df.iloc[i,-1] = df.iloc[i,-1] + "hydrogen bond, "
                c+=1

        if dist<3.5 or dist==3.5:
            if(atom1 == 'O' and atom2 == 'N') or (atom2 == 'O' and atom1 == 'N'):
                df.iloc[i,-1] = df.iloc[i,-1] + "hydrogen bond, "
                c+=1
        if dist<3.5 or dist==3.5:
            if(atom1 == 'H' and atom2 == 'N') or (atom2 == 'O' and atom1 == 'H') or (atom1 == 'N' and atom2 == 'H') or (atom2 == 'H' and atom1 == 'O'):
                df.iloc[i,-1] = df.iloc[i,-1] + "hydrogen bond, "
                c+=1

        if dist<2.04 or dist==2.04:
            if aa1 == 'CYS' and aa2 == 'CYS' and atom1 =='S' and atom2 == 'S':
                df.iloc[i,-1] = df.iloc[i,-1] + "disulphide bond, "
                c+=1

        if c == 0:
            df.iloc[i,-1] = "van der Waals interaction, "
    int_type = df['IntType']

    return df

inttable = inttype(df)
inttable.to_csv("interactions-" + options.f1 + ".csv", index = False)
print("--- %s seconds ---" % (time.time() - start_time))
