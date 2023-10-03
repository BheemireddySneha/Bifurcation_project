import os
import pandas as pd

#------------------------------------------------------------------------------------------

int_list = []
interactions = {}

os.chdir("/home/user/revathy/bifur/programs/new_code/test_files/rbi/")
for file in os.listdir("/home/user/revathy/bifur/programs/new_code/test_files/rbi/"):
    filename = os.fsdecode(file)                        # precautionary : decoding filenames if encoded (like in Windows)
    if filename.endswith(".csv"):
        df = pd.read_csv(file)       
        inttype = list(df['inttype'])

        for element in inttype:
            ele = str(element).split(', ')
            ele = [i for i in ele if len(i) > 3]        # ensuring no spaces or commas are in the list
            int_list.extend(ele)

for i in int_list:
    if i not in interactions:
        interactions[i] = 1
    else:
        interactions[i] += 1


pd.DataFrame.from_dict(interactions, orient='index').to_csv('/home/user/revathy/bifur/programs/new_code/test_files/final/nofilter_rbi_int_type.csv')