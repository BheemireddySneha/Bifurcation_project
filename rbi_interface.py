# check if RBI file's and interface file's PDB IDs match, and copy to new folder if they do

#------------------------------------------------------------------
import os
#------------------------------------------------------------------

matches = []

l1 = os.listdir('/home/user/revathy/bifur/programs/new_code/test_files/rbi/')
l2 = os.listdir('/home/user/revathy/bifur/programs/new_code/test_files/interface/Interface_interactions/')

# match filenames from the 2 folders and store it in a list if the PDB ID matches
for file in l1:
    pdb = file[0:4]
    for file2 in l2:
        if pdb in file2:
            matches.append(file2)

# copy the corresponding files to a new folder
for file in matches:
    os.system('cp /home/user/revathy/bifur/programs/new_code/test_files/interface/Interface_interactions/' + file + ' /home/user/revathy/bifur/programs/new_code/test_files/interface/rbi')