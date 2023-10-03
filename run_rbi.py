import os

l = os.listdir("/home/user/revathy/bifur/Programs/new_code/test_files/interface/Interface_interactions")

for i in l:
    os.system("python3 find_rbi_no_hydrogen_filter.py --f1 /home/user/revathy/bifur/Programs/new_code/test_files/interface/Interface_interactions/" + i)

# for i in l:
#     os.system("python3 find_rbi.py --f1 /home/user/revathy/bifur/Programs/new_code/test_files/interface/Interface_interactions/" + i)