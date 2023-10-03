import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import collections

#------------------------------------------------------------
# for RBIs

# df = pd.read_csv("/home/user/revathy/bifur/programs/new_code/test_files/all_rbi_combined.csv")
# # edges = list(df['trbi'].astype(str) + " " + df['prbi'].astype(str))
# # edges = list(df['trbi'].astype(str) + df['res_num1'].astype(str) + " " + df['prbi'].astype(str)+ df['res_num2'].astype(str))
# # edges = list(df['PDB_ID'].astype(str) + df['chain1'].astype(str) + df['trbi'].astype(str) + df['res_num1'].astype(str) + " " + df['PDB_ID'].astype(str) + df['chain2'].astype(str) + df['prbi'].astype(str)+ df['res_num2'].astype(str))
# edges = list(df['PDB_ID'].astype(str) + df['chain1'].astype(str) + df['trbi'].astype(str) + df['res_num1'].astype(str) + df['aa1'].astype(str) + " " + df['PDB_ID'].astype(str) + df['chain2'].astype(str) + df['prbi'].astype(str)+ df['res_num2'].astype(str) + df['aa2'].astype(str))
# edge_list = [x for x in edges if x != 'PDB_IDchain1trbires_num1 PDB_IDchain2prbiresnum2']

# # Load the protein structure network data
# G = nx.read_edgelist(edge_list)

# # Calculate the degree distribution
# degree_sequence = sorted([d for n, d in G.degree()], reverse=True)
# degreeCount = collections.Counter(degree_sequence)
# deg, cnt = zip(*degreeCount.items())

# # Plot the degree distribution
# # fig, ax = plt.subplots()
# # plt.bar(deg, cnt, width=0.80, color="b")
# # plt.title("Degree Histogram")
# # plt.ylabel("Count")
# # plt.xlabel("Degree")
# # ax.set_xticks([d + 0.4 for d in deg])
# # ax.set_xticklabels(deg)
# # plt.show()
# # plt.savefig("rbi_deg_dis4.pdf", dpi = 300)

# # Calculate centrality measures
# # degree_centrality = nx.degree_centrality(G)
# # closeness_centrality = nx.closeness_centrality(G)
# # betweenness_centrality = nx.betweenness_centrality(G)
# # eigenvector_centrality = nx.eigenvector_centrality(G)
# # # Plot the degree centrality values
# # plt.bar(range(len(degree_centrality)), list(degree_centrality.values()), align='center')
# # plt.xticks(range(len(degree_centrality)), list(degree_centrality.keys()))
# # plt.savefig("rbi_deg_centrality.pdf", dpi = 300)
# # print("rbi_deg_centrality")
# # # Plot the closeness centrality values
# # plt.bar(range(len(closeness_centrality)), list(closeness_centrality.values()), align='center')
# # plt.xticks(range(len(closeness_centrality)), list(closeness_centrality.keys()))
# # plt.savefig("rbi_closeness_centrality.pdf", dpi = 300)
# # print("rbi_closeness_centrality.pdf")
# # # Plot the betweenness centrality values
# # plt.bar(range(len(betweenness_centrality)), list(betweenness_centrality.values()), align='center')
# # plt.xticks(range(len(betweenness_centrality)), list(betweenness_centrality.keys()))
# # plt.savefig("rbi_betweenness_centrality.pdf", dpi = 300)
# # print("rbi_betweenness_centrality.pdf")
# # # Plot the eigenvector centrality values
# # plt.bar(range(len(eigenvector_centrality)), list(eigenvector_centrality.values()), align='center')
# # plt.xticks(range(len(eigenvector_centrality)), list(eigenvector_centrality.keys()))
# # plt.savefig("rbi_eig_centrality.pdf", dpi = 300)
# # print("rbi_eig_centrality.pdf")
# # # Calculate the clustering coefficient
# # clustering_coefficient = nx.clustering(G)
# # average_clustering_coefficient = nx.average_clustering(G)
# # # Plot the clustering coefficient values
# # plt.bar(range(len(clustering_coefficient)), list(clustering_coefficient.values()), align='center')
# # plt.xticks(range(len(clustering_coefficient)), list(clustering_coefficient.keys()))
# # plt.savefig("rbi_clustering_coeff.pdf", dpi = 300)
# # print("rbi_clustering_coeff.pdf")
# # Find the connected components
# connected_component_sizes = nx.connected_components(G)
# print("rbi_connected_component: ", connected_component_sizes)
# # Plot the connected component sizes
# # plt.bar(range(len(connected_component_sizes)), connected_component_sizes, align='center')
# # plt.xticks(range(len(connected_component_sizes)), range(1, len(connected_component_sizes)+1))
# # plt.savefig("rbi_connected_component.pdf", dpi = 300)
# # print("rbi_connected_component.pdf")
# # Calculate the degree assortativity
# degree_assortativity = nx.degree_assortativity_coefficient(G)
# print("RBI | degree assortativity: ", degree_assortativity)
# RBI | degree assortativity:  -0.06249109814877064

#------------------------------------------------------
# for interface residues

df = pd.read_csv("/home/revathy/lab/bifur/bifur/programs/new_code/test_files/interface/Interface_interactions/combined_interface.csv")
df.drop_duplicates()
df['AA1Res_num1'] = df['AA1Res_num1'].str.replace(' ', '')
df['AA2Res_num2'] = df['AA2Res_num2'].str.replace(' ', '')
# df['res1'] = df['AA1Res_num1'].str[0:3]
# df['res2'] = df['AA2Res_num2'].str[0:3]

# edges = list(df['res1'].astype(str) + " " + df['res2'].astype(str))
edges = list(df['PDB_ID'].astype(str) + df['Chain1'].astype(str) + df['AA1Res_num1'].astype(str) + df['Atom_type1'] + " " + df['PDB_ID'].astype(str) + df['Chain2'].astype(str) + df['AA2Res_num2'].astype(str) + df['Atom_type2'].astype(str))
edge_list = [x for x in edges if x != 'PDB_IDChain1AA1Res_num1 PDB_IDChain2AA2Res_num2']

# Load the protein structure network data
G = nx.read_edgelist(edge_list)

# Calculate the degree distribution
# degree_sequence = sorted([d for n, d in G.degree()], reverse=True)
# degreeCount = collections.Counter(degree_sequence)
# deg, cnt = zip(*degreeCount.items())

# # Plot the degree distribution
# # fig, ax = plt.subplots()
# # plt.bar(deg, cnt, width=0.80, color="b")
# # plt.title("Degree Histogram")
# # plt.ylabel("Count")
# # plt.xlabel("Degree")
# # ax.set_xticks([d + 0.4 for d in deg])
# # ax.set_xticklabels(deg)
# # plt.show()
# # plt.savefig("interface_deg_dis4.pdf", dpi = 300)

# # Calculate centrality measures
# degree_centrality = nx.degree_centrality(G)
# closeness_centrality = nx.closeness_centrality(G)
# betweenness_centrality = nx.betweenness_centrality(G)
# eigenvector_centrality = nx.eigenvector_centrality(G)
# # Plot the degree centrality values
# plt.bar(range(len(degree_centrality)), list(degree_centrality.values()), align='center')
# plt.xticks(range(len(degree_centrality)), list(degree_centrality.keys()))
# plt.savefig("interface_deg_centrality.pdf", dpi = 300)
# print("interface_deg_centrality.pdf")
# # Plot the closeness centrality values
# plt.bar(range(len(closeness_centrality)), list(closeness_centrality.values()), align='center')
# plt.xticks(range(len(closeness_centrality)), list(closeness_centrality.keys()))
# plt.savefig("interface_closeness_centrality.pdf", dpi = 300)
# print("interface_closeness_centrality.pdf")
# # Plot the betweenness centrality values
# plt.bar(range(len(betweenness_centrality)), list(betweenness_centrality.values()), align='center')
# plt.xticks(range(len(betweenness_centrality)), list(betweenness_centrality.keys()))
# plt.savefig("interface_betweenness_centrality.pdf", dpi = 300)
# print("interface_betweenness_centrality.pdf")
# # Plot the eigenvector centrality values
# plt.bar(range(len(eigenvector_centrality)), list(eigenvector_centrality.values()), align='center')
# plt.xticks(range(len(eigenvector_centrality)), list(eigenvector_centrality.keys()))
# plt.savefig("interface_eigen_centrality.pdf", dpi = 300)
# print("interface_eigen_centrality.pdf")
# # Calculate the clustering coefficient
# clustering_coefficient = nx.clustering(G)
# average_clustering_coefficient = nx.average_clustering(G)
# # Plot the clustering coefficient values
# plt.bar(range(len(clustering_coefficient)), list(clustering_coefficient.values()), align='center')
# plt.xticks(range(len(clustering_coefficient)), list(clustering_coefficient.keys()))
# plt.savefig("interface_clustering_coefficient.pdf", dpi = 300)
# print("interface_clustering_coefficient.pdf")
# # Find the connected components
# connected_component_sizes = nx.connected_components(G)
# print("connected components :", connected_components_sizes)
# Plot the connected component sizes
# plt.bar(range(len(connected_component_sizes)), connected_component_sizes, align='center')
# plt.xticks(range(len(connected_component_sizes)), range(1, len(connected_component_sizes)+1))
# plt.savefig("interface_connected_components.pdf", dpi = 300)
#print("interface_connected_components: ")
# Calculate the degree assortativity
degree_assortativity = nx.degree_assortativity_coefficient(G)
print("Interface | degree assortativity: ", degree_assortativity)
