# code for the Kruskal-Wallis and Wilcoxon tests

import numpy as np
import pandas as pd
from scipy.stats import kruskal, wilcoxon, shapiro, ranksums


df1 = pd.read_csv("/home/revathy/lab/bifur/bifur/new_code/interface_files/results/interface_deg.csv")
# inter = list(df1['Degree'])
# i = [x.strip() for x in inter if x != 'Degree']
# print(shapiro(i))

df2 = pd.read_csv("/home/revathy/lab/bifur/bifur/new_code/interface_files/results/rbi_deg.csv")
rbi = list(df2['Degree'])
# print(shapiro(rbi))
# ShapiroResult(statistic=0.9658846259117126, pvalue=2.33442026153401e-15)
# not normally distributed

df3 = pd.read_csv("/home/revathy/lab/bifur/bifur/new_code/interface_files/results/prbi_deg.csv")
prbi = list(df3['Degree'])
# print(shapiro(prbi))
# ShapiroResult(statistic=0.9644635319709778, pvalue=1.987084137810915e-23)
# not normally distributed
# print("Wilcoxon, RBI ~ pRBI : ", wilcoxon(rbi, prbi))

# statistic, pvalue = ranksums(rbi, prbi)
# print(statistic, pvalue)
# 6.229381438317112 4.682803855197416e-10

df = pd.read_csv("~/Downloads/degrees.csv")
i = list(df['interface'])
statistic, pvalue = ranksums(rbi, i)
# print(statistic, pvalue)
# 12.658493644611639 1.0040749336909716e-36
statistic, pvalue = ranksums(prbi, i)
# print(statistic, pvalue)
# 10.63300557991136 2.0925785902584102e-26