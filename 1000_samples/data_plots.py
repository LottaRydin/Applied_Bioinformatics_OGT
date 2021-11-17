import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import random
from collections import Counter

#load data
os.chdir(os.path.dirname(os.path.realpath(__file__)))
df = pd.read_csv("data_whole.tsv", sep='\t') #compile_data_medium.py output
df = df.set_index(["OTU_ID"]) #set OTU col as the index
df.sort_index(inplace=True) #sort dataframe based on increasing index values
tempura_otu_df = pd.read_csv("matched_temp_multi_v2.tsv", sep='\t')
tempura_otu_df = tempura_otu_df.set_index(["OTU_id"])
metadata = pd.read_csv('temp_samples.tsv', sep='\t') #specify metadata file
metadata = metadata.drop(columns=['Unnamed: 0']) #drop unnecessary col
tempura = pd.read_csv("TEMPURA.csv") #load all tempura data

df.index.unique("OTU_ID") #unique OTU ID indices from the appended_df
df.groupby(level="OTU_ID").size() #count data samples per OTU, pandas series output
Counter(df.groupby(level="OTU_ID").size().tolist()) #how many OTUs have 1 data pt, 2 etc

#filtering
#filter out eukaryotes
df = df[~df.taxonomy.str.contains("Eukaryota")]

#Filter for number of samples per OTU
#df.drop(df[df["abundance"] < 3].index, inplace = True) #drop sample with less that 3 reads

#filter out OTUs with less than 10 observations
OTU_size = df.groupby(level="OTU_ID").size()
df = df.drop(OTU_size.index[OTU_size < 10].tolist())

##############################################################################
#plt violin
# for otu in random.sample(df.index.unique("OTU_ID").tolist(), 3):
#         
# otu_4572 = df.loc[4572, :].temp #species
# otu_232020 = df.loc[232020, :].temp
# otu_68134 = df.loc[68134, :].temp
# 
# labels = ['otu_4572', 'otu_232020', 'otu_68134']
# # Extract Figure and Axes instance
# fig, ax = plt.subplots()
# # Create a plot
# ax.violinplot([otu_4572, otu_232020, otu_68134], showmeans=True)
# ax.set_xticks(np.arange(1, len(labels) + 1))
# ax.set_xticklabels(labels)
# # Add title
# #ax.set_title('Violin Plot')
# plt.show()
# 
##############################################################################
#seaborn plot
plot_df = []
rand_otu = random.sample(df.index.unique("OTU_ID").tolist(), 10)
for indx in rand_otu:
    df_otu = df.loc[indx, :]
    plot_df.append(df_otu)

plot_df = pd.concat(plot_df)
plot_df.reset_index(inplace=True)

ax = sns.violinplot(x="taxonomy", y="temp", data=plot_df, scale="count", inner="stick") #If count, the width of the violins will be scaled by the number of observations in that bin. 
ax.tick_params(axis='x', rotation=90)

##############################################################################

len(list(set(df.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id"))))

def tempura_and_prediction(OTU_nr):
    plt.figure()
    plt.scatter(df.loc[OTU_nr, :].temp, df.loc[OTU_nr, :].rel_abundance, alpha=0.5)
    tempura_match = tempura_otu_df.loc[OTU_nr, :][['Tmin', 'Topt_ave', 'Tmax']]
    tempura_match = tempura_match.to_frame()
    tempura_match.rename(columns={ tempura_match.columns[0]: "temp" }, inplace = True)
    tempura_match["rel_abundance"] = [0, max(df.loc[OTU_nr, :].rel_abundance), 0]
    plt.scatter(tempura_match.temp, tempura_match.rel_abundance, alpha=0.5, c='red')
    plt.xlabel("Temperature")
    plt.ylabel("Relative abundance")
    plt.title(df.loc[OTU_nr, :].taxonomy.values[0])
    plt.show()
    del tempura_match

print("Out of", len(df.index.unique("OTU_ID")),"OTUs present in MGnify data", \
      len(list(set(df.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id")))), "were present in TEMPURA" )

rand_otu = random.sample(list(set(df.index.unique("OTU_ID")) & \
                              set(tempura_otu_df.index.unique("OTU_id"))), 5)
for otu in rand_otu:
    tempura_and_prediction(otu)
 
################################################################################ 
#plot how ratio of how close Topt is to Tmax in TEMPURA (but not all TEMPURA data, we use the file for OTU ID matching)
# tempura_tmax = tempura_otu_df[tempura_otu_df.taxonomy.str.contains("sk__Bacteria")]["Tmax"]
# tempura_tmin = tempura_otu_df[tempura_otu_df.taxonomy.str.contains("sk__Bacteria")]["Tmin"]
# tempura_topt = tempura_otu_df[tempura_otu_df.taxonomy.str.contains("sk__Bacteria")]["Topt_ave"] #get almost same as below if ignore archea
tempura_tmax = tempura_otu_df["Tmax"] #when a table with all TEMPURA data is used, get about same mean ratio
tempura_tmin = tempura_otu_df["Tmin"]
tempura_topt = tempura_otu_df["Topt_ave"]
tempura_ratio = (tempura_topt - tempura_tmin) /(tempura_tmax - tempura_tmin)

tempura_ratio.plot(kind='hist', bins=12,rwidth=0.9,color='#658b38', title = "TEMPURA mean ratio " \
                   + str("{:.5f}".format(tempura_ratio.mean()))+ " with std "+\
                       str("{:.5f}".format(tempura_ratio.std())))
plt.xlabel("ratio")
plt.show()
#tempura_ratio[tempura_ratio < 0] #how many ratios below 0, so that Topt is below Tmax
#tempura_ratio[(tempura_ratio > 0) & (tempura_ratio < 0.1)] #Topt very close to Tmin
###############################################################################

metadata.temperature.plot(kind='hist', bins=20, rwidth=0.9, color='#ae7181', \
                          title = "Temperatures of all samples" )
plt.xlabel("temperature")
plt.show()

df.temp.plot(kind='hist', bins=20, rwidth=0.9, color='#607c8e', \
             title = "Temperatures of all observations" )
plt.xlabel("temperature")
plt.show()

OTU_size.plot(kind='hist', bins=20, rwidth=0.9, color='#014d4e', \
              title = "Number of bservations of OTU" )
plt.xlabel("observation nr")
plt.show()
###############################################################################
temp_df = pd.DataFrame(columns=('OTU_ID', 'Tmin', 'Topt', 'Tmax', 'temp_Tmin', 'temp_Topt','temp_Tmax',))
for otu in df.index.unique("OTU_ID"):
    optimum_est = min(df.loc[199, :].temp) + (max(df.loc[199, :].temp) - min(df.loc[199, :].temp)) *0.66
    temps_est = [max(df.loc[199, :].temp), optimum_est, max(df.loc[199, :].temp)]
    








