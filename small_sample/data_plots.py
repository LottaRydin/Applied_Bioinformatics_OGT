import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import random

os.chdir(os.path.dirname(os.path.realpath(__file__)))

df = pd.read_csv("data_whole.tsv", sep='\t') #app_bioin_data_summary.py output
df = df.set_index(["OTU_ID"]) #set OTU col as the index
df.sort_index(inplace=True) #sort dataframe based on increasing index values

tempura_df = pd.read_csv("C:/Users/Gretka/Documents/MA_Bioinformatics/Applied_Bioinformatics/data_analysis/matched_temp.tsv", sep='\t')
tempura_df = tempura_df.set_index(["OTU_id"])

##############################################################################
#plt violin
otu_4572 = df.loc[4572, :].temp #species
otu_232020 = df.loc[232020, :].temp
otu_68134 = df.loc[68134, :].temp

labels = ['otu_4572', 'otu_232020', 'otu_68134']
# Extract Figure and Axes instance
fig, ax = plt.subplots()
# Create a plot
ax.violinplot([otu_4572, otu_232020, otu_68134], showmeans=True)
ax.set_xticks(np.arange(1, len(labels) + 1))
ax.set_xticklabels(labels)
# Add title
#ax.set_title('Violin Plot')
plt.show()

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

len(list(set(df.index.unique("OTU_ID")) & set(tempura_df.index.unique("OTU_id"))))

def tempura_and_prediction(OTU_nr):
    plt.figure()
    plt.scatter(df.loc[OTU_nr, :].temp, df.loc[OTU_nr, :].rel_abundance, alpha=0.5)
    tempura = tempura_df.loc[OTU_nr, :][['Tmin', 'Topt_ave', 'Tmax']]
    tempura = tempura.to_frame()
    tempura.rename(columns={ tempura.columns[0]: "temp" }, inplace = True)
    tempura["rel_abundance"] = [0, max(df.loc[OTU_nr, :].rel_abundance), 0]
    plt.scatter(tempura.temp, tempura.rel_abundance, alpha=0.5, c='red')
    plt.xlabel("Temperature")
    plt.ylabel("Relative abundance")
    plt.title(df.loc[OTU_nr, :].taxonomy.values[0])
    plt.show()
    del tempura

print("Out of", len(df.index.unique("OTU_ID")),"OTUs present in MGnify data", \
      len(list(set(df.index.unique("OTU_ID")) & set(tempura_df.index.unique("OTU_id")))), "were present in TEMPURA" )

rand_otu = random.sample(list(set(df.index.unique("OTU_ID")) & \
                              set(tempura_df.index.unique("OTU_id"))), 5)
for otu in rand_otu:
    tempura_and_prediction(otu)
 
################################################################################  
tempura_tmax = tempura_df[tempura_df.taxonomy.str.contains("sk__Bacteria")]["Tmax"]
tempura_tmin = tempura_df[tempura_df.taxonomy.str.contains("sk__Bacteria")]["Tmin"]
tempura_topt = tempura_df[tempura_df.taxonomy.str.contains("sk__Bacteria")]["Topt_ave"] 
tempura_ratio = (tempura_topt - tempura_tmin) /(tempura_tmax - tempura_tmin)

tempura_ratio.plot(kind='hist', bins=12, title = "TEMPURA mean ratio " \
                   + str("{:.5f}".format(tempura_ratio.mean()))+ " with std "+\
                       str("{:.5f}".format(tempura_ratio.std())) )

tempura_ratio[tempura_ratio < 0] #how many ratios below 0, so that Topt is below Tmax
tempura_ratio[(tempura_ratio > 0) & (tempura_ratio < 0.1)] #Topt very close to Tmin


