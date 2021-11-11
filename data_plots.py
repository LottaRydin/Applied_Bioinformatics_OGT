import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import random

os.chdir(os.path.dirname(os.path.realpath(__file__)))

df = pd.read_csv("data_whole.tsv", sep='\t')
df = df.set_index(["OTU_ID"]) #set OTU col as the index
df.sort_index(inplace=True) #sort dataframe based on increasing index values

#for indx in appended_df.index.unique("OTU_ID")[:10]: #just for the first 10 OTUs, grab data
#    print(appended_df.loc[indx, :])

#Filter for number of samples per OTU
#appended_df = appended_df.drop(appended_df[appended_df.rel_abundance < 0.00001].index)

#Filter for number of samples per OTU

otu_113966 = df.loc[113966, :].temp
otu_232020 = df.loc[232020, :].temp
otu_68134 = df.loc[68134, :].temp

labels = ['otu_113966', 'otu_232020', 'otu_68134']
# Extract Figure and Axes instance
fig, ax = plt.subplots()
# Create a plot
ax.violinplot([otu_113966, otu_232020, otu_68134], showmeans=True)
ax.set_xticks(np.arange(1, len(labels) + 1))
ax.set_xticklabels(labels)
# Add title
ax.set_title('Violin Plot')
plt.show()


##############################################################################

plot_df = []
rand_otu = random.sample(df.index.unique("OTU_ID").tolist(), 10)
for indx in rand_otu:
    df_otu = df.loc[indx, :]
    plot_df.append(df_otu)

plot_df = pd.concat(plot_df)
plot_df.reset_index(inplace=True)

ax = sns.violinplot(x="OTU_ID", y="temp", data=plot_df, scale="count", inner="stick") #If count, the width of the violins will be scaled by the number of observations in that bin. 




