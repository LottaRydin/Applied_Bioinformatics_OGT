import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import random
from collections import Counter

##load data
os.chdir(os.path.dirname(os.path.realpath(__file__))) #set current dir to script location
tempura_otu_df = pd.read_csv("matched_temp_multi_v2.tsv", sep='\t')
tempura_otu_df = tempura_otu_df.set_index(["OTU_id"])
tempura = pd.read_csv("TEMPURA.csv") #load entire tempura db
metadata = pd.read_csv('temp_samples.tsv', sep='\t') #specify metadata file
metadata = metadata.drop(columns=['Unnamed: 0']) #drop unnecessary col
pipeline_map = pd.read_csv('both_pipes.tsv', sep='\t')
df = pd.read_csv("data_whole.tsv", sep='\t') #compile_data_medium.py output

## converting v5 otu id to v4.1
#which v5 otu not in pipeline_map, remove this data
list_1 = pipeline_map.OTU_id_5.tolist() #all v5 otu id in map file
list_2 = df.loc[df.pipeline == 5, :].OTU_ID.unique() #list of v5 otu from tsv files
missing_v5_id = list(set(list_2).difference(list_1))
df = df.drop(df[(df['pipeline'] == 5) & (df['OTU_ID'].isin(missing_v5_id))].index) #removed 400 rows because no mtching v5 otu id in the map
#now convert v5 otu id to v4.1
otu_dict = dict(zip(pipeline_map.OTU_id_5, pipeline_map.OTU_id_4_1)) #dict which links v5 id t ov4.1 id
df.loc[df.pipeline == 5, 'OTU_ID'] = df['OTU_ID'].map(otu_dict) #in rows where pipeline == 5, convert otu id based on dict
df["OTU_ID"] = df["OTU_ID"].astype(int, errors = 'raise')
df = df.set_index(["OTU_ID"]) #set OTU col as the index
df.sort_index(inplace=True) #sort dataframe based on increasing index values
del list_1
del list_2
del otu_dict

##fitering
#filter out eukaryotes
df = df[~df.taxonomy.str.contains("Eukaryota")]
#Filter for number of samples per OTU
#df.drop(df[df["abundance"] < 3].index, inplace = True) #drop sample with less that 3 reads
#filter out OTUs with less than 10 observations
OTU_size = df.groupby(level="OTU_ID").size()
df = df.drop(OTU_size.index[OTU_size < 10].tolist())

#data summary
df.index.unique("OTU_ID") #unique OTU ID indices df
df.groupby(level="OTU_ID").size() #count data samples per OTU, pandas series output
Counter(df.groupby(level="OTU_ID").size().tolist()) #how many OTUs have 1 data pt, 2 etc
##############################################################################
#plt violin for specific otu s
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
#seaborn violin plot
plot_df = []
rand_otu = random.sample(df.index.unique("OTU_ID").tolist(), 10) #pick 10 random otu ids
for indx in rand_otu:
    df_otu = df.loc[indx, :]
    plot_df.append(df_otu)

plot_df = pd.concat(plot_df)
plot_df.reset_index(inplace=True)

ax = sns.violinplot(x="taxonomy", y="temp", data=plot_df, scale="count", inner="stick") #If count, the width of the violins will be scaled by the number of observations in that bin. 
ax.tick_params(axis='x', rotation=90)

##############################################################################
#plot scatter plot to compare otu data to tempura
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
                              set(tempura_otu_df.index.unique("OTU_id"))), 5) #pick random ots which have data also in tempura
for otu in rand_otu:
    tempura_and_prediction(otu)
 
################################################################################ 
#plot ratio of how close Topt is to Tmax in TEMPURA (but not all TEMPURA data, we use the file with 4.1 otu matches)
# tempura_topt = tempura_otu_df[tempura_otu_df.taxonomy.str.contains("sk__Bacteria")]["Topt_ave"] #get almost same as below if ignore archea
tempura_tmax = tempura_otu_df["Tmax"] #when a table with all TEMPURA data is used, get about same mean ratio
tempura_tmin = tempura_otu_df["Tmin"]
tempura_topt = tempura_otu_df["Topt_ave"]
tempura_ratio = (tempura_topt - tempura_tmin) /(tempura_tmax - tempura_tmin)
del tempura_tmax
del tempura_tmin
del tempura_topt

tempura_ratio.plot(kind='hist', bins=12,rwidth=0.9,color='#658b38', title = "TEMPURA mean ratio " \
                   + str("{:.5f}".format(tempura_ratio.mean()))+ " with std "+\
                       str("{:.5f}".format(tempura_ratio.std())))
plt.xlabel("ratio")
plt.show()
#tempura_ratio[tempura_ratio < 0] #how many ratios below 0, so that Topt is below Tmax
#tempura_ratio[(tempura_ratio > 0) & (tempura_ratio < 0.1)] #Topt very close to Tmin
###############################################################################
#some summary histograms
df.groupby(level="OTU_ID").size().plot(kind='hist', bins=30, rwidth=0.9, color='#014d4e', \
              title = "Number of observations of OTU" ) #df.groupby(level="OTU_ID")
plt.xlabel("observation nr")
plt.show()

metadata.temperature.plot(kind='hist', bins=20, rwidth=0.9, color='#ae7181', \
                          title = "Temperatures of all samples" )
plt.xlabel("temperature")
plt.show()

df.temp.plot(kind='hist', bins=20, rwidth=0.9, color='#607c8e', \
             title = "Temperatures of all observations" )
plt.xlabel("temperature")
plt.show()

###############################################################################
#combining MGnify and TEMPURA temp data
cols = ['OTU_ID', 'Tmin', 'Topt', 'Tmax']
temps_pred = []
for otu in df.index.unique("OTU_ID"):
    optimum_est = min(df.loc[otu, :].temp) + (max(df.loc[otu, :].temp) - min(df.loc[otu, :].temp)) *0.66
    temps_est = [otu, min(df.loc[otu, :].temp), optimum_est, max(df.loc[otu, :].temp)]
    temps_pred.append(temps_est)

temps_pred = pd.DataFrame(temps_pred, columns=cols)
del cols
del temps_est
temps_pred = temps_pred.set_index(["OTU_ID"]) #set OTU col as the index
temps_pred.sort_index(inplace=True) #sort dataframe based on increasing index values

tempura_matches = list(set(df.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id")))
temps_compare = temps_pred.loc[tempura_matches]
temps_compare = temps_compare.reindex(columns=['Tmin', 'Topt', 'Tmax', 'temp_Tmin', 'temp_Topt', 'temp_Tmax'])
temps_compare.loc[:,'temp_Tmin':'temp_Tmax'] = tempura_otu_df.loc[temps_compare.index,'Tmin':'Tmax'].values

from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score

def regression_plot(x, y):
    X = x.values.reshape(-1, 1)
    Y = y.values.reshape(-1, 1)
    linear_regressor = LinearRegression()  # create object for the class
    linear_regressor.fit(X, Y)  # perform linear regression
    Y_pred = linear_regressor.predict(X)  # make predictions
    plt.scatter(X, Y)
    plt.plot(X, Y_pred, color='grey')
    mse = mean_squared_error(y_true=Y, y_pred=Y_pred, squared=True)
    rmse = mean_squared_error(y_true=Y, y_pred=Y_pred, squared=False)
    r_sq = r2_score(Y, Y_pred)
    plt.figtext(0.15,0.8, "RMSE: "+str("{:.3f}".format(rmse))+\
        "\n"+"{}\u00b2".format("R")+": "+str("{:.3f}".format(r_sq)))
    plt.xlabel("predicted " + x.name)
    plt.ylabel("TEMPURA " + x.name)
    plt.show()

regression_plot(temps_compare.Topt, temps_compare.temp_Topt)
regression_plot(temps_compare.Tmin, temps_compare.temp_Tmin)
regression_plot(temps_compare.Tmax, temps_compare.temp_Tmax)



