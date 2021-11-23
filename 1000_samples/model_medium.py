import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import random
from collections import Counter
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score

##load data
os.chdir(os.path.dirname(os.path.realpath(__file__))) #set current dir to script location
tempura_otu_df = pd.read_csv("matched_temp_multi.tsv", sep='\t')
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
##fitering
#filter out eukaryotes
df = df[~df.taxonomy.str.contains("Eukaryota")]

#Filter for number of samples per OTU
df = df[df.abundance > 1] #drop sample with less that 3 reads
#filter out OTUs with less than 10 observations
OTU_size = df.groupby(level="OTU_ID").size()
df = df.drop(OTU_size.index[OTU_size < 10].tolist())
del list_1
del list_2
del otu_dict
del missing_v5_id
del OTU_size

cols = ['OTU_ID', 'Tmin', 'Topt', 'Tmax']
temps_pred = []
for otu in df.index.unique("OTU_ID"):
    optimum_est = min(df.loc[otu, :].temp) + (max(df.loc[otu, :].temp) - min(df.loc[otu, :].temp)) *0.66
    temps_est = [otu, min(df.loc[otu, :].temp), optimum_est, max(df.loc[otu, :].temp)]
    temps_pred.append(temps_est)

temps_pred = pd.DataFrame(temps_pred, columns=cols)
del cols
del temps_est
del otu
del optimum_est
temps_pred = temps_pred.set_index(["OTU_ID"]) #set OTU col as the index
temps_pred.sort_index(inplace=True) #sort dataframe based on increasing index values

tempura_matches = list(set(df.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id")))
temps_compare = temps_pred.loc[tempura_matches]
temps_compare = temps_compare.reindex(columns=['Tmin', 'Topt', 'Tmax', 'temp_Tmin', 'temp_Topt', 'temp_Tmax'])
temps_compare.loc[:,'temp_Tmin':'temp_Tmax'] = tempura_otu_df.loc[temps_compare.index,'Tmin':'Tmax'].values
del tempura_matches
###############################################################################


def regression_plot(x, y, xlab, ylab, title):
    X = x.values.reshape(-1, 1)
    Y = y.values.reshape(-1, 1)
    linear_regressor = LinearRegression()  # create object for the class
    linear_regressor.fit(X, Y)  # perform linear regression
    Y_pred = linear_regressor.predict(X)  # make predictions
    plt.scatter(X, Y)
    plt.plot(X, Y_pred, color='grey')
#    mse = mean_squared_error(y_true=Y, y_pred=Y_pred, squared=True)
    rmse = mean_squared_error(y_true=Y, y_pred=Y_pred, squared=False)
    r_sq = r2_score(Y, Y_pred)
    plt.figtext(0.15,0.8, "RMSE: "+str("{:.3f}".format(rmse))+\
        "\n"+"{}\u00b2".format("R")+": "+str("{:.3f}".format(r_sq)))
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.show()
    print("y-intercept "+ str(linear_regressor.intercept_))
    print("slope " + str(linear_regressor.coef_))

regression_plot(temps_compare.temp_Topt, temps_compare.Topt, "TEMPURA Topt", "predicted Topt", "Topt prediction")
regression_plot(temps_compare.temp_Tmin, temps_compare.Tmin, "TEMPURA Tmin", "predicted Tmin", "Tmin prediction")
regression_plot(temps_compare.temp_Tmax, temps_compare.Tmax, "TEMPURA Tmax", "predicted Tmax", "Tmax prediction")

#https://stats.stackexchange.com/questions/214886/regression-what-is-the-utility-of-r-squared-compared-to-rmse




