from numpy.core.fromnumeric import mean
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from scipy.stats import zscore
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_val_score
import numpy as np
import math
import sklearn.metrics
import warnings
warnings.filterwarnings('error')

tempura_otu_df = pd.read_csv("matched_temp_multi.tsv", sep='\t')
tempura_otu_df = tempura_otu_df.set_index(["OTU_id_5"])

df = pd.read_csv("data_whole_combined.tsv", sep='\t') #compile_data_medium.py output

# Parameters
p_ratio = 0.55
p_min_samp = 5
p_min_read = 5
p_z_max = 2
p_z_min = 0.1
min_val = 10

df = df[df.temp < 110]
df = df[df.temp > -21]
df = df[~df.taxonomy.str.contains("Eukaryota")]

df = df.set_index(["OTU_ID"]) #set OTU col as the index
df.sort_index(inplace=True) #sort dataframe based on increasing index values
print("Out of", len(df.index.unique("OTU_ID")),"OTUs present in unfiltered MGnify data", \
      len(list(set(df.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id_5")))), "were present in TEMPURA" )

#tempura_matches = list(set(df_filt.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id_5")))
#cols = ['OTU_ID', 'Tmin', 'Topt', 'Tmax']

#Filter for number of samples per OTU
df_filt = df[df.abundance > p_min_read] #drop sample with less that 3 reads
#filter out OTUs with less than 10 observations
OTU_size = df_filt.groupby(level="OTU_ID").size()
df_filt = df_filt.drop(OTU_size.index[OTU_size < p_min_samp].tolist())

# Cutoff for z score
#z_cutoff = max_z

 # Removing outliers
tempura_matches = list(set(df_filt.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id_5")))
tempura_matches_1 = list(set(df_filt.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id_5")))
cols = ['OTU_ID', 'Tmin', 'Topt', 'Tmax']
temps_pred = []
#print(tempura_matches)
for otu in tempura_matches_1:
    # remove outliers
    otu_data = df_filt.loc[otu, :]
    otu_temp = df_filt.loc[otu, :].temp
    z_scores = zscore(otu_temp)
    abs_z_scores = np.abs(z_scores, dtype=object)
    filtered_entries_max = (abs_z_scores < p_z_max) #.all(axis=1)
    filtered_entries_min = (abs_z_scores < p_z_min) #.all(axis=1)
    max_df = otu_data[filtered_entries_max]
    min_df = otu_data[filtered_entries_min]
    #print(f_df.temp)

    if len(max_df.temp) > 0 and len(min_df.temp) > 0:
        min_est = min(min_df.temp) - min_val
        optimum_est = min_est + (max(max_df.temp) - min(min_df.temp)) * p_ratio
        temps_est = [otu, min_est, optimum_est, max(max_df.temp)]
        temps_pred.append(temps_est)
    else:
        #print(f'OTU failed: {otu}')
        if otu in tempura_matches:
            tempura_matches.remove(otu) 
    
    
temps_pred = pd.DataFrame(temps_pred, columns=cols) 
temps_pred = temps_pred.set_index(["OTU_ID"]) #set OTU col as the index
temps_pred = temps_pred.reindex(columns=['Tmin', 'Topt', 'Tmax', 'TEMP_Tmin', 'TEMP_Topt', 'TEMP_Tmax'])
temps_pred.loc[:,'TEMP_Tmin':'TEMP_Tmax'] = tempura_otu_df.loc[temps_pred.index,'Tmin':'Tmax'].values
col_dict = {"opt":[4, 1], "min":[3, 0], "max":[5, 2]}
#print(temps_pred)

X= temps_pred["Topt"].values.reshape(-1,1)
Y= temps_pred["TEMP_Topt"].values.reshape(-1,1)
lr = LinearRegression()
lr.fit(X,Y)

opt_pred = lr.predict(X)
opt_mse = mean_squared_error(Y,opt_pred)
opt_rmse = np.sqrt(opt_mse)
print(opt_rmse)

cvl = LeaveOneOut()
scores = cross_val_score(lr, X, Y, scoring="neg_mean_squared_error", cv=cvl)
lr_rmse_opt_scores = np.sqrt(-scores)

def display_scores(scores):
    print("Scores:", scores)
    print("Mean:", scores.mean())
    print("Standard deviation:", scores.std())

display_scores(lr_rmse_opt_scores)

    