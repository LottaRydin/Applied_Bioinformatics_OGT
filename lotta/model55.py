from numpy.core.fromnumeric import mean
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from scipy.stats import zscore
import numpy as np
import math
import sklearn.metrics
import warnings
warnings.filterwarnings('error')

tempura_otu_df = pd.read_csv("matched_temp_multi.tsv", sep='\t')
tempura_otu_df = tempura_otu_df.set_index(["OTU_id_5"])

# Load data that will be used for predicition
df_for_pred = pd.read_csv("data_whole_combined.tsv", sep='\t') # Data for prediction of growth tempertures
df_result = pd.DataFrame(columns=[ 'opt', 'min', 'max', 'taxonomy']) # Dataframe for predicted values


# Load training data
df = pd.read_csv("data_whole_combined.tsv", sep='\t') # Data for training model

# Parameters
ratio = 0.55
min_samp = 5
min_read = 5
max_z = 0.5
min_z = 0.5

# ----------------- Preprocess training data ------------------

# Filter for temperature and organisms
df = df[df.temp < 110]
df = df[df.temp > -21]
df = df[~df.taxonomy.str.contains("Eukaryota")]

df = df.set_index(["OTU_ID"]) #set OTU col as the index
df.sort_index(inplace=True) #sort dataframe based on increasing index values

#Filter for number of samples per OTU
df_filt = df[df.abundance > min_read] #drop sample with less that 3 reads

#filter out OTUs with less than 10 observations
OTU_size = df_filt.groupby(level="OTU_ID").size()
df_filt = df_filt.drop(OTU_size.index[OTU_size < min_samp].tolist())

# Find matches in tempura for training
tempura_matches = list(set(df_filt.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id_5")))
tempura_matches_1 = list(set(df_filt.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id_5")))
cols = ['OTU_ID', 'Tmin', 'Topt', 'Tmax']
temps_pred_train = []

# Make initial predicitions for training data
for otu in tempura_matches_1:
    # remove outliers
    otu_data = df_filt.loc[otu, :]
    otu_temp = df_filt.loc[otu, :].temp
    z_scores = zscore(otu_temp)
    abs_z_scores = np.abs(z_scores, dtype=object)
    filtered_entries_max = (abs_z_scores < max_z) #.all(axis=1)
    filtered_entries_min = (abs_z_scores < min_z) #.all(axis=1)
    max_df = otu_data[filtered_entries_max]
    min_df = otu_data[filtered_entries_min]

    if len(max_df.temp) > 0 and len(min_df.temp) > 0:
        # Make initial temperature estimations
        optimum_est = min(min_df.temp) + (max(max_df.temp) - min(min_df.temp)) * ratio
        max_est = max(max_df.temp)
        temps_est = [otu, min(min_df.temp), optimum_est, max_est]
        temps_pred_train.append(temps_est)
    else:
        if otu in tempura_matches:
            tempura_matches.remove(otu)

# Reform predictions into dataframe
temps_pred_train = pd.DataFrame(temps_pred_train, columns=cols) 
temps_pred_train = temps_pred_train.set_index(["OTU_ID"]) #set OTU col as the index
temps_pred_train = temps_pred_train.reindex(columns=['Tmin', 'Topt', 'Tmax', 'TEMP_Tmin', 'TEMP_Topt', 'TEMP_Tmax'])
temps_pred_train.loc[:,'TEMP_Tmin':'TEMP_Tmax'] = tempura_otu_df.loc[temps_pred_train.index,'Tmin':'Tmax'].values
col_dict = {"opt":[4, 1], "min":[3, 0], "max":[5, 2]}

# ----------- preprocessing data for predicting ---------------

# Filter for temperature and organisms
df_for_pred = df_for_pred[df_for_pred.temp < 110]
df_for_pred = df_for_pred[df_for_pred.temp > -21]
df_for_pred = df_for_pred[~df_for_pred.taxonomy.str.contains("Eukaryota")]

df_for_pred = df_for_pred.set_index(["OTU_ID"]) #set OTU col as the index
df_for_pred.sort_index(inplace=True) #sort dataframe based on increasing index values

#Filter for number of samples per OTU
df_for_pred_filt = df_for_pred[df_for_pred.abundance > min_read] #drop sample with less that 3 reads

#filter out OTUs with less than 5 observations
OTU_size = df_for_pred_filt.groupby(level="OTU_ID").size()
df_for_pred_filt = df_for_pred_filt.drop(OTU_size.index[OTU_size < min_samp].tolist())

# Make initial predicitions for prediction data
temps_pred = []
otu_list = list(df_for_pred_filt.index.unique("OTU_ID"))
otu_list_1 = list(df_for_pred_filt.index.unique("OTU_ID"))

cols_pred = ['OTU_ID', 'Tmin', 'Topt', 'Tmax', 'taxonomy']
for otu in otu_list_1:
    # remove outliers
    otu_data = df_filt.loc[otu, :]
    otu_temp = df_filt.loc[otu, :].temp
    z_scores = zscore(otu_temp)
    abs_z_scores = np.abs(z_scores, dtype=object)
    filtered_entries_max = (abs_z_scores < max_z) #.all(axis=1)
    filtered_entries_min = (abs_z_scores < min_z) #.all(axis=1)
    max_df = otu_data[filtered_entries_max]
    min_df = otu_data[filtered_entries_min]

    if len(max_df.temp) > 0 and len(min_df.temp) > 0:
        # Make initial temperature estimations
        optimum_est = min(min_df.temp) + (max(max_df.temp) - min(min_df.temp)) * ratio
        max_est = max(max_df.temp)
        temps_est = [otu, min(min_df.temp), optimum_est, max_est, min(min_df.taxonomy)]
        temps_pred.append(temps_est)
    else:
        if otu in tempura_matches:
            otu_list.remove(otu)

# Reform predictions into dataframe
temps_pred = pd.DataFrame(temps_pred, columns=cols_pred) 
temps_pred = temps_pred.set_index(["OTU_ID"]) #set OTU col as the index

for type_T in ["opt", "min", "max"]:

    # Training and evaluating model for all subsets
    train_set = temps_pred_train

    # Training model: find information about regression for tranformation
    X = temps_pred_train.iloc[:,col_dict[type_T][0]].values.reshape(-1, 1)
    Y = temps_pred_train.iloc[:,col_dict[type_T][1]].values.reshape(-1, 1)
    linear_regressor = LinearRegression()  # create object for the class
    linear_regressor.fit(X, Y)  # perform linear regression
    a = linear_regressor.coef_[0][0]; b = linear_regressor.intercept_[0]

    # Make prediction
    preds = ((temps_pred.iloc[:,col_dict[type_T][1]]-b)/a) #.values.reshape(-1, 1)
    df_result[type_T]=preds

df_result['taxonomy']=temps_pred.iloc[:,3]
print(df_result)

df_result.to_csv('predicterd_temperatures.tsv', sep='\t')
print(f'tempura matches used for training: {len(tempura_matches)}')


