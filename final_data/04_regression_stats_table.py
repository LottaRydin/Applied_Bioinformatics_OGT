import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score

tempura_otu_df = pd.read_csv("matched_temp_multi.tsv", sep='\t')
tempura_otu_df = tempura_otu_df.set_index(["OTU_id_5"])

df = pd.read_csv("data_whole_combined.tsv", sep='\t') #compile_data_medium.py output

df = df[df.temp < 110]
df = df[df.temp > -21]
df = df[~df.taxonomy.str.contains("Eukaryota")]

df = df.set_index(["OTU_ID"]) #set OTU col as the index
df.sort_index(inplace=True) #sort dataframe based on increasing index values
print("Out of", len(df.index.unique("OTU_ID")),"OTUs present in unfiltered MGnify data", \
      len(list(set(df.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id_5")))), "were present in TEMPURA" )

#function which returns regression stats given x and y inputs
def regression_stats(x, y):
    X = x.values.reshape(-1, 1)
    Y = y.values.reshape(-1, 1)
    linear_regressor = LinearRegression()  # create object for the class
    linear_regressor.fit(X, Y)  # perform linear regression
    Y_pred = linear_regressor.predict(X)  # make predictions
    rmse = mean_squared_error(y_true=Y, y_pred=Y_pred, squared=False)
    r_sq = r2_score(Y, Y_pred)
    return([r_sq, rmse, linear_regressor.intercept_[0],linear_regressor.coef_[0][0], \
             len(df_filt.index.unique("OTU_ID")), len(list(set(df_filt.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id_5"))))])

#two big loops below generate 3 data frames whih describe how regression stats change 
#when you filter data using different parameters. After you generate them you can
#save them as tsv files to use later as they take a long time to generate
max_opt = []
min_opt = []
for min_read in list(range(1, 11)):
    for min_sampl_nr in list(range(2, 50)):
            #Filter for number of samples per OTU
            df_filt = df[df.abundance > min_read] #drop sample with less that 3 reads
            #filter out OTUs with less than 10 observations
            OTU_size = df_filt.groupby(level="OTU_ID").size()
            df_filt = df_filt.drop(OTU_size.index[OTU_size < min_sampl_nr].tolist())
            tempura_matches = list(set(df_filt.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id_5")))
            cols = ['OTU_ID', 'Tmin', 'Topt', 'Tmax']
            temps_pred = []
            for otu in tempura_matches:
                optimum_est = min(df_filt.loc[otu, :].temp) + (max(df_filt.loc[otu, :].temp) - min(df_filt.loc[otu, :].temp)) * 0.6
                temps_est = [otu, min(df_filt.loc[otu, :].temp), optimum_est, max(df_filt.loc[otu, :].temp)]
                temps_pred.append(temps_est)
            temps_pred = pd.DataFrame(temps_pred, columns=cols)
            temps_pred = temps_pred.set_index(["OTU_ID"]) #set OTU col as the index
            temps_pred = temps_pred.reindex(columns=['Tmin', 'Topt', 'Tmax', 'TEMP_Tmin', 'TEMP_Topt', 'TEMP_Tmax'])
            temps_pred.loc[:,'TEMP_Tmin':'TEMP_Tmax'] = tempura_otu_df.loc[temps_pred.index,'Tmin':'Tmax'].values
            loop_values_max = [min_read, min_sampl_nr]
            loop_values_min = [min_read, min_sampl_nr]
            loop_values_max.extend(regression_stats(temps_pred.TEMP_Tmax, temps_pred.Tmax))
            loop_values_min.extend(regression_stats(temps_pred.TEMP_Tmin, temps_pred.Tmin))
            max_opt.append(loop_values_max)
            min_opt.append(loop_values_min)
cols = ['min_reads', 'min_sampl_nr', 'r_sq', 'rmse', "y_intercept", "slope", "OTU_nr", "TEMPURA_matches"]
max_opt = pd.DataFrame(max_opt, columns=cols)
min_opt = pd.DataFrame(min_opt, columns=cols)

opt_opt = []
for min_read in list(range(1, 11)):
    for min_sampl_nr in list(range(2, 50)):
        for ratio in list(np.arange(0, 1.05, 0.05)):
            #Filter for number of samples per OTU
            df_filt = df[df.abundance > min_read] #drop sample with less that 3 reads
            #filter out OTUs with less than 10 observations
            OTU_size = df_filt.groupby(level="OTU_ID").size()
            df_filt = df_filt.drop(OTU_size.index[OTU_size < min_sampl_nr].tolist())
            tempura_matches = list(set(df_filt.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id_5")))
            cols = ['OTU_ID', 'Tmin', 'Topt', 'Tmax']
            temps_pred = []
            for otu in tempura_matches:
                optimum_est = min(df_filt.loc[otu, :].temp) + (max(df_filt.loc[otu, :].temp) - min(df_filt.loc[otu, :].temp)) * ratio
                temps_est = [otu, min(df_filt.loc[otu, :].temp), optimum_est, max(df_filt.loc[otu, :].temp)]
                temps_pred.append(temps_est)
            temps_pred = pd.DataFrame(temps_pred, columns=cols)
            temps_pred = temps_pred.set_index(["OTU_ID"]) #set OTU col as the index
            temps_pred = temps_pred.reindex(columns=['Tmin', 'Topt', 'Tmax', 'TEMP_Tmin', 'TEMP_Topt', 'TEMP_Tmax'])
            temps_pred.loc[:,'TEMP_Tmin':'TEMP_Tmax'] = tempura_otu_df.loc[temps_pred.index,'Tmin':'Tmax'].values
            loop_values = [min_read, min_sampl_nr, ratio]
            loop_values.extend(regression_stats(temps_pred.TEMP_Topt, temps_pred.Topt))
            opt_opt.append(loop_values)
cols = ['min_reads', 'min_sampl_nr', 'ratio','r_sq', 'rmse', "y_intercept", "slope", "OTU_nr", "TEMPURA_matches"]
opt_opt = pd.DataFrame(opt_opt, columns=cols)

max_opt.to_csv("max_opt.tsv", sep='\t')
min_opt.to_csv("min_opt.tsv", sep='\t')
opt_opt.to_csv("opt_opt.tsv", sep='\t')

