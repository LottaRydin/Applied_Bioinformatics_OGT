import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score

##load data
os.chdir(os.path.dirname(os.path.realpath(__file__))) #set current dir to script location

#script can generate these files from line 100 (tables which show how regresion coeff 
#changes using different filtering values), can take 10 min to generate them so can load if you already have it
#opt_opt = pd.read_csv("opt_opt.tsv", sep='\t')
#opt_opt = opt_opt.drop(columns=['Unnamed: 0']) #drop unnecessary col
#max_opt = pd.read_csv("max_opt.tsv", sep='\t')
#max_opt = max_opt.drop(columns=['Unnamed: 0']) #drop unnecessary col
#min_opt = pd.read_csv("min_opt.tsv", sep='\t')
#min_opt = min_opt.drop(columns=['Unnamed: 0']) #drop unnecessary col

tempura_otu_df = pd.read_csv("matched_temp_multi.tsv", sep='\t')
tempura_otu_df = tempura_otu_df.set_index(["OTU_id"])
pipeline_map = pd.read_csv('both_pipes.tsv', sep='\t')

df = pd.read_csv("data_whole.tsv", sep='\t') #compile_data_medium.py output
#filter out eukaryotes
df = df[~df.taxonomy.str.contains("Eukaryota")]
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

###############################################################################
#plot regression for Topt, Tmin or Tmax with data filtered using different parameters
#type_T: "opt", "min", "max", xy_transform if you want the regression to use raw MGnify values (False, default)
#or regression where MGnify values are transformed so regression intercept is 0 and slope is 1 (True)
def regression_filter_plot(type_T, min_read, min_samp, ratio, xy_transform = False):
            #Filter for number of samples per OTU
            df_filt = df[df.abundance > min_read] #drop sample with less that 3 reads
            #filter out OTUs with less than 10 observations
            OTU_size = df_filt.groupby(level="OTU_ID").size()
            df_filt = df_filt.drop(OTU_size.index[OTU_size < min_samp].tolist())
            tempura_matches = list(set(df_filt.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id")))
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
            col_dict = {"opt":[4, 1], "min":[3, 0], "max":[5, 2]}
            if xy_transform == False:
                a = 1; b = 0
            else:
                X = temps_pred.iloc[:,col_dict[type_T][0]].values.reshape(-1, 1)
                Y = temps_pred.iloc[:,col_dict[type_T][1]].values.reshape(-1, 1)
                linear_regressor = LinearRegression()  # create object for the class
                linear_regressor.fit(X, Y)  # perform linear regression
                Y_pred = linear_regressor.predict(X)  # make predictions
                a = linear_regressor.coef_[0][0]; b = linear_regressor.intercept_[0]
            X = temps_pred.iloc[:,col_dict[type_T][0]].values.reshape(-1, 1)
            Y = ((temps_pred.iloc[:,col_dict[type_T][1]]-b)/a).values.reshape(-1, 1)
            linear_regressor = LinearRegression()  # create object for the class
            linear_regressor.fit(X, Y)  # perform linear regression
            Y_pred = linear_regressor.predict(X)  # make predictions
            plt.scatter(X, Y, )
            plt.plot(X, Y_pred, color='grey')
            plt.plot()
            plt.plot([min([min(X), min(Y)]), max(max(X), max(Y))],\
                     [min([min(X), min(Y)]),max(max(X), max(Y))], color = 'black', linestyle='dashed')
            rmse = mean_squared_error(y_true=Y, y_pred=Y_pred, squared=False)
            r_sq = r2_score(Y, Y_pred)
            plt.figtext(0.28,0.73,"RMSE: "+str("{:.3f}".format(rmse))+\
                "\n"+"{}\u00b2".format("R")+": "+str("{:.3f}".format(r_sq))+"\n"+\
                    "y-intercept: "+ str("{:.3f}".format(linear_regressor.intercept_[0]))+"\n"+\
                        "slope: " + str("{:.3f}".format(linear_regressor.coef_[0][0])), size = "small")
            ax = plt.gca()
            ax.set_aspect('equal')
            plt.title("T"+type_T)
            plt.ylabel("prediction")
            plt.xlabel("TEMPURA")            
            plt.show()

regression_filter_plot("opt", 5, 5, 1, xy_transform=False) #opt
regression_filter_plot("opt", 5, 5, 1, xy_transform=True) #opt
regression_filter_plot("min", 7, 10, .9) #min
regression_filter_plot("max", 2, 16, .9, xy_transform=True) #max

###############################################################################
##generate tables which show how regression statistics change when filtering parameter are varied
#only need to run once on the data and save tables, they take about 10 min to generate for 1400 samples
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
             len(df_filt.index.unique("OTU_ID")), len(list(set(df_filt.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id"))))])

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
            tempura_matches = list(set(df_filt.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id")))
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
            tempura_matches = list(set(df_filt.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id")))
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