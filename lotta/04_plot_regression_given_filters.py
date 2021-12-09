import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from scipy.stats import zscore
import numpy as np

tempura_otu_df = pd.read_csv("matched_temp_multi.tsv", sep='\t')
tempura_otu_df = tempura_otu_df.set_index(["OTU_id_5"])

df = pd.read_csv("data_whole_combined.tsv", sep='\t') #compile_data_medium.py output

# Parameters
p_ratio = 0.25
p_min_samp = 48
p_min_read = 7
p_z_cutoff = 1

df = df[df.temp < 110]
df = df[df.temp > -21]
df = df[~df.taxonomy.str.contains("Eukaryota")]

df = df.set_index(["OTU_ID"]) #set OTU col as the index
df.sort_index(inplace=True) #sort dataframe based on increasing index values
print("Out of", len(df.index.unique("OTU_ID")),"OTUs present in unfiltered MGnify data", \
      len(list(set(df.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id_5")))), "were present in TEMPURA" )

def regression_filter_plot(type_T, min_read, min_samp, ratio, max_z, xy_transform = False):
            #Filter for number of samples per OTU
            df_filt = df[df.abundance > min_read] #drop sample with less that 3 reads
            #filter out OTUs with less than 10 observations
            OTU_size = df_filt.groupby(level="OTU_ID").size()
            df_filt = df_filt.drop(OTU_size.index[OTU_size < min_samp].tolist())

            # Cutoff for z score
            z_cutoff = max_z

            # Not removing outliers
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
            col_dict = {"opt":[4, 1], "min":[3, 0], "max":[5, 2]}

            # Removing outliers
            tempura_matches_f = list(set(df_filt.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id_5")))
            cols_f = ['OTU_ID', 'Tmin', 'Topt', 'Tmax']
            temps_pred_f = []
            for otu in tempura_matches:
                # remove outliers
                otu_data = df_filt.loc[otu, :]
                otu_temp = df_filt.loc[otu, :].temp
                z_scores = zscore(otu_temp)
                abs_z_scores = np.abs(z_scores)
                filtered_entries = (abs_z_scores < z_cutoff) #.all(axis=1)
                f_df = otu_data[filtered_entries]
                #print(f_df.temp)

                if len(f_df.temp) > 0:
                    optimum_est = min(f_df.temp) + (max(f_df.temp) - min(f_df.temp)) * ratio
                    temps_est = [otu, min(f_df.temp), optimum_est, max(f_df.temp)]
                    temps_pred_f.append(temps_est)
                else:
                    print(f'OTU failed: {otu}')
                    if otu in tempura_matches_f:
                        tempura_matches_f.remove(otu) 
            temps_pred_f = pd.DataFrame(temps_pred_f, columns=cols_f) 
            temps_pred_f = temps_pred_f.set_index(["OTU_ID"]) #set OTU col as the index
            temps_pred_f = temps_pred_f.reindex(columns=['Tmin', 'Topt', 'Tmax', 'TEMP_Tmin', 'TEMP_Topt', 'TEMP_Tmax'])
            temps_pred_f.loc[:,'TEMP_Tmin':'TEMP_Tmax'] = tempura_otu_df.loc[temps_pred_f.index,'Tmin':'Tmax'].values
            col_dict_f = {"opt":[4, 1], "min":[3, 0], "max":[5, 2]}

            if xy_transform == False:
                a = 1; b = 0
                a_f = 1; b_f = 0
            else:
                # For outliers not removed
                X = temps_pred.iloc[:,col_dict[type_T][0]].values.reshape(-1, 1)
                Y = temps_pred.iloc[:,col_dict[type_T][1]].values.reshape(-1, 1)
                linear_regressor = LinearRegression()  # create object for the class
                linear_regressor.fit(X, Y)  # perform linear regression
                Y_pred = linear_regressor.predict(X)  # make predictions
                a = linear_regressor.coef_[0][0]; b = linear_regressor.intercept_[0]

                # For outliers removed
                X_f = temps_pred_f.iloc[:,col_dict_f[type_T][0]].values.reshape(-1, 1)
                Y_f = temps_pred_f.iloc[:,col_dict_f[type_T][1]].values.reshape(-1, 1)
                linear_regressor_f = LinearRegression()  # create object for the class
                linear_regressor_f.fit(X_f, Y_f)  # perform linear regression
                Y_pred_f = linear_regressor_f.predict(X_f)  # make predictions
                a_f = linear_regressor_f.coef_[0][0]; b_f = linear_regressor_f.intercept_[0]

            # For outliers not removed
            X = temps_pred.iloc[:,col_dict[type_T][0]].values.reshape(-1, 1)
            Y = ((temps_pred.iloc[:,col_dict[type_T][1]]-b)/a).values.reshape(-1, 1)
            linear_regressor = LinearRegression()  # create object for the class
            linear_regressor.fit(X, Y)  # perform linear regression
            Y_pred = linear_regressor.predict(X)  # make predictions
            plt.scatter(X, Y, )
            plt.plot(X, Y_pred, color='grey')

            # For outliers removed
            X_f = temps_pred_f.iloc[:,col_dict_f[type_T][0]].values.reshape(-1, 1)
            Y_f = ((temps_pred_f.iloc[:,col_dict_f[type_T][1]]-b_f)/a_f).values.reshape(-1, 1)
            linear_regressor_f = LinearRegression()  # create object for the class
            linear_regressor_f.fit(X_f, Y_f)  # perform linear regression
            Y_pred_f = linear_regressor_f.predict(X_f)  # make predictions
            plt.scatter(X_f, Y_f, c='darkorange')
            plt.plot(X_f, Y_pred_f, color='peru')

            plt.plot()
            plt.plot([min([min(X), min(Y)]), 120],\
                     [min([min(X), min(Y)]),120], color = 'black', linestyle='dashed')
            rmse = mean_squared_error(y_true=Y, y_pred=Y_pred, squared=False)
            r_sq = r2_score(Y, Y_pred)
            rmse_f = mean_squared_error(y_true=Y_f, y_pred=Y_pred_f, squared=False)
            r_sq_f = r2_score(Y_f, Y_pred_f)
            plt.figtext(0.28,0.73,"RMSE: "+str("{:.3f}".format(rmse))+\
                "\n"+"{}\u00b2".format("R")+": "+str("{:.3f}".format(r_sq))+"\n"+\
                "y-intercept: "+ str("{:.3f}".format(linear_regressor.intercept_[0]))+"\n"+\
                "slope: " + str("{:.3f}".format(linear_regressor.coef_[0][0]))+"\n"+\
                f"Z cutoff: {max_z}"+"\n"+\
                "RMSE_f: "+str("{:.3f}".format(rmse_f))+\
                "\n"+"{}\u00b2".format("R_f")+": "+str("{:.3f}".format(r_sq_f))+"\n"+\
                "y-intercept_f: "+ str("{:.3f}".format(linear_regressor_f.intercept_[0]))+"\n"+\
                "slope_f: " + str("{:.3f}".format(linear_regressor_f.coef_[0][0])), size = "small")
            ax = plt.gca()
            ax.set_aspect('equal')
            plt.title("T"+type_T)
            plt.ylabel("prediction")
            plt.xlabel("TEMPURA")
            if xy_transform:
                plt.savefig("T"+type_T+f"z{max_z}_trans"+".png")
            else:
                 plt.savefig("T"+type_T+f"z{max_z}"+".png")
            plt.show()
#            print("OTUs plotted", list(set(df_filt.index.unique("OTU_ID")) &\
#                                       set(tempura_otu_df.index.unique("OTU_id_5"))), "for", type_T, "regression")
            print("OTUs plotted", len(list(set(df_filt.index.unique("OTU_ID")))))
                                       
regression_filter_plot("opt", p_min_read, p_min_samp, p_ratio, p_z_cutoff, xy_transform=False) #opt
regression_filter_plot("opt", p_min_read, p_min_samp, p_ratio, p_z_cutoff, xy_transform=True) #opt
regression_filter_plot("min", p_min_read, p_min_samp, p_ratio, p_z_cutoff,) #min
regression_filter_plot("min", p_min_read, p_min_samp, p_ratio, p_z_cutoff, xy_transform=True) #min
regression_filter_plot("max", p_min_read, p_min_samp, p_ratio, p_z_cutoff,) #max
regression_filter_plot("max", p_min_read, p_min_samp, p_ratio, p_z_cutoff, xy_transform=True) #max


