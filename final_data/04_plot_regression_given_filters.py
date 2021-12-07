import pandas as pd
import matplotlib.pyplot as plt
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

def regression_filter_plot(type_T, min_read, min_samp, ratio, xy_transform = False):
            #Filter for number of samples per OTU
            df_filt = df[df.abundance > min_read] #drop sample with less that 3 reads
            #filter out OTUs with less than 10 observations
            OTU_size = df_filt.groupby(level="OTU_ID").size()
            df_filt = df_filt.drop(OTU_size.index[OTU_size < min_samp].tolist())
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
            plt.plot([min([min(X), min(Y)])-15, max([max(X), max(Y)])+10],\
                     [min([min(X), min(Y)])-15,max([max(X), max(Y)])+10], color = 'black', linestyle='dashed')
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
            plt.savefig("T"+type_T+".png")
            plt.show()
            print(len(list(set(df_filt.index.unique("OTU_ID")) &\
                                       set(tempura_otu_df.index.unique("OTU_id_5")))),"OTUs plotted",  "for", type_T, "regression")
#            print("OTUs plotted", len(list(set(df_filt.index.unique("OTU_ID")))))
                                       
regression_filter_plot("opt", 6, 25, .35, xy_transform=False) #opt
regression_filter_plot("min", 6, 25, .35) #min
regression_filter_plot("max", 6, 25, .35) #max

#if you want the transformed plot, set xy_transform=True (default in False)
#regression_filter_plot("opt", 6, 25, .35, xy_transform=True) #opt

