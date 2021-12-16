#0 as min, halfway between 0 and max
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score

from numpy.core.fromnumeric import mean
from scipy.stats import zscore
import math
import sklearn.metrics
import warnings
warnings.filterwarnings('error')

from scipy import stats
##load data

tempura_otu_df = pd.read_csv("matched_temp_multi.tsv", sep='\t')
tempura_otu_df = tempura_otu_df.set_index(["OTU_id_5"])

df = pd.read_csv("data_whole_combined.tsv", sep='\t') #compile_data_medium.py output

minimum = 10
min_reads = 5
min_samples = 5
ratio = .55
min_z = 0.1
max_z = 2
df = df[df.temp < 110]
df = df[df.temp > -21]
df = df[~df.taxonomy.str.contains("Eukaryota")]

df = df.set_index(["OTU_ID"]) #set OTU col as the index
df.sort_index(inplace=True) #sort dataframe based on increasing index values


print("Out of", len(df.index.unique("OTU_ID")),"OTUs present in unfiltered MGnify data", \
      len(list(set(df.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id_5")))), "were present in TEMPURA" )

def basic_regression_filter_plot(type_T, min_read, min_samp, ratio, xy_transform = False):
            #Filter for number of samples per OTU
            df_filt = df[df.abundance > min_read] #drop sample with less that 3 reads
            #filter out OTUs with less than 10 observations
            OTU_size = df_filt.groupby(level="OTU_ID").size()
            df_filt = df_filt.drop(OTU_size.index[OTU_size < min_samp].tolist())
            df_filt = df_filt[df_filt.groupby(['OTU_ID'])['temp'].transform(stats.zscore).abs() < 2]
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

116740
145929
def regression_filter_plot(type_T, min_read, min_samp, ratio, min_z, max_z, xy_transform = False):
            #Filter for number of samples per OTU
            df_filt = df[df.abundance > min_read] #drop sample with less that 3 reads
            #filter out OTUs with less than 10 observations
            OTU_size = df_filt.groupby(level="OTU_ID").size()
            df_filt = df_filt.drop(OTU_size.index[OTU_size < min_samp].tolist())
            cols = ['OTU_ID', 'Tmin', 'Topt', 'Tmax']
            temps_pred = []
#            min_z = 2
#            max_z = 2
            min_df = df_filt[df_filt.groupby(['OTU_ID'])['temp'].transform(stats.zscore).abs() < min_z]
            max_df = df_filt[df_filt.groupby(['OTU_ID'])['temp'].transform(stats.zscore).abs() < max_z]
            tempura_matches = list(set(df_filt.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id_5")))

            otu_list = list(set(min_df.index.unique("OTU_ID")) & set(max_df.index.unique("OTU_ID")))
            otu_plot = list(set(otu_list) & set(tempura_matches))
            for otu in otu_plot:
                
                if type(min_df.loc[otu, :].temp) == np.float64:
                    minimum_est = min_df.loc[otu, :].temp
                else:
                    minimum_est = min(min_df.loc[otu, :].temp)
                    
                if type(max_df.loc[otu, :].temp) == np.float64:
                    maximum_est = max_df.loc[otu, :].temp
                else:
                    maximum_est = max(max_df.loc[otu, :].temp)
                    
                optimum_est = minimum_est + (maximum_est - minimum_est) * ratio
                temps_est = [otu, minimum_est, optimum_est, maximum_est]
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
            plt.figtext(0.28,0.68,"RMSE: "+str("{:.3f}".format(rmse))+\
                "\n"+"{}\u00b2".format("R")+": "+str("{:.3f}".format(r_sq))+"\n"+\
                    "y-intercept: "+ str("{:.3f}".format(linear_regressor.intercept_[0]))+"\n"+\
                        "slope: " + str("{:.3f}".format(linear_regressor.coef_[0][0]))\
                        +"\n"+"min_z: "+ str(min_z)+"\n"+"max_z: "+ str(max_z), size = "small")
            ax = plt.gca()
            ax.set_aspect('equal')
            plt.title("T"+type_T)
            plt.ylabel("prediction")
            plt.xlabel("TEMPURA")  
            plt.savefig("T"+type_T+".png")
            plt.show()
#            print(len(list(set(df_filt.index.unique("OTU_ID")) &\
#                                       set(tempura_otu_df.index.unique("OTU_id_5")))),"OTUs plotted",  "for", type_T, "regression")
            print("OTUs plotted", len(temps_pred))
basic_regression_filter_plot("max", min_reads, min_samples, ratio, xy_transform=False) #opt
regression_filter_plot("max", min_reads, min_samples, ratio, min_z, max_z, xy_transform=False) #opt

def leave_1_out(df, type_T, min_read, min_samp, ratio, max_z, xy_transform = False):

    #Filter for number of samples per OTU
    df_filt = df[df.abundance > min_read] #drop sample with less that 3 reads
    #filter out OTUs with less than 10 observations
    OTU_size = df_filt.groupby(level="OTU_ID").size()
    df_filt = df_filt.drop(OTU_size.index[OTU_size < min_samp].tolist())
    df_filt = df_filt[df_filt.groupby(['OTU_ID'])['temp'].transform(stats.zscore).abs() < max_z]
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
    #print(temps_pred)

    # Statistics for all iterations
    r_sq_all = []
    rmse_all = []
    true_all = []
    pred_all = []

    # Training and evaluating model for all subsets
    for otu in tempura_matches:
        #print(f'--------------OTU: {otu}-------------')
        #print(temps_pred)
        train_set = temps_pred.iloc[lambda x: x.index != otu]
        #print(df_filt.iloc[lambda x: x.index == otu])
        #print(temps_pred.iloc[lambda x: x.index == otu])
        #print(f'temp: {otu_temp}')

        if xy_transform == False:
            a = 1; b = 0
        else:
            # Find information abour regression for tranformation
            X = train_set.iloc[:,col_dict[type_T][0]].values.reshape(-1, 1)
            Y = train_set.iloc[:,col_dict[type_T][1]].values.reshape(-1, 1)
            linear_regressor = LinearRegression()  # create object for the class
            linear_regressor.fit(X, Y)  # perform linear regression
            Y_pred = linear_regressor.predict(X)  # make predictions
            a = linear_regressor.coef_[0][0]; b = linear_regressor.intercept_[0]

        # Values for test otu
        otu_TEMP = temps_pred.iloc[lambda x: x.index == otu,col_dict[type_T][0]].values.reshape(-1, 1)
        otu_pred = ((temps_pred.iloc[lambda x: x.index == otu,col_dict[type_T][1]]-b)/a).values.reshape(-1, 1)
        otu_mmo = temps_pred.iloc[lambda x: x.index == otu,col_dict[type_T][1]].values.reshape(-1, 1)
        
        # preform linear regression
        X = train_set.iloc[:,col_dict[type_T][0]].values.reshape(-1, 1)
        Y = ((train_set.iloc[:,col_dict[type_T][1]]-b)/a).values.reshape(-1, 1)
        linear_regressor = LinearRegression()  # create object for the class
        linear_regressor.fit(X, Y)  # perform linear regression
        #Y_pred_otu = linear_regressor.predict(otu_TEMP)
        Y_pred = linear_regressor.predict(X) # make predictions
        #rmse = mean_squared_error(y_true=otu_TEMP, y_pred=otu_pred, squared=False)
        # print(f'OTU: {otu}')
        # print(otu_TEMP)
        # print(otu_pred)
        rmse = float(np.abs(otu_TEMP-otu_pred, dtype=object))
        true_all.append(otu_TEMP[0])
        pred_all.append(otu_pred[0])
        rmse_all.append(rmse)

        plt.scatter(otu_TEMP, otu_pred )
        plt.annotate(round(rmse, 1), (otu_TEMP, otu_pred ))
        plt.plot(X, Y_pred, color='grey')
        #print(f'OTU: {otu}, min/max/opt: {otu_mmo} prediction: {otu_pred}, TEMPURA: {otu_TEMP}')

        #plot.plot()
        X = [x[0] for x in X]
        Y = [x[0] for x in Y]
        # print(X)
        # print(type(Y))
        plt.plot([min([min(X), min(Y)]), 120],\
                    [min([min(X), min(Y)]),120], color = 'black', linestyle='dashed')

        #r_sq = r2_score(otu_TEMP, otu_pred)
        #abs = np.abs(otu_TEMP-otu_pred)
        #rmse_lr = mean_squared_error(y_true=Y, y_pred=Y_pred_otu, squared=False)
        #r_sq_lr = r2_score(Y, Y_pred_otu)
    # print(true_all)
    # print(pred_all)
    mse = sklearn.metrics.mean_squared_error(true_all, pred_all)
    r2 = sklearn.metrics.r2_score(true_all, pred_all)
    true_rmse = rmse = math.sqrt(mse)
    plt.figtext(0.28,0.73,"mean error: "+str("{:.3f}".format(mean(rmse_all)))+"\n"+\
        "rmse: "+str("{:.3f}".format(mean(true_rmse)))+"\n"+\
        "\n"+"{}\u00b2".format("R")+": "+str("{:.3f}".format(r2))+"\n"+\
        #f"Abs diff: {abs}"+"\n"+\
        "y-intercept: "+ str("{:.3f}".format(linear_regressor.intercept_[0]))+"\n"+\
        "slope: " + str("{:.3f}".format(linear_regressor.coef_[0][0]))+"\n"+\
        f"Z cutoff: {max_z}", size = "small")
        # "RMSE_f: "+str("{:.3f}".format(rmse_lr))+\
        # "\n"+"{}\u00b2".format("R_f")+": "+str("{:.3f}".format(r_sq_lr))+"\n"+\
        # "y-intercept_f: "+ str("{:.3f}".format(linear_regressor.intercept_[0]))+"\n"+\
        # "slope_f: " + str("{:.3f}".format(linear_regressor.coef_[0][0])), size = "small")
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.title("T"+type_T)
    plt.ylabel("prediction")
    plt.xlabel("TEMPURA")
    plt.show()

    #print(rmse_all)
    print(f'mean rmse: {mean(rmse_all)}')
    plt.hist(rmse_all, bins=40)
    plt.show()

leave_1_out(df, "opt", min_reads, min_samples, ratio, 2, xy_transform=True) #opt
#leave_1_out(df, "max", 5, 5, .55, 2, xy_transform=True) #opt
#leave_1_out(df, "min", 5, 5, .55, 2, xy_transform=True) #opt

###############################################################################
###############################################################################

def leave_1_out(df, type_T, min_read, min_samp, ratio, min_z, max_z, xy_transform = False):

    #Filter for number of samples per OTU
    df_filt = df[df.abundance > min_read] #drop sample with less that 3 reads
    #filter out OTUs with less than 10 observations
    OTU_size = df_filt.groupby(level="OTU_ID").size()
    df_filt = df_filt.drop(OTU_size.index[OTU_size < min_samp].tolist())
    cols = ['OTU_ID', 'Tmin', 'Topt', 'Tmax']
    temps_pred = []
    min_df = df_filt[df_filt.groupby(['OTU_ID'])['temp'].transform(stats.zscore).abs() < min_z]
    max_df = df_filt[df_filt.groupby(['OTU_ID'])['temp'].transform(stats.zscore).abs() < max_z]
    tempura_matches = list(set(df_filt.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id_5")))

    otu_list = list(set(min_df.index.unique("OTU_ID")) & set(max_df.index.unique("OTU_ID")))
    otu_plot = list(set(otu_list) & set(tempura_matches))
    for otu in otu_plot:
        
        if type(min_df.loc[otu, :].temp) == np.float64:
            minimum_est = min_df.loc[otu, :].temp
        else:
            minimum_est = min(min_df.loc[otu, :].temp)
            
        if type(max_df.loc[otu, :].temp) == np.float64:
            maximum_est = max_df.loc[otu, :].temp
        else:
            maximum_est = max(max_df.loc[otu, :].temp)
            
        optimum_est = minimum_est + (maximum_est - minimum_est) * ratio
        temps_est = [otu, minimum_est, optimum_est, maximum_est]
        temps_pred.append(temps_est)
    temps_pred = pd.DataFrame(temps_pred, columns=cols) 
    temps_pred = temps_pred.set_index(["OTU_ID"]) #set OTU col as the index
    temps_pred = temps_pred.reindex(columns=['Tmin', 'Topt', 'Tmax', 'TEMP_Tmin', 'TEMP_Topt', 'TEMP_Tmax'])
    temps_pred.loc[:,'TEMP_Tmin':'TEMP_Tmax'] = tempura_otu_df.loc[temps_pred.index,'Tmin':'Tmax'].values
    col_dict = {"opt":[4, 1], "min":[3, 0], "max":[5, 2]}
    #print(temps_pred)

    # Statistics for all iterations
    r_sq_all = []
    rmse_all = []
    true_all = []
    pred_all = []

    # Training and evaluating model for all subsets
    for otu in otu_plot:
        #print(f'--------------OTU: {otu}-------------')
        #print(temps_pred)
        train_set = temps_pred.iloc[lambda x: x.index != otu]
        #print(df_filt.iloc[lambda x: x.index == otu])
        #print(temps_pred.iloc[lambda x: x.index == otu])
        #print(f'temp: {otu_temp}')

        if xy_transform == False:
            a = 1; b = 0
        else:
            # Find information abour regression for tranformation
            X = train_set.iloc[:,col_dict[type_T][0]].values.reshape(-1, 1)
            Y = train_set.iloc[:,col_dict[type_T][1]].values.reshape(-1, 1)
            linear_regressor = LinearRegression()  # create object for the class
            linear_regressor.fit(X, Y)  # perform linear regression
            Y_pred = linear_regressor.predict(X)  # make predictions
            a = linear_regressor.coef_[0][0]; b = linear_regressor.intercept_[0]

        # Values for test otu
        otu_TEMP = temps_pred.iloc[lambda x: x.index == otu,col_dict[type_T][0]].values.reshape(-1, 1)
        otu_pred = ((temps_pred.iloc[lambda x: x.index == otu,col_dict[type_T][1]]-b)/a).values.reshape(-1, 1)
        otu_mmo = temps_pred.iloc[lambda x: x.index == otu,col_dict[type_T][1]].values.reshape(-1, 1)
        
        # preform linear regression
        X = train_set.iloc[:,col_dict[type_T][0]].values.reshape(-1, 1)
        Y = ((train_set.iloc[:,col_dict[type_T][1]]-b)/a).values.reshape(-1, 1)
        linear_regressor = LinearRegression()  # create object for the class
        linear_regressor.fit(X, Y)  # perform linear regression
        #Y_pred_otu = linear_regressor.predict(otu_TEMP)
        Y_pred = linear_regressor.predict(X) # make predictions
        #rmse = mean_squared_error(y_true=otu_TEMP, y_pred=otu_pred, squared=False)
        # print(f'OTU: {otu}')
        # print(otu_TEMP)
        # print(otu_pred)
        rmse = float(np.abs(otu_TEMP-otu_pred, dtype=object))
        true_all.append(otu_TEMP[0])
        pred_all.append(otu_pred[0])
        rmse_all.append(rmse)

        plt.scatter(otu_TEMP, otu_pred )
        plt.annotate(round(rmse, 1), (otu_TEMP, otu_pred ))
        plt.plot(X, Y_pred, color='grey')
        #print(f'OTU: {otu}, min/max/opt: {otu_mmo} prediction: {otu_pred}, TEMPURA: {otu_TEMP}')

        #plot.plot()
        X = [x[0] for x in X]
        Y = [x[0] for x in Y]
        # print(X)
        # print(type(Y))
        plt.plot([min([min(X), min(Y)]), 120],\
                    [min([min(X), min(Y)]),120], color = 'black', linestyle='dashed')

        #r_sq = r2_score(otu_TEMP, otu_pred)
        #abs = np.abs(otu_TEMP-otu_pred)
        #rmse_lr = mean_squared_error(y_true=Y, y_pred=Y_pred_otu, squared=False)
        #r_sq_lr = r2_score(Y, Y_pred_otu)
    # print(true_all)
    # print(pred_all)
    mse = sklearn.metrics.mean_squared_error(true_all, pred_all)
    r2 = sklearn.metrics.r2_score(true_all, pred_all)
    true_rmse = rmse = math.sqrt(mse)
    plt.figtext(0.28,0.73,"mean error: "+str("{:.3f}".format(mean(rmse_all)))+"\n"+\
        "rmse: "+str("{:.3f}".format(mean(true_rmse)))+"\n"+\
        "\n"+"{}\u00b2".format("R")+": "+str("{:.3f}".format(r2))+"\n"+\
        #f"Abs diff: {abs}"+"\n"+\
        "y-intercept: "+ str("{:.3f}".format(linear_regressor.intercept_[0]))+"\n"+\
        "slope: " + str("{:.3f}".format(linear_regressor.coef_[0][0]))+"\n"+\
        f"Z cutoff: {max_z}", size = "small")
        # "RMSE_f: "+str("{:.3f}".format(rmse_lr))+\
        # "\n"+"{}\u00b2".format("R_f")+": "+str("{:.3f}".format(r_sq_lr))+"\n"+\
        # "y-intercept_f: "+ str("{:.3f}".format(linear_regressor.intercept_[0]))+"\n"+\
        # "slope_f: " + str("{:.3f}".format(linear_regressor.coef_[0][0])), size = "small")
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.title("T"+type_T)
    plt.ylabel("prediction")
    plt.xlabel("TEMPURA")
    plt.show()

    #print(rmse_all)
    print(f'mean rmse: {mean(rmse_all)}')
    plt.hist(rmse_all, bins=40)
    plt.show()

leave_1_out(df, "opt", min_reads, min_samples, ratio, min_z, max_z, xy_transform=True) #opt
#leave_1_out(df, "max", 5, 5, .55, 2, xy_transform=True) #opt
#leave_1_out(df, "min", 5, 5, .55, 2, xy_transform=True) #opt


