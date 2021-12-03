import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import random
from collections import Counter
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score

##load data
df = pd.read_csv("/crex/proj/snic2021-23-617/work/02_compile_data/engineered/data_whole.tsv", sep='\t') #compile_data_medium.py output
df = df.set_index(["OTU_ID"]) #set OTU col as the index
metadata = pd.read_csv('/crex/proj/snic2021-23-617/work/01_convert_biom_to_tsv/engineered/temp_samples.tsv', sep='\t') #specify metadata file
metadata = metadata.drop(columns=['Unnamed: 0']) #drop unnecessary col

tempura_otu_df = pd.read_csv("/crex/proj/snic2021-23-617/raw_data/tempura/matched_temp_multi.tsv", sep='\t')
tempura_otu_df = tempura_otu_df.set_index(["OTU_id_5"])
tempura = pd.read_csv("/crex/proj/snic2021-23-617/raw_data/tempura/TEMPURA.csv") #load entire tempura db

##set filters
read_nr = 2
otu_size = 10
ratio = 0.5

df = df[df.temp < 110]
df = df[df.temp > -21]
df = df[~df.taxonomy.str.contains("Eukaryota")]

metadata = metadata[metadata.temperature < 110]
metadata = metadata[metadata.temperature > -21]

#Filter for number of samples per OTU
df = df[df.abundance > read_nr] #drop sample with less that 3 reads
#filter out OTUs with less than 10 observations
OTU_size = df.groupby(level="OTU_ID").size()
df = df.drop(OTU_size.index[OTU_size < otu_size].tolist())

##############################################################################
#seaborn violin plot
violin_df = []
rand_otu = random.sample(df.index.unique("OTU_ID").tolist(), 10) #pick 10 random otu ids
for indx in rand_otu:
    df_otu = df.loc[indx, :]
    violin_df.append(df_otu)

violin_df = pd.concat(violin_df)
violin_df.reset_index(inplace=True)

ax = sns.violinplot(x="taxonomy", y="temp", data=violin_df, scale="count", inner="stick") #If count, the width of the violins will be scaled by the number of observations in that bin. 
ax.tick_params(axis='x', rotation=90)
fig = ax.get_figure()
fig.savefig('violin_plot.png')

##############################################################################
#scatter plots to compare otu data to tempura
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
    plt.savefig(str(OTU_nr)+"_otu.png")
    plt.show()

print("Out of", len(df.index.unique("OTU_ID")),"OTUs present in MGnify data", \
      len(list(set(df.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id_5")))), "were present in TEMPURA" )

rand_otu = random.sample(list(set(df.index.unique("OTU_ID")) & \
                              set(tempura_otu_df.index.unique("OTU_id_5"))), 4) #pick random ots which have data also in tempura
for otu in rand_otu:
    tempura_and_prediction(otu)
 
################################################################################ 
#plot ratio of how close Topt is to Tmax in TEMPURA (but not all TEMPURA data, we use the file with 4.1 otu matches)
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
plt.savefig("TEMPURA_ratio_hist.png")
plt.show()
#tempura_ratio[tempura_ratio < 0] #how many ratios below 0, so that Topt is below Tmax
#tempura_ratio[(tempura_ratio > 0) & (tempura_ratio < 0.1)] #Topt very close to Tmin
###############################################################################
#some summary histograms
df.groupby(level="OTU_ID").size().plot(kind='hist', bins=30, rwidth=0.9, color='#014d4e', \
              title = "Number of samples per OTU" ) #df.groupby(level="OTU_ID")
plt.xlabel("sample nr")
plt.savefig("sample_per_otu.png")
plt.show()

metadata.temperature.plot(kind='hist', bins=20, rwidth=0.9, color='#ae7181', \
                          title = "Temperatures of samples" )
plt.xlabel("temperature")
plt.savefig("temperatures.png")
plt.show()

df.temp.plot(kind='hist', bins=20, rwidth=0.9, color='#607c8e', \
             title = "Temperatures of all observations" )
plt.xlabel("temperature")
plt.savefig("temperature_samples.png")
plt.show()

###############################################################################
#combining MGnify and TEMPURA temp data
cols = ['OTU_ID', 'Tmin', 'Topt', 'Tmax']
temps_pred = []
for otu in df.index.unique("OTU_ID"):
    optimum_est = min(df.loc[otu, :].temp) + (max(df.loc[otu, :].temp) - min(df.loc[otu, :].temp)) * ratio
    temps_est = [otu, min(df.loc[otu, :].temp), optimum_est, max(df.loc[otu, :].temp)]
    temps_pred.append(temps_est)

temps_pred = pd.DataFrame(temps_pred, columns=cols)
del cols
del temps_est
temps_pred = temps_pred.set_index(["OTU_ID"]) #set OTU col as the index
temps_pred.sort_index(inplace=True) #sort dataframe based on increasing index values

tempura_matches = list(set(df.index.unique("OTU_ID")) & set(tempura_otu_df.index.unique("OTU_id_5")))
temps_compare = temps_pred.loc[tempura_matches]
temps_compare = temps_compare.reindex(columns=['Tmin', 'Topt', 'Tmax', 'temp_Tmin', 'temp_Topt', 'temp_Tmax'])
temps_compare.loc[:,'temp_Tmin':'temp_Tmax'] = tempura_otu_df.loc[temps_compare.index,'Tmin':'Tmax'].values

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
    plt.figtext(0.15,0.73, "RMSE: "+str("{:.3f}".format(rmse))+\
        "\n"+"{}\u00b2".format("R")+": "+str("{:.3f}".format(r_sq))+"\n"+\
            "y-intercept: "+ str("{:.3f}".format(linear_regressor.intercept_[0]))+"\n"+\
                "slope: " + str("{:.3f}".format(linear_regressor.coef_[0][0])), size = "small", fontweight='bold')
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.savefig(title+".png")
    plt.show()

regression_plot(temps_compare.temp_Topt, temps_compare.Topt, "TEMPURA Topt", "predicted Topt", "Topt_prediction",)
regression_plot(temps_compare.temp_Tmin, temps_compare.Tmin, "TEMPURA Tmin", "predicted Tmin", "Tmin_prediction")
regression_plot(temps_compare.temp_Tmax, temps_compare.Tmax, "TEMPURA Tmax", "predicted Tmax", "Tmax_prediction")

regression_plot(tempura.Tmin, tempura.Tmax, "TEMPURA Tmin", "TEMPURA Tmax", "TEMPURA_tmax_tmin")
regression_plot(temps_pred.Tmin, temps_pred.Tmax, "prediction Tmin", "prediction Tmax", "Prediction_tmax_tmin")

regression_plot(tempura.Tmin, tempura.Topt_ave, "TEMPURA Tmin", "TEMPURA Topt", "TEMPURA_topt_tmin")
regression_plot(temps_pred.Tmin, temps_pred.Topt, "prediction Tmin", "prediction Topt", "Prediction_topt_tmin")

###############################################################################
#comparing temperture ranges
regression_plot(temps_compare.temp_Tmax-temps_compare.temp_Tmin,\
                temps_compare.Tmax-temps_compare.Tmin, \
                    "TEMPURA Tmax-Tmin", "prediction Tmax-Tmin", "Compare temperature ranges")





