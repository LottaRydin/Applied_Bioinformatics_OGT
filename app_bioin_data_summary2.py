import os
import pandas as pd

os.chdir("C:/Users/Gretka/Documents/MA_Bioinformatics/Applied_Bioinformatics/data_summary/") #metadada directory
data_dir = "small_data_salt_marsh/" #tsv file directory
metadata = pd.read_csv('temp_samples_salt_marsh.tsv', sep='\t') #modify metadata.tsv name

appended_df = []
for filename in os.listdir(data_dir):
    sample_acc = filename.split(".", 1)[0] #get sample acc from file name
    df = pd.read_csv(data_dir + filename, sep='\t') #load df for an analysis
    runID = df.columns[1] #save run ID just for reference
    df = df.rename(columns={runID: "abundance", "OTU ID": "OTU_ID"}) #rename OTU, abundance col
    df = df[~df.taxonomy.str.contains("sk__Eukaryota")] #remove OTU classified as eukaryotes
    df['rel_abundance'] = df.iloc[:, 1] / sum(df.iloc[:, 1]) #add relative abundance col
    df['run_ID'] = runID #add col with run ID
    df['temp'] = metadata.loc[metadata['sample_accession'] == sample_acc, 'temperature'].values[0] #add temp col
    df = df.reindex(columns=["OTU_ID", "run_ID", "temp", "abundance", "rel_abundance", "taxonomy"]) #reorder cols for future hierarchical indexing
    appended_df.append(df)
    
appended_df = pd.concat(appended_df)
appended_df = appended_df.set_index(["OTU_ID"]) #set OTU col as the index
appended_df.sort_index(inplace=True) #sort dataframe based on increasing index values

for indx in appended_df.index.unique("OTU_ID")[:10]: #just for the first 10 OTUs, grab data
    print(appended_df.loc[indx, :])


