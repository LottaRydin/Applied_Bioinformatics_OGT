import os
import pandas as pd

#the location of the script should also contain tsv_files/ folder with MGnify tsv
#files (files converted using biom_to_tsv.sh) and temp_samples.tsv metadata file

os.chdir(os.path.dirname(os.path.realpath(__file__))) #set current dir to script location
data_dir = "tsv_files/" #tsv sample files' directory
metadata = pd.read_csv('temp_samples.tsv', sep='\t') #specify metadata file
metadata = metadata.drop(columns=['Unnamed: 0']) #drop unnecessary col
      
#adding metadata to tsv files and storing them as df s in a list
appended_df = []
for filename in os.listdir(data_dir):
    sample_acc = filename.split(".", 1)[0] #get sample acc from the tsv file name
    df = pd.read_csv(data_dir + filename, sep='\t') #load a tsv as df
    df = df.rename(columns={df.columns[1]: "abundance", "OTU ID": "OTU_ID"}) #rename OTU, abundance col
    df['rel_abundance'] = df.iloc[:, 1] / sum(df.iloc[:, 1]) #add relative abundance col
    df['sample_acc'] = sample_acc #add col with run ID
    df['temp'] = metadata.loc[metadata['sample_accession'] == sample_acc, 'temperature'].values[0] #add metadata columns
    df['pipeline'] = metadata.loc[metadata['sample_accession'] == sample_acc, 'pipeline_version'].values[0]
    df['environment'] = metadata.loc[metadata['sample_accession'] == sample_acc, 'enviroment'].values[0]
    df['seq_platform'] = metadata.loc[metadata['sample_accession'] == sample_acc, 'sequencing_platform'].values[0]
    df = df.reindex(columns=["OTU_ID", "sample_acc", "temp", "abundance", "rel_abundance", \
                             "taxonomy", 'environment', 'seq_platform', 'pipeline']) #reorder cols for future hierarchical indexing
    df = df[df.taxonomy.str.contains("g__")] #keep only OTUs resolved to at least genus level
    appended_df.append(df) #add df to the appended_df list

#converting the list of df s into a single df
appended_df = pd.concat(appended_df)
appended_df = appended_df.set_index(["OTU_ID"]) #set OTU col as the index
appended_df.sort_index(inplace=True) #sort dataframe based on increasing index values

#write appended_df as tsv file
appended_df.to_csv("data_whole.tsv", sep='\t')

