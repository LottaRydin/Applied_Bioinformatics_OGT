import pandas as pd

environmental = pd.read_csv("/crex/proj/snic2021-23-617/work/02_compile_data/environmental/data_whole.tsv", sep='\t')
engineered = pd.read_csv("/crex/proj/snic2021-23-617/work/02_compile_data/engineered/data_whole.tsv", sep='\t')
host = pd.read_csv("/crex/proj/snic2021-23-617/work/02_compile_data/host-associated/data_whole.tsv", sep='\t')
mixed = pd.read_csv("/crex/proj/snic2021-23-617/work/02_compile_data/mixed/data_whole.tsv", sep='\t')

combined_whole = environmental.append([engineered, host, mixed])
combined_whole = combined_whole.set_index(["OTU_ID"])
combined_whole.sort_index(inplace=True) #sort dataframe based on increasing index values
combined_whole.to_csv("data_whole_combined.tsv", sep='\t')

environmental = pd.read_csv("/crex/proj/snic2021-23-617/work/01_convert_biom_to_tsv/environmental/temp_samples.tsv", sep='\t')
engineered = pd.read_csv("/crex/proj/snic2021-23-617/work/01_convert_biom_to_tsv/engineered/temp_samples.tsv", sep='\t')
host = pd.read_csv("/crex/proj/snic2021-23-617/work/01_convert_biom_to_tsv/host-associated/temp_samples.tsv", sep='\t')
mixed = pd.read_csv("/crex/proj/snic2021-23-617/work/01_convert_biom_to_tsv/mixed/temp_samples.tsv", sep='\t')

environmental = environmental.drop(columns=['Unnamed: 0']) 
engineered = engineered.drop(columns=['Unnamed: 0'])
host = host.drop(columns=['Unnamed: 0'])
mixed = mixed.drop(columns=['Unnamed: 0']) 

combined_whole = environmental.append([engineered, host, mixed])
combined_whole.to_csv("temp_samples_combined.tsv", sep='\t')


