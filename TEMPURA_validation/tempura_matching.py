import pandas as pd #import pandas for dataframes

#Read in sample OTUs
otus = pd.read_csv("ssu_4_1.otu", sep='\t')

otus

#Read in 5.0 pipeline ref database OTUs
otus_5 = pd.read_csv("ssu_5.otu", sep='\t')

otus_5

#Converts taxonomy_id from float to sting and remove decimal point
otus_5['taxonomy_id'] = otus_5['taxonomy_id'].astype(str).replace('\.0', '', regex=True) #removes the decimal points from the IDs

#Rename OTU column
otus_5.rename(columns={'OTU_id': 'OTU_id_5'}, inplace=True)

otus_5

#Read in pipeline 4.1 OTUs with taxonomy_id added
translation = pd.read_csv("final_otu_rounded.tsv", sep= '\t')

#Converts taxonomy_id from float to sting and remove decimal point
translation['taxonomy_id'] = translation['taxonomy_id'].astype(str).replace('\.0', '', regex=True) #removes the decimal points from the IDs

translation

#Rename OTU column
translation.rename(columns={'OTU_id': 'OTU_id_4_1'}, inplace=True)

translation

#Joins sample OTU table with the 4.1 ref database on OTU_id
translated=otus.merge(translation[['OTU_id_4_1', 'taxonomy_id']], left_on='OTU_id',right_on='OTU_id_4_1', how='left')

translated

#For multi pipeline use translated table is merged with pipeline 5.0 ref database on taxonomy hierarchy
translated_multi=translated.merge(otus_5[['OTU_id_5','taxonomy']], left_on='taxonomy',right_on='taxonomy')

translated_multi

#Changes order of columns
translated_multi_v2 = translated_multi[['OTU_id','taxonomy','OTU_id_4_1','OTU_id_5','taxonomy_id']]

translated_multi_v2

#Reads in the TEMPURA database
tempura = pd.read_csv("TEMPURA.csv")

tempura

#Changes column data type to string
ttempura = tempura.astype({"taxonomy_id": str})

ttempura

#Merges the single-pipeline 4.1 OTU-table with matching TEMPURA data
matched_temp =translated.merge(ttempura[['taxonomy_id', 'strain','Tmin', 'Topt_ave', 'Tmax','Topt_low','Topt_high','Tmax_Tmin']], left_on='taxonomy_id',right_on='taxonomy_id')

matched_temp

#Writes the matched table to a .tsv file
matched_temp.to_csv('matched_temp.tsv', sep='\t', index=False)

#Merges the multi-pipeline OTU-table with matching TEMPURA data
matched_temp_multi =translated_multi_v2.merge(ttempura[['taxonomy_id', 'strain','Tmin', 'Topt_ave', 'Tmax','Topt_low','Topt_high','Tmax_Tmin']], left_on='taxonomy_id',right_on='taxonomy_id')

matched_temp_multi

#Writes the matched table to a .tsv file
matched_temp_multi.to_csv('matched_temp_multi.tsv', sep='\t', index=False)