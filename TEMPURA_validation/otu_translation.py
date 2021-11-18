
#import pandas for dataframes
import pandas as pd

fourone = pd.read_csv("ssu_4_1.otu", sep='\t', header=0) #read in the 4.1 OTU file

fourone

five = pd.read_csv("ssu_5.otu", sep='\t', header=0) #read in the 5.0 OTU file

four_test= fourone.copy() #make a copy to mess with

# does a left join of the two tables with taxonomy as key in order to match NCBI id from 5.0 with OTUs in 4.1
final_test=four_test.merge(five[['taxonomy', 'taxonomy_id']], left_on='taxonomy',right_on='taxonomy', how='left')

final_test

final_test.to_csv('final_otu_test.tsv', sep='\t', index=False) #write the table to a .tsv file

final_test['taxonomy_id'] = final_test['taxonomy_id'].astype(str).replace('\.0', '', regex=True) #removes the decimal points from the IDs

final_test

final_rounded = final_test.copy() #makes a copy to round

final_rounded

final_rounded.to_csv('final_otu_rounded.tsv', sep='\t', index=False) #write the version without decimal points to a .tsv file