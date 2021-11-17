import argparse
import csv
import os
import pandas

try:
    from urllib import urlencode
except ImportError:
    from urllib.parse import urlencode

import requests
from jsonapi_client import Session, Filter

# https://www.ebi.ac.uk/metagenomics/api/v1/samples?accession=&experiment_type=&biome_name=&lineage=root%3AEnvironmental%3AAquatic%3AMarine%3AIntertidal+zone%3AEstuary&geo_loc_name=&latitude_gte=&latitude_lte=&longitude_gte=&longitude_lte=&species=&instrument_model=&instrument_platform=&metadata_key=temperature&metadata_value_gte=&metadata_value_lte=&metadata_value=&environment_material=&environment_feature=&study_accession=&include=
# For testing: root:enviromental:Aquatic:Marine:Intertial Zone:Estuary (8 st)
# for real: root:enviromental:Aquatic:Marine:Intertial Zone:salt marsh (104 st)

API_FILTERED = "/metagenomics/api/v1/samples?accession=&experiment_type=&biome_name=&lineage=root%3AEnvironmental%3AAquatic%3AMarine%3AIntertidal+zone%3AEstuary&geo_loc_name=&latitude_gte=&latitude_lte=&longitude_gte=&longitude_lte=&species=&instrument_model=&instrument_platform=&metadata_key=temperature&metadata_value_gte=&metadata_value_lte=&metadata_value=&environment_material=&environment_feature=&study_accession=&include="
API_BASE = "https://www.ebi.ac.uk/metagenomics/api/v1"

# Make folder for tsv files. Should be changed to correct folder (Not in git repository)
if not os.path.exists("tsv_files"):
    os.mkdir("tsv_files")


# Making dataframe for storing temperature data:
temp_df = pandas.DataFrame(columns=['sample_accession', 'temperature', 'enviroment', 'sequencing_platform', 'pipeline_version', 'analysis_accession'])

# Extracting infor from MGnify:
with Session(API_BASE) as session:
    count = 0
    params = {
        'metadata_key': 'temperature', # Filter for temperature data
    }

    f = Filter(urlencode(params)) # Create filter

    # Iterating through samples from salt marshes (enviroment) with temperature data availible
    for sample in session.iterate('biomes/root:Environmental:Aquatic:Marine:Coastal/samples', f): #biomes/root:Environmental:Aquatic:Marine:Intermediate Zone/samples
        if count == 1000:
            break
        count += 1
        #print("sample:", sample.accession, ", count: ", count)
    # if True:
    #     sample = session.get('samples', "SRS373032").resource    

        # Finding the best analysis for this sample to use:
        best_run = None
        best_run_val = 0
        best_analysis = None
        pipeline_v = None
        instrument_platform = None

        # If keys are missing for analysis, but pipeline and experiment type still are good:
        backup_analysis = None
        backup_instrument_platform = None
        backup_pipeline_v = None
        backup_exists = False

        for run in sample.runs:
            #print("run acc: ", run.accession, "pipeline: ", [p.release_version for p in run.pipelines], "experiment type: ", run.experiment_type)
            
            pipelines_used = [p.release_version for p in run.pipelines]
 
            if ('4.1' in pipelines_used or '5.0' in pipelines_used) and run.experiment_type == 'amplicon':

                # Backup if keys are wrong. Done once per sample and takes first analysis.
                if run.analyses[0].pipeline_version in ['4.1', '5.0'] and not backup_exists:
                    backup_analysis = run.analyses[0].accession
                    backup_instrument_platform = run.instrument_platform
                    backup_pipeline_v = run.analyses[0].pipeline_version
                    backup_exists = True

                for analysis in run.analyses:
                    #print(f'ananlysis acc: {analysis.accession}, analysis pipeline: {analysis.pipeline_version}') #, analysis sum: {analysis.analysis_summary}')
                    
                    for e in analysis.analysis_summary:
                        if e['key'] == 'Reads with predicted RNA' or e['key'] == 'Predicted SSU sequences': 
                            if analysis.pipeline_version == '4.1':
                                if int(e['value']) >= best_run_val:
                                    best_run = run.accession
                                    best_run_val = int(e['value'])
                                    best_analysis = analysis.accession
                                    pipeline_v = '4.1'
                                    instrument_platform = run.instrument_platform
                            elif analysis.pipeline_version == '5.0':
                                if int(e['value']) >= best_run_val and (pipeline_v == None or pipeline_v == '5.0'):
                                    best_run = run.accession
                                    best_run_val = int(e['value'])
                                    best_analysis = analysis.accession
                                    pipeline_v = '5.0'
                                    instrument_platform = run.instrument_platform
            #else:
               # print("------------------NOT A GOOD PIPELINE/EXP TYPE---------------------")
            
        #print(f"best run: {best_run}, {best_run_val} in analysis: {best_analysis}")

        # Download TSV file for the best analysis:
        if best_analysis == None and not backup_exists:
            #print("-------------------no analysis was usable---------------------")
            continue
        else:
            RESOURCE = "analyses"
            if best_analysis != None:
                ANALYSIS_ACCESSION = best_analysis 
            else:
                ANALYSIS_ACCESSION = backup_analysis
                print('******BACKUP USED*******')

            analysis = session.get(RESOURCE, ANALYSIS_ACCESSION).resource

            # this will iterate over the download files from the API for the analysis
            # and will download only the otu tsv file.
            #print(API_BASE + "/" + RESOURCE + "/" + ANALYSIS_ACCESSION + "/downloads")
            for download in session.iterate(RESOURCE + "/" + ANALYSIS_ACCESSION + "/downloads"):
                # print(f'download alias: {download.alias}')

                if "FASTQ_SSU_OTU.tsv" in download.alias: # Only download tsv file
                    output_file = os.path.join("tsv_files", sample.accession + ".tsv")
                    
                    # if os.path.exists(output_file):
                    #     #print("Already downloaded: " + download.alias)
                    #     continue
                    #print("Downloading: " + download.alias)

                    with requests.get(download.links.self) as response:
                        response.raise_for_status()
                        with open(output_file, "wb") as f:
                            f.write(response.content)
                    break # Break loop, file is found
            
            # Add sample metadata to data frame after tsv file is downloaded

            env = sample.biome.id #sample.nvironment_material
            tmp = 'foo'
            for meta in sample.sample_metadata:
                if meta['key'] == 'temperature':
                    tmp = meta['value']
                # elif meta['key'] == 'environment (material)':
                #     env = meta['value'] 
            temp_df = temp_df.append({'sample_accession': sample.accession, 'temperature': tmp, 'enviroment': env, 'sequencing_platform': instrument_platform, 'pipeline_version': pipeline_v, 'analysis_accession': ANALYSIS_ACCESSION }, ignore_index=True)
temp_df.to_csv('temp_samples.tsv', sep='\t')

            


