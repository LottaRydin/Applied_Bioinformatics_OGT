import argparse
import csv
import os
import pandas
import time

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


# Making dataframes for storing temperature data and data from connection errors:
temp_df = pandas.DataFrame(columns=['sample_accession', 'temperature', 'enviroment', 'sequencing_platform', 'pipeline_version', 'analysis_accession'])
samples_error = pandas.DataFrame(columns=['sample_accession'])
analysis_error = pandas.DataFrame(columns=['analysis_accession'])
backup_record = pandas.DataFrame(columns=['sample_accession'])

# Extracting infor from MGnify:
with Session(API_BASE) as session:
    count = 0
    add_count = 0
    params = {
        'metadata_key': 'temperature', # Filter for temperature data
    }

    f = Filter(urlencode(params)) # Create filter

    # sample_list=session.iterate('biomes/root:Environmental:Aquatic:Marine/samples', f)
    # print('STARTING')
    # time.sleep(5)
    
    # Iterating through samples from salt marshes (enviroment) with temperature data availible
    try:
        for sample in session.iterate('biomes/root:Environmental/samples', f): #biomes/root:Environmental:Aquatic:Marine:Intermediate Zone/samples
        #for sample in sample_list:
            if count == 20000 or add_count == 20000:
                print(f'FINAL: Explored: {count}, downloaded: {add_count}')
                break
            elif add_count%50 == 0:
                print(f'Explored: {count}, downloaded: {add_count}')
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
            
            max_attempts = 3
            for attempt in range(1, max_attempts+1):
                try:   
                    for run in sample.runs:
                        #print("run acc: ", run.accession, "pipeline: ", [p.release_version for p in run.pipelines], "experiment type: ", run.experiment_type)
                        try:
                            pipelines_used = [p.release_version for p in run.pipelines]
                        except requests.exceptions.SSLError:
                            print('Error rad 68!')
            
                        if ('4.1' in pipelines_used or '5.0' in pipelines_used) and run.experiment_type == 'amplicon':

                            # Backup if keys are wrong. Done once per sample and takes first analysis.
                            if run.analyses[0].pipeline_version in ['4.1', '5.0'] and not backup_exists:
                                backup_analysis = run.analyses[0].accession
                                backup_instrument_platform = run.analyses[0].instrument_platform
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
                                                instrument_platform = analysis.instrument_platform
                                        elif analysis.pipeline_version == '5.0':
                                            if int(e['value']) >= best_run_val and (pipeline_v == None or pipeline_v == '5.0'):
                                                best_run = run.accession
                                                best_run_val = int(e['value'])
                                                best_analysis = analysis.accession
                                                pipeline_v = '5.0'
                                                instrument_platform = analysis.instrument_platform
                        #else:
                            #print("------------------NOT A GOOD PIPELINE/EXP TYPE---------------------")
                except (requests.exceptions.SSLError, requests.exceptions.ConnectionError):
                    print(f'Connection error when retrieving information about sample. Attempt: {attempt}, Sample: {sample.accession}')
                    #time.sleep(5)
                    if attempt == max_attempts:
                        print(f'Maximum attempts for retrieving information about sample reached. {sample.accession} added to excess file')
                        #add file
                        samples_error = samples_error.append({'sample_accession': sample.accession}, ignore_index=True)
                else:
                    break
                
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
                    backup_record = backup_record.append({'sample_accession': ANALYSIS_ACCESSION}, ignore_index=True)
                    print('******BACKUP USED*******')
                
                
                for attempt in range(1, max_attempts+1):
                    try:
                        #print("DOWNLOAD")
                        #time.sleep(5)
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
                                print("Downloading: " + download.alias)
                                
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
                        add_count += 1
                        #print("DOWNLOAD FINISHED")
                    except (requests.exceptions.SSLError, requests.exceptions.ConnectionError):
                        print(f'Connection error when dowloading TSV file for sample. Attempt: {attempt}, Analysis: {ANALYSIS_ACCESSION}')
                        #time.sleep(5)
                        if attempt == max_attempts:
                            print(f'Maximum attempts for dowloading TSV file for sample reached. {ANALYSIS_ACCESSION} added to excess file')
                            analysis_error = analysis_error.append({'analysis_accession': ANALYSIS_ACCESSION }, ignore_index=True)
                            best_analysis = None
                            backup_exists = False
                    else:
                        break
    except (requests.exceptions.SSLError, requests.exceptions.ConnectionError):
        print('Long lasting connection errors. Download interrupted.')
        temp_df.to_csv('temp_samples.tsv', sep='\t')
        samples_error.to_csv('samples_not_included', sep='\t')
        analysis_error.to_csv('analyses_not_downloaded', sep='\t')
        backup_record.to_csv('backups', sep='\t')
temp_df.to_csv('temp_samples.tsv', sep='\t')
samples_error.to_csv('samples_not_included', sep='\t')
analysis_error.to_csv('analyses_not_downloaded', sep='\t')
backup_record.to_csv('backups', sep='\t')
print(f'FINAL: Explored: {count}, downloaded: {add_count}')