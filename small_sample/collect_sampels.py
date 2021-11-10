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
temp_df = pandas.DataFrame(columns=['sample_accession', 'temperature'])

# Extracting infor from MGnify:
with Session(API_BASE) as session:
    count = 0
    params = {
        'metadata_key': 'temperature', # Filter for temperature data
    }

    f = Filter(urlencode(params)) # Create filter

    # Iterating through samples from salt marshes (enviroment) with temperature data availible
    for sample in session.iterate('biomes/root:Environmental:Aquatic:Marine:Intertidal zone:Salt marsh/samples', f):
        count += 1
        print("sample:", sample.accession, ", count: ", count)
    # if True:
    #     sample = session.get('samples', "SRS373032").resource    

        # Finding the best analysis for this sample to use:
        best_run = None
        best_run_val = 0
        best_analysis = None

        for run in sample.runs:
            print("run acc: ", run.accession, "pipeline: ", [p.release_version for p in run.pipelines], "experiment type: ", run.experiment_type)
            pipelines_used = [p.release_version for p in run.pipelines]
            if '4.1' in pipelines_used and run.experiment_type == 'amplicon':
                for analysis in run.analyses:
                    print(f'ananlysis acc: {analysis.accession}, analysis pipeline: {analysis.pipeline_version}') #, analysis sum: {analysis.analysis_summary}')
                    for e in analysis.analysis_summary:
                        if e['key'] == 'Reads with predicted RNA' and analysis.pipeline_version == '4.1':
                            if int(e['value']) >= best_run_val:
                                best_run = run.accession
                                best_run_val = int(e['value'])
                                best_analysis = analysis.accession
            else:
                print("------------------NO GOOD PIPELINE/EXP TYPE---------------------")
            
        print(f"best run: {best_run}, {best_run_val} in analysis: {best_analysis}")

        # Download TSV file for the best analysis:
        if best_analysis == None:
            print("-------------------no analysis was usable---------------------")
        else:
            RESOURCE = "analyses"
            ANALYSIS_ACCESSION = best_analysis 

            analysis = session.get(RESOURCE, ANALYSIS_ACCESSION).resource

            # this will iterate over the download files from the API for the analysis
            # and will download only the otu tsv file.
            print(API_BASE + "/" + RESOURCE + "/" + ANALYSIS_ACCESSION + "/downloads")
            for download in session.iterate(RESOURCE + "/" + ANALYSIS_ACCESSION + "/downloads"):
                # print(f'download alias: {download.alias}')

                if "FASTQ_SSU_OTU.tsv" in download.alias: # Only download tsv file
                    output_file = os.path.join("tsv_files", sample.accession + "_" + download.alias + ".tsv")
                    if os.path.exists(output_file):
                        print("Already downloaded: " + download.alias)
                        continue
                    print("Downloading: " + download.alias)
                    with requests.get(download.links.self) as response:
                        response.raise_for_status()
                        with open(output_file, "wb") as f:
                            f.write(response.content)
            
            # Add sample temperature to data frame after tsv file is downloaded
            for meta in sample.sample_metadata:
                if meta['key'] == 'temperature':
                    temp_df = temp_df.append({'sample_accession': sample.accession, 'temperature': meta['value']}, ignore_index=True)
                    break
temp_df.to_csv('temp_samples.tsv', sep='\t')

            


