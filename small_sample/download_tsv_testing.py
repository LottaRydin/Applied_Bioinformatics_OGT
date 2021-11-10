#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import csv
import os

import requests
from jsonapi_client import Session

API_BASE = "https://www.ebi.ac.uk/metagenomics/api/v1"
RESOURCE = "analyses"
ANALYSIS_ACCESSION = "MGYA00031511"

with Session(API_BASE) as session:

    analysis = session.get(RESOURCE, ANALYSIS_ACCESSION).resource

    print(analysis.accession)

    # Make folder for file. Should be changed for correct folder
    if not os.path.exists(ANALYSIS_ACCESSION):
        os.mkdir(ANALYSIS_ACCESSION)

    # this will iterate over the download files from the API for the analysis
    # and will download only the otu tsv file.
    print(API_BASE + "/" + RESOURCE + "/" + ANALYSIS_ACCESSION + "/downloads")
    for download in session.iterate(RESOURCE + "/" + ANALYSIS_ACCESSION + "/downloads"):
        print(f'download alias: {download.alias}')

        if "_MERGED_FASTQ_otu.tsv" in download.alias:
            output_file = os.path.join(ANALYSIS_ACCESSION, download.alias)
            if os.path.exists(output_file):
                print("Already downloaded: " + download.alias)
                continue
            print("Downloading: " + download.alias)
            with requests.get(download.links.self) as response:
                response.raise_for_status()
                with open(output_file, "wb") as f:
                    f.write(response.content)

    # # we can also download the functional annotation that is
    # # stored in the API
    # with open(os.path.join(GENOME_ACCESSION, "API_COG.csv"), "w") as cog_file:
    #     print("Getting the COG annotations from the API")
    #     fields = [
    #         "name",
    #         "description",
    #         "genome-count",
    #         "pangenome-count",
    #     ]
    #     writer = csv.DictWriter(cog_file, fieldnames=fields)
    #     writer.writeheader()
    #     for cog in session.iterate(RESOURCE + "/" + GENOME_ACCESSION + "/cogs"):
    #         writer.writerow(
    #             {
    #                 "name": cog.name,
    #                 "description": cog.description,
    #                 "genome-count": cog.genome_count,
    #                 "pangenome-count": cog.pangenome_count,
    #             }
    #         )
