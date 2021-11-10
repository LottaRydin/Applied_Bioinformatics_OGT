import argparse
import csv
import os

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

with Session(API_BASE) as session:
    params = {
        'metadata_key': 'temperature',
    }

    f = Filter(urlencode(params))
    for sample in session.iterate('biomes/root:Environmental:Aquatic:Marine:Intertidal zone:Estuary/samples', f):
        print("sample:", sample.accession)