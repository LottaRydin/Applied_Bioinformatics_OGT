from pandas import DataFrame

try:
    from urllib import urlencode
except ImportError:
    from urllib.parse import urlencode

from jsonapi_client import Session, Filter

API_BASE = 'https://www.ebi.ac.uk/metagenomics/api/latest/'

#https://www.ebi.ac.uk/metagenomics/api/v1/samples/SRS875056
with Session(API_BASE) as s:
    sample = s.get('samples', "SRS373032").resource
    print('Sample name:', sample.sample_name)
    print('sample enviroment:', sample.environment_material)
    for meta in sample.sample_metadata:
        print('Biome:', meta["key"], meta["value"])
    
    for run in sample.runs:
        print("acc: ", run.accession, "pipeline: ", [p.release_version for p in run.pipelines])
        for analysis in run.analyses:
            print(f'ananlysis acc: {analysis.accession}')
            print(analysis.relationships.downloads.links)
            for download in analysis.relationships.downloads:
                if download.description["label"]=="OTUs, counts and taxonomic assignments":
                    print(download.description["description"])
