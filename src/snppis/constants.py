import os

HERE = os.path.abspath(os.path.dirname(__file__))
RESOURCES = os.path.join(HERE, 'mappings')

DATABASES = [
    'kegg',
    'wikipathways',
    'reactome',
]
