# -*- coding: utf-8 -*-

"""Look up the SNP from each pathway in MyVariant."""

import json
import os

import click

from snppis.api import batched_query

HERE = os.path.abspath(os.path.dirname(__file__))
RESOURCES = os.path.join(HERE, 'mappings')


@click.command()
def main():
    """Generate scores for each pathway."""
    for db in ('kegg', 'wikipathways', 'reactome'):
        lookup_predictions(db)


def lookup_predictions(db):
    """Look up the SNP from each pathway in MyVariant."""
    with open(os.path.join(RESOURCES, f'{db}.json')) as file:
        pathway_to_snp = json.load(file)

    pathway_to_snp_keyed = {
        (
            entry['pathway']['namespace'],
            entry['pathway']['identifier'],
        ): [
            snp['identifier']
            for snp in entry['snps']
        ]
        for entry in pathway_to_snp
    }

    dbsnp_ids = {
        snp
        for pathway, snps in pathway_to_snp_keyed.items()
        for snp in snps
    }

    click.echo(f'Got {len(dbsnp_ids)} SNPs from {db}')

    results = list(batched_query(dbsnp_ids))
    path = os.path.join(RESOURCES, f'{db}_snp_scores.json')
    with open(path, 'w') as file:
        json.dump(results, file, indent=2)


if __name__ == '__main__':
    main()
