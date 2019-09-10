# -*- coding: utf-8 -*-

"""Look up the SNP from each pathway in MyVariant."""

import json
import os

import click

from snppis.api import batched_query
from snppis.constants import DATABASES, RESOURCES


@click.command()
def main():
    """Generate scores for each pathway."""
    for db in DATABASES:
        lookup_predictions(db)


def load_pathway_to_snps(db):
    with open(os.path.join(RESOURCES, f'{db}.json')) as file:
        pathway_to_snp = json.load(file)

    return {
        (
            entry['pathway']['namespace'],
            entry['pathway']['identifier'],
        ): [
            snp['identifier']
            for snp in entry['snps']
        ]
        for entry in pathway_to_snp
    }


def lookup_predictions(db):
    """Look up the SNP from each pathway in MyVariant."""
    pathway_to_snps = load_pathway_to_snps(db)

    dbsnp_ids = {
        snp
        for pathway, snps in pathway_to_snps.items()
        for snp in snps
    }

    click.echo(f'Got {len(dbsnp_ids)} SNPs from {db}')

    results = list(batched_query(dbsnp_ids))
    path = os.path.join(RESOURCES, f'{db}_snp_scores.json')
    with open(path, 'w') as file:
        json.dump(results, file, indent=2)


if __name__ == '__main__':
    main()
